
import numpy as np
import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from irspec.spec_helpers import trim_spec, cut_line, fit_continuum, find_nans
from astropy.io import ascii
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region
from matplotlib.ticker import MultipleLocator

from irspec.fitfuncs import * 
from lmfit.models import GaussianModel
from lmfit import Parameters
from astropy.table import Table
from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import astropy.constants as const
from matplotlib.colors import LogNorm, Normalize
from astropy.visualization.wcsaxes import add_scalebar
import matplotlib.font_manager as fm
import scipy.special as sp
from scipy import ndimage as ndi
from scipy.stats import f as f_dist

from irspec.datacube import Datacube
from irspec.plotparams import PlotParams
from irspec.emission_io import read_line_params
from irspec import paths

# Importing this module styles nothing globally; each Spaxelcube carries its
# own PlotParams (see `__init__`) and re-asserts it before every figure.
DEFAULT_PALETTE = PlotParams.DEFAULT_PALATTE  # "light"
DEFAULT_SCALING = "presentation"

# Imaging whose contours a map can be overlaid with, instead of the cube's
# own footprint outline. Paths are relative to the image-data root, which is
# resolved lazily through irspec.paths so the maps render the same from
# wherever they are launched.
# Keyed by filter, since that is what identifies the image scientifically;
# the tag a map's filename carries is the key.
CONTOUR_IMAGES = {
    "f200w": "nircam/IR23128_f200w.fits",
    "f356w": "nircam/IR23128_f356w.fits",
    "f560w": "miri/jw03368-o047_t031_miri_f560w-brightsky_i2d.fits",
    "f770w": "miri/jw03368-o047_t031_miri_f770w-brightsky_i2d.fits",
    "f1500w": "miri/jw03368-o047_t031_miri_f1500w-brightsky_i2d.fits",
}

# Log-spaced contour levels, and the percentile of the in-frame pixels the
# lowest one sits at. Starting well up the distribution keeps the faint end
# -- noise at this depth -- from burying the map under contours that trace
# the read pattern rather than the galaxy.
CONTOUR_LEVELS = 8
CONTOUR_FLOOR_PERCENTILE = 92.0


def log_contour_levels(values, n_levels=CONTOUR_LEVELS,
                       floor_percentile=CONTOUR_FLOOR_PERCENTILE):
    """Log-spaced contour levels spanning a percentile floor to the maximum.

    `values` should already be restricted to the finite, positive pixels of
    the region being drawn: a log spacing is undefined at or below zero, and
    levels set from the whole mosaic would not suit the small patch shown.
    Returns an empty array when there is nothing to contour.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.array([])
    low = np.percentile(values, floor_percentile)
    high = values.max()
    if not high > low:
        low, high = values.min(), values.max()
    if not high > low > 0:
        return np.array([])
    return np.geomspace(low, high, n_levels)


def contour_field(data, n_levels=CONTOUR_LEVELS,
                  floor_percentile=CONTOUR_FLOOR_PERCENTILE):
    """Prepares an image for contouring, around its saturated cores.

    The bright nuclei in this imaging are saturated, and the pipeline writes
    those pixels as exactly 0. A log stretch cannot hold zero, so masking
    them opens a hole -- up to 0.8" across in MIRI F1500W -- at the very
    peak of the source. Contour levels above the rim of that hole then trace
    its edge instead of closing on the peak, which reads as the source being
    centred somewhere it is not. The hole centres sit ~0.2" off the true
    peak in the MIRI bands, so the artefact is comparable to the real
    cross-instrument astrometric offsets it would be confused with.

    Two things are needed to fix it, and neither suffices alone. Enclosed
    holes are filled, so the level sets close around them rather than
    ringing them. And the top level is capped at the faintest rim value of
    any filled hole, so no contour is ever drawn inside one: every line on
    the figure still traces measured data, and the value used for the fill
    cannot influence where any contour goes.

    Only holes lying wholly inside the contoured range are touched. A gap in
    the faint outskirts has no contour near it to distort, and filling it
    would invent a bright source rather than repair one.

    Returns (field, levels); field is None when there is nothing to draw.
    """
    data = np.asarray(data, dtype=float)
    valid = np.isfinite(data) & (data > 0)
    values = data[valid]
    levels = log_contour_levels(values, n_levels, floor_percentile)
    if len(levels) == 0:
        return None, levels

    field = np.where(valid, data, np.nan)
    holes = ndi.binary_fill_holes(valid) & ~valid
    if not holes.any():
        return field, levels

    labels, n_holes = ndi.label(holes)
    index = np.arange(1, n_holes + 1)
    # Tag each hole's rim -- the valid pixels touching it -- with the hole's
    # own label, so every hole is judged by its own surroundings.
    rim = np.where(valid, ndi.grey_dilation(labels, size=3), 0)
    rim_min = np.array(ndi.minimum(np.where(valid, data, np.inf), rim, index))

    embedded = rim_min >= levels[0]
    if not embedded.any():
        return field, levels

    field[np.isin(labels, index[embedded])] = values.max()
    # Drop the levels saturation prevents us from tracing, rather than
    # re-spacing the set into the range that survives. Keeping the original
    # spacing leaves every outer contour exactly where it would be on an
    # unsaturated image; re-spacing would pull the whole set down and crowd
    # extra levels into the noise, changing the figure far from the core to
    # fix a defect at it.
    levels = levels[levels <= rim_min[embedded].min()]
    if len(levels) == 0:
        return None, levels
    return field, levels


# Loaded contour images, keyed by path. These mosaics run to ~90 MB, and a
# rendering run opens the same one for every map of every line, so the read
# happens once per process rather than 175 times.
_CONTOUR_IMAGE_CACHE = {}


def load_contour_image(path):
    """Reads (and remembers) the science array and WCS of a contour image."""
    if path not in _CONTOUR_IMAGE_CACHE:
        with fits.open(path) as hdul:
            _CONTOUR_IMAGE_CACHE[path] = (hdul["SCI"].data,
                                          WCS(hdul["SCI"].header).celestial)
    return _CONTOUR_IMAGE_CACHE[path]


MODES = ["cube", "spaxel", "spaxels"]
TEST_SPAXELS = [(32, 10), (26, 16), (31, 19), (26, 36), (25, 21), (28, 38),
                (25, 22), (35, 12), (29, 10), (29, 32), (22, 26), (26, 22),
                (19, 26), (26, 19), (31, 14)]

# Number of table columns per spaxel after XPIX/YPIX:
# G1AMP..G3SIGMA (9) + REDCHI + NCOMP
_ROW_LENGTH = 11


### Per-spaxel fitting machinery
#
# These are module-level (rather than Spaxelcube methods) so that
# ProcessPoolExecutor can pickle them by reference and ship individual
# spaxel tasks to worker processes. Each task is self-contained: plain
# numpy arrays plus a config dict of per-line constants precomputed once
# by Spaxelcube._build_fit_config().


def _mask_spaxel_regions(spaxel, cfg):
    """Extracts the continuum and line spectral regions from a spaxel."""
    unit = u.Unit(cfg["wv_unit"])
    lower_continuum_region = SpectralRegion(cfg["lower_continuum"][0] * unit,
                                            cfg["lower_continuum"][1] * unit)
    upper_continuum_region = SpectralRegion(cfg["upper_continuum"][0] * unit,
                                            cfg["upper_continuum"][1] * unit)
    continuum_region = lower_continuum_region + upper_continuum_region
    line_region = SpectralRegion(cfg["line_region"][0] * unit,
                                 cfg["line_region"][1] * unit)
    spaxel_continuum = extract_region(spaxel, continuum_region, return_single_spectrum=True)
    spaxel_line = extract_region(spaxel, line_region)
    return spaxel_continuum, spaxel_line


def _background_subtraction(spaxel_continuum, spaxel_line):
    """Performs background subtraction on a spaxel.

    If `spaxel_line` carries an uncertainty, it is propagated through
    the subtraction (the continuum model itself is treated as exact)
    and returned as `line_flux_errors`; otherwise `line_flux_errors`
    is None.
    """
    nanidx = np.where(np.isfinite(spaxel_continuum.flux.value))[0]
    popt, pcov = fit_continuum(spaxel_continuum.spectral_axis.value[nanidx],
                               spaxel_continuum.flux.value[nanidx])
    spaxel_continuum_fit = Spectrum1D(OneDPolynomial(spaxel_line.spectral_axis.value, popt[0], popt[1]) * u.MJy / u.steradian,
                                        spaxel_line.spectral_axis)
    spaxel_continuum_sub = spaxel_line - spaxel_continuum_fit
    line_fluxes = spaxel_continuum_sub.flux.value
    line_wavelengths = spaxel_continuum_sub.spectral_axis.value
    if spaxel_continuum_sub.uncertainty is not None:
        line_flux_errors = spaxel_continuum_sub.uncertainty.quantity.to(spaxel_continuum_sub.flux.unit).value
    else:
        line_flux_errors = None
    return line_wavelengths, line_fluxes, spaxel_continuum_fit.flux.value, line_flux_errors


def _empirical_noise_scale(spaxel_continuum, min_points=8):
    """Estimates a per-spaxel noise floor from the line-free continuum
    bands: the factor by which the pipeline ERR extension underestimates
    the actual scatter of the data around the fitted continuum.

    JWST pipeline errors commonly underestimate the effective per-channel
    noise (correlated noise, fringing residuals), which leaves the
    absolute chi-square scale of the fits uncalibrated. Because the
    continuum bands contain no line, their standardized residuals
    (data - linear continuum) / ERR should scatter with a standard
    deviation of 1 if ERR is correct; the excess factor is measured here
    with a robust MAD estimator (insensitive to stray weak lines or hot
    pixels inside the band) and used to inflate the errors before
    fitting.

    Returns a scale factor >= 1 (ERR is never deflated -- if the pipeline
    overestimates the noise, the conservative value is kept), or 1.0 when
    the estimate is not possible (too few points, missing uncertainty).
    """
    if spaxel_continuum.uncertainty is None:
        return 1.0
    flux = spaxel_continuum.flux.value
    wave = spaxel_continuum.spectral_axis.value
    err = spaxel_continuum.uncertainty.quantity.to(spaxel_continuum.flux.unit).value
    good = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    if good.sum() < min_points:
        return 1.0
    try:
        popt, pcov = fit_continuum(wave[good], flux[good])
    except Exception:
        return 1.0
    standardized = (flux[good] - OneDPolynomial(wave[good], popt[0], popt[1])) / err[good]
    scale = 1.4826 * np.nanmedian(np.abs(standardized - np.nanmedian(standardized)))
    if not np.isfinite(scale):
        return 1.0
    return max(1.0, float(scale))


def _amplitude_scale(line_fluxes):
    """Returns a positive flux scale for seeding amplitude guesses/bounds."""
    peak = np.nanmax(line_fluxes)
    if not np.isfinite(peak) or peak <= 0:
        peak = np.nanmax(np.abs(line_fluxes))
    if not np.isfinite(peak) or peak <= 0:
        peak = 1.0
    return peak


def _fit_single_gaussian(line_wavelengths, line_fluxes, weights, cfg):
    """Fits the single-Gaussian model. Amplitude guesses and bounds are
    scaled from the data (area of a Gaussian with the observed peak height)
    rather than hardcoded."""
    line_center = cfg["line_center"]
    scale = _amplitude_scale(line_fluxes)
    amp_guess = np.sqrt(2 * np.pi) * scale * cfg["narrow_guess_sigma"]
    amp_max = np.sqrt(2 * np.pi) * scale * cfg["maximum_sigma"] * 4
    model = GaussianModel(prefix="g1_")
    params = model.make_params(
        center=dict(value=line_center,
                    min=line_center - cfg["center_offset"],
                    max=line_center + cfg["center_offset"]),
        sigma=dict(value=cfg["narrow_guess_sigma"],
                   min=cfg["minimum_sigma"], max=cfg["maximum_sigma"]),
        amplitude=dict(value=amp_guess, min=0.0, max=amp_max))
    return model.fit(line_fluxes, params, x=line_wavelengths, weights=weights,
                     max_nfev=cfg["max_nfev"],
                     scale_covar=cfg.get("scale_covar", True))


def _fit_double_gaussian_ladder(line_wavelengths, line_fluxes, weights, single_result, cfg):
    """Fits the double-Gaussian model from several qualitatively distinct
    starting guesses and returns the lowest chi-square result (or None if
    every start failed).

    The chi-square surface of a bounded two-Gaussian model has a small
    number of distinct basins corresponding to different line topologies,
    so one seed per topology (plus one data-driven seed) protects against
    local minima at a fraction of the cost of random multistart or global
    optimizers:

    1. residual-seeded: the second component starts at the peak of the
       single-Gaussian fit's |residual|, i.e. wherever the one-component
       model most underperforms;
    2. narrow core + broad wing at the line center (outflow wing);
    3. narrow + narrow, second component redshifted;
    4. narrow + narrow, second component blueshifted.

    The primary component's center is bounded within +/-`center_offset`
    of the catalog line center (systemic), while the secondary may roam
    within the wider +/-`center_cutoff` (shifted outflow components).
    `max_nfev` bounds the cost of degenerate fits: when the data holds no
    real second component, the chi-square valley is nearly flat and the
    optimizer can otherwise crawl for thousands of evaluations -- exactly
    the spaxels the F-test will reject anyway.
    """
    line_center = cfg["line_center"]
    offset = cfg["center_offset"]
    cutoff = max(cfg["center_cutoff"], offset)
    narrow_sigma = cfg["narrow_guess_sigma"]
    broad_sigma = cfg["broad_guess_sigma"]
    scale = _amplitude_scale(line_fluxes)
    narrow_amp = np.sqrt(2 * np.pi) * scale * narrow_sigma
    broad_amp = np.sqrt(2 * np.pi) * scale * broad_sigma
    amp_max = np.sqrt(2 * np.pi) * scale * cfg["maximum_sigma"] * 4

    # (g1_center, g1_sigma, g1_amp, g2_center, g2_sigma, g2_amp) per start
    residual = line_fluxes - single_result.best_fit
    residual_idx = int(np.nanargmax(np.abs(residual)))
    residual_center = float(np.clip(line_wavelengths[residual_idx],
                                    line_center - cutoff, line_center + cutoff))
    residual_amp = np.sqrt(2 * np.pi) * max(np.abs(residual[residual_idx]), 0.1 * scale) * narrow_sigma
    starts = [
        (line_center, narrow_sigma, narrow_amp, residual_center, narrow_sigma, residual_amp),
        (line_center, narrow_sigma, narrow_amp, line_center, broad_sigma, broad_amp),
        (line_center, narrow_sigma, narrow_amp, line_center + offset, narrow_sigma, narrow_amp),
        (line_center, narrow_sigma, narrow_amp, line_center - offset, narrow_sigma, narrow_amp),
    ]

    best_result = None
    for g1_center, g1_sigma, g1_amp, g2_center, g2_sigma, g2_amp in starts:
        gauss_1 = GaussianModel(prefix="g1_")
        gauss_2 = GaussianModel(prefix="g2_")
        params = gauss_1.make_params(
            center=dict(value=g1_center, min=line_center - offset, max=line_center + offset),
            sigma=dict(value=g1_sigma, min=cfg["minimum_sigma"], max=cfg["maximum_sigma"]),
            amplitude=dict(value=min(g1_amp, amp_max), min=0.0, max=amp_max))
        params.update(gauss_2.make_params(
            center=dict(value=g2_center, min=line_center - cutoff, max=line_center + cutoff),
            sigma=dict(value=g2_sigma, min=cfg["minimum_sigma"], max=cfg["maximum_sigma"]),
            amplitude=dict(value=min(g2_amp, amp_max), min=0.0, max=amp_max)))
        try:
            result = (gauss_1 + gauss_2).fit(line_fluxes, params, x=line_wavelengths,
                                             weights=weights, max_nfev=cfg["max_nfev"],
                                             scale_covar=cfg.get("scale_covar", True))
        except Exception:
            continue
        if best_result is None or result.chisqr < best_result.chisqr:
            best_result = result
    return best_result


def _double_model_preferred_test(single_result, double_result, alpha=0.05):
    """F-test on the nested single/double-Gaussian fits; see
    Spaxelcube._double_model_preferred for the full rationale."""
    dof_single = single_result.nfree
    dof_double = double_result.nfree
    chi2_single = single_result.chisqr
    chi2_double = double_result.chisqr

    extra_params = dof_single - dof_double
    if extra_params <= 0 or dof_double <= 0:
        return False
    if chi2_double >= chi2_single:
        return False

    f_stat = ((chi2_single - chi2_double) / extra_params) / (chi2_double / dof_double)
    p_value = f_dist.sf(f_stat, extra_params, dof_double)
    return p_value < alpha


def _validate_components(fit_result, snr_threshold=3.0):
    """Keeps fitted components whose amplitude / amplitude-uncertainty
    ratio exceeds `snr_threshold` (meaningful because the fits are
    weighted by the real per-pixel flux errors, so lmfit's stderr comes
    from a calibrated covariance matrix).

    Surviving components are sorted narrow-to-broad (ascending sigma) so
    the G1/G2/G3 output columns keep a consistent physical identity across
    spaxels; with multiple fit starts, which lmfit prefix a component
    lands in is otherwise arbitrary (label switching).
    """
    components = []
    for idx in range(3):
        prefix = f"g{idx+1}_"
        if prefix + "amplitude" not in fit_result.best_values:
            continue
        amp_param = fit_result.params[prefix + "amplitude"]
        if amp_param.stderr is None or amp_param.stderr <= 0:
            continue
        if amp_param.value / amp_param.stderr > snr_threshold:
            components.append((fit_result.best_values[prefix + "amplitude"],
                               fit_result.best_values[prefix + "center"],
                               fit_result.best_values[prefix + "sigma"]))
    components.sort(key=lambda component: component[2])
    amplitudes = [component[0] for component in components]
    centers = [component[1] for component in components]
    sigmas = [component[2] for component in components]
    return amplitudes, centers, sigmas


def _fit_spaxel_arrays(flux, flux_err, dq_flag, wvs, cfg, return_details=False):
    """Fits one spaxel from plain arrays. Returns the 11-element table row
    (G1AMP..G3SIGMA, REDCHI, NCOMP), all-NaN if the spaxel has bad data or
    no significant components. With `return_details=True`, also returns a
    dict of intermediates for diagnostic plotting (None if fitting bailed
    out early)."""
    def _done(row, details=None):
        return (row, details) if return_details else row

    nan_row = [np.nan] * _ROW_LENGTH
    if dq_flag or len(flux) == 0:
        return _done(nan_row)

    flux_unit = u.Unit(cfg.get("flux_unit", "MJy/sr"))
    wv_unit = u.Unit(cfg["wv_unit"])
    good = np.isfinite(flux) & np.isfinite(flux_err) & (flux_err > 0)
    good_idx = np.where(good)[0]
    if len(good_idx) < 6:
        return _done(nan_row)
    spaxel = Spectrum1D(flux[good_idx] * flux_unit, wvs[good_idx] * wv_unit,
                        uncertainty=StdDevUncertainty(flux_err[good_idx] * flux_unit))

    # Extraction/subtraction can fail on pathological spaxels (e.g. an
    # FOV-edge spaxel whose continuum band is entirely NaN, leaving an
    # empty sub-region); such spaxels are unfittable, so report NaN.
    try:
        spaxel_continuum, spaxel_line = _mask_spaxel_regions(spaxel, cfg)
        if len(spaxel_line.flux) < 6 or len(spaxel_continuum.flux) < 6:
            return _done(nan_row)
        line_wavelengths, line_fluxes, background_flux, line_flux_errors = _background_subtraction(spaxel_continuum, spaxel_line)
    except Exception:
        return _done(nan_row)
    # Empirical noise floor: inflate the propagated errors by the excess
    # scatter measured in the line-free continuum bands, so the fit
    # weights (and hence chi-square and parameter uncertainties) are
    # calibrated against the actual noise rather than the pipeline ERR.
    if cfg.get("noise_floor", True):
        noise_scale = _empirical_noise_scale(spaxel_continuum)
    else:
        noise_scale = 1.0
    line_flux_errors = line_flux_errors * noise_scale
    finite_idx = np.where(np.isfinite(line_fluxes) & np.isfinite(line_flux_errors) & (line_flux_errors > 0))[0]
    line_wavelengths = line_wavelengths[finite_idx]
    line_fluxes = line_fluxes[finite_idx]
    background_flux = background_flux[finite_idx]
    line_flux_errors = line_flux_errors[finite_idx]
    if len(line_fluxes) < 6:
        return _done(nan_row)
    weights = 1.0 / line_flux_errors

    try:
        single_result = _fit_single_gaussian(line_wavelengths, line_fluxes, weights, cfg)
    except Exception:
        return _done(nan_row)
    double_result = _fit_double_gaussian_ladder(line_wavelengths, line_fluxes,
                                                weights, single_result, cfg)

    prefer_double = False
    if double_result is not None:
        try:
            prefer_double = _double_model_preferred_test(single_result, double_result,
                                                         cfg["alpha"])
        except Exception:
            prefer_double = False
    best_result = double_result if prefer_double else single_result

    amplitudes, centers, sigmas = _validate_components(best_result, cfg["snr_threshold"])
    row = list(nan_row)
    if len(amplitudes) > 0:
        for comp_idx in range(len(amplitudes)):
            row[3 * comp_idx] = amplitudes[comp_idx]
            row[3 * comp_idx + 1] = centers[comp_idx]
            row[3 * comp_idx + 2] = sigmas[comp_idx]
        row[9] = best_result.redchi
        row[10] = len(amplitudes)

    details = None
    if return_details:
        details = {"spaxel_line": spaxel_line, "background_flux": background_flux,
                   "line_wavelengths": line_wavelengths, "noise_scale": noise_scale,
                   "single_result": single_result, "double_result": double_result}
    return _done(row, details)


def _fit_spaxel_task(task):
    """Worker entry point for multiprocessing: unpacks one spaxel task and
    returns (x_pix, y_pix, row)."""
    x_pix, y_pix, flux, flux_err, dq_flag, wvs, cfg = task
    return x_pix, y_pix, _fit_spaxel_arrays(flux, flux_err, dq_flag, wvs, cfg)


class Spaxelcube(Datacube):
    """
    A Datacube that performs multicomponent Gaussian fits to its own IFS
    spectra and exposes the fitted parameter maps as a derived product.

    A Spaxelcube adopts the loaded state of a source Datacube (headers,
    data, WCS, wavelength grid, unit helpers), so all of the inherited
    Datacube methods -- `spaxel_values`, `vel_to_wv`, `mrs_resolving_power`,
    `wcs`, etc. -- operate on this object directly rather than on a wrapped
    attribute.
    """

    def __init__(self, datacube, name, output, plot_output=None,
                 palette=DEFAULT_PALETTE, scaling=DEFAULT_SCALING,
                 mode="cube", test_spaxel=(20,20), test_spaxels=None,
                 distance=194.99 * u.Mpc, contour_image=None):
        """Inits Spaxelcube.

        Args:
            datacube (Datacube): The source datacube to fit. Its loaded state
                (headers, data, WCS, wavelength grid, unit helpers) is adopted
                into this Spaxelcube so the inherited Datacube methods operate
                on the fitted cube directly.
            name (str): The name of the emission line being fitted.
            output (str): The path to save data products (the fit table, the
                FITS parameter cube, q3dfit's per-spaxel files).
            plot_output (str, optional): The directory rendered figures are
                written to, kept separate from ``output`` so the data products
                and the plots can live in different trees (the repository
                splits them into ``spaxs/`` and ``plots/line_maps/``). Defaults
                to ``output``, i.e. figures land alongside the data.
            palette (str): ``light`` (default) or ``dark``; the palette every
                figure this object renders is drawn in. It is re-asserted onto
                the global rcParams before each figure, so it holds regardless
                of what else in the process has styled a plot, and it tags the
                saved filenames (``..._light.png`` / ``..._dark.png``) so both
                palettes of a map coexist in one directory.
            scaling (str): ``paper``, ``presentation`` (default), or ``poster``;
                the font/DPI scaling passed to ``PlotParams``.
            mode (str): This should take on either ``cube``, ``spaxel`` or
                ``spaxels``. ``cube`` is the default and fits the entire
                datacube. ``spaxel`` fits a designated spaxel. ``spaxels`` fits
                a list of predetermined spaxels.
            test_spaxel (tuple): The coordinates of the spaxel to fit in
                ``spaxel`` mode. The default is (20, 20).
            test_spaxels (list): A list of tuples corresponding to the spaxels
                to trial in ``spaxels`` mode.
            distance (Quantity): The distance to the target, used for the
                physical scalebar on rendered maps (default: 194.99 Mpc,
                IR 23128).
            contour_image (str, optional): A key of ``CONTOUR_IMAGES``
                (``f200w``, ``f356w``) or a path to an imaging FITS file.
                When given, the maps are overlaid with that image's intensity
                contours *instead of* the cube's footprint outline, and their
                filenames are tagged with it so the two variations coexist in
                one run directory. Default None: outline the footprint, which
                is the standard map.
        """

        if mode not in MODES:
            raise ValueError("Mode must be 'cube', 'spaxel', or 'spaxels'")
        if test_spaxels == None:
            self.test_spaxels = TEST_SPAXELS

        # Adopt the source cube's loaded state so the inherited Datacube
        # helper methods (spaxel_values, vel_to_wv, wcs, ...) operate on
        # this object directly instead of a wrapped attribute.
        self.__dict__.update(datacube.__dict__)
        self.name = name
        self.output = output
        self.plot_output = output if plot_output is None else plot_output
        self.pltparams = PlotParams(palatte=palette, scaling=scaling)
        self.mode = mode
        self.test_spaxel = test_spaxel
        self.distance = distance
        # Resolved once here so a bad name or path fails when the object is
        # built, not seven maps into a rendering run.
        self.contour_image = self._resolve_contour_image(contour_image)
        self.line_dict = read_line_params()
        self.line_center = self.line_dict[self.name][4]
        self.label = None

    @staticmethod
    def _resolve_contour_image(contour_image):
        """Turns a contour-image name or path into (tag, path), or None.

        The tag is what the saved filenames carry, so a map overlaid with
        F200W contours cannot overwrite the footprint version of itself, nor
        the F356W one.
        """
        if contour_image is None:
            return None
        key = str(contour_image).lower()
        if key in CONTOUR_IMAGES:
            return key, str(paths.image_data_dir() / CONTOUR_IMAGES[key])
        if not os.path.exists(contour_image):
            raise FileNotFoundError(
                f"No contour image {contour_image!r}; expected a path or one "
                f"of {sorted(CONTOUR_IMAGES)}")
        tag = os.path.splitext(os.path.basename(contour_image))[0].lower()
        return tag, contour_image


    # I/O methods
    
    def load_fit(self, filepath):
        """Loads a saved instance of Spaxelcube"""
        self.fitparams = ascii.read(filepath, format="ipac") 
    
    
    #def jwst_miri_broadening(self, wavelength):
    #    """Estimates the JWST instrumental broadening of the MIRI MRS
    #    intrument at a given wavelength (in microns). Returns as a tuple:
    #    (FWHM, velocity dispersion)."""
    #    fwhm = wavelength / (4603 - 128 * wavelength) / 2.3548
    #    vel_disp = self.sigma_to_disp(fwhm, wavelength)
    #    return fwhm, vel_disp
    
    
    
    
    
    def _mask_spaxel(self, spaxel):
        """Develops the masks necessary to perform gaussian fitting."""
        line_params = self.line_dict[self.name]
        cfg = {"wv_unit": str(self._wv_unit),
               "lower_continuum": (line_params[2][0], line_params[3][0]),
               "upper_continuum": (line_params[2][1], line_params[3][1]),
               "line_region": (line_params[2][0], line_params[2][1])}
        return _mask_spaxel_regions(spaxel, cfg)
    
    
    def clean_spaxel(self, spaxel):
        """Removes troublesome nan values from a spaxel."""
        nanidx = np.where(np.isfinite(spaxel.flux.value))[0]
        #print(spaxel.flux.value)
        #print(spaxel.flux.value[nanidx])

    def background_subtraction(self, spaxel_continuum, spaxel_line):
        """Performs background subtraction on a spaxel.

        If `spaxel_line` carries an uncertainty, it is propagated through
        the subtraction (the continuum model itself is treated as exact)
        and returned as `line_flux_errors`; otherwise `line_flux_errors`
        is None.
        """
        return _background_subtraction(spaxel_continuum, spaxel_line)
    
    
    
    '''def two_gaussian_fit(self, diagnose=False):
        """ 
        A simple fitting routine which fits two gaussian 
        components to a line.
        """
        
        # Bookkeeping variables and data tracking.
        names = ["amplitude", "center", "sigma"]
        prefixes = ["g1_", "g2_", "g3_"]
        best_fit_values = [[], [], [], [], [], [], [], [], [], [], [], [], []]
        
        # Iterate through all spaxels
        for y_pix in tqdm(range(self.im_shape[1])):
            for x_pix in tqdm(range(self.im_shape[0]), leave=False):
                
                # Record X and Y pixels
                best_fit_values[0].append(x_pix)
                best_fit_values[1].append(y_pix)
                
                # Test for mode case
                if self.mode == "spaxels":
                    if (x_pix, y_pix) not in self.test_spaxels:
                        continue
                if self.mode == "spaxel":
                    if (x_pix, y_pix) != self.test_spaxel:
                        continue
                
                # Record null fit values for spaxels with bad data
                (flux, flux_err, dq) = self.spaxel_values(x_pix, y_pix)
                if 513 in np.unique(dq[0]):
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue
                
                #if x_pix == 10 and y_pix == 10:
                #    fig, ax = plt.subplots()
                #    ax.plot(self.wvs.value, flux.value)
                #    plt.show()
                """# Added for line fitting case
                if np.isnan(flux).any():
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue"""
                #print(len(flux))
                #print(len(self.wvs))
                
                # Report null values for data with zero-length arrays
                if len(flux) == 0:
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue
                nanidx = np.where(np.isfinite(flux))[0]
                spaxel = Spectrum1D(flux[nanidx], self.wvs[nanidx])
                if len(spaxel.flux) < 6:
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue
                
                
                # Process spaxel
                spaxel_continuum, spaxel_line = self._mask_spaxel(spaxel)
                #self.clean_spaxel(spaxel_line)
                if len(spaxel_line.flux) < 6 or len(spaxel_continuum.flux) < 6:
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue
                line_wavelengths, line_fluxes, background_flux = self.background_subtraction(spaxel_continuum, spaxel_line)
                nanidx = np.where(np.isfinite(line_fluxes))[0]
                #line_wavelengths = line_wavelengths[nanidx]
                #line_fluxes = line_fluxes[nanidx]
                #background_flux = background_flux[nanidx]
                
                """# Estimate line parameters
                
                line_center = self.line_dict[self.name][4]
                resolving_power = self.mrs_resolving_power(line_center)
                #
                center_offset = self.vel_to_wv(self.line_dict[self.name][5], line_center).value
                center_cutoff = self.vel_to_wv(self.line_dict[self.name][6], line_center).value
                
                #print(line_center / resolving_power)
                #minimum_sigma = self.disp_to_vel(self.line_dict[self.name][7], line_center).value
                minimum_sigma = line_center / resolving_power
                narrow_sigma = self.disp_to_vel(self.line_dict[self.name][8], line_center).value
                broad_sigma = self.disp_to_vel(self.line_dict[self.name][9], line_center).value
                maximum_sigma = self.disp_to_vel(self.line_dict[self.name][10], line_center).value
                #print(minimum_sigma, narrow_sigma, broad_sigma, maximum_sigma)
                #print(self.sigma_to_disp(minimum_sigma, line_center), self.line_dict[self.name][8], self.line_dict[self.name][9], self.line_dict[self.name][10])
                narrow_broad_sigma = (minimum_sigma + maximum_sigma) / 2
                minimum_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*minimum_sigma
                narrow_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*narrow_sigma
                broad_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*broad_sigma
                maximum_amp = np.sqrt(2*np.pi)*np.max(line_fluxes)*maximum_sigma * 4
                narrow_broad_amp = (minimum_amp + maximum_amp) / 2
                
                #test_relvel = 800 * u.kilometer / u.second
                #test_sigma = 500 * u.kilometer / u.second
                #print(self.vel_to_wv(test_relvel, line_center))
                #print(self.disp_to_vel(test_sigma, line_center))"""
                
                # Estimate line parameters
                
                line_center = self.line_dict[self.name][4]
                resolving_power = self.mrs_resolving_power(line_center)
                
                # Wavelength parameters
                center_offset = self.vel_to_wv(self.line_dict[self.name][5], line_center).value
                center_cutoff = self.vel_to_wv(self.line_dict[self.name][6], line_center).value
                
                # Line dispersion parameters
                minimum_sigma = line_center / resolving_power
                cutoff_sigma = self.vel_to_sigma(self.line_dict[self.name][8], line_center).value
                maximum_sigma = self.vel_to_sigma(self.line_dict[self.name][10], line_center).value
                narrow_guess_sigma = (minimum_sigma + cutoff_sigma) / 2
                broad_guess_sigma = (maximum_sigma + cutoff_sigma) / 2

                # Amplitude parameters
                minimum_amp = 0.001#np.sqrt(2*np.pi)*np.max(line_fluxes)*minimum_sigma
                narrow_amp = 10#np.sqrt(2*np.pi)*np.max(line_fluxes)*narrow_guess_sigma
                broad_amp = 200#np.sqrt(2*np.pi)*np.max(line_fluxes)*broad_guess_sigma
                maximum_amp = 10000#np.sqrt(2*np.pi)*np.max(line_fluxes)*maximum_sigma * 4
                
                ### Run lmfit fitting routine ###
                
                # Single Gaussian
                single_g1 = GaussianModel(prefix="g1_")
                single_params = single_g1.guess(line_fluxes, x=line_wavelengths)
                single_params.update(single_g1.make_params(center=dict(value=line_center, min=line_center-center_offset, max=line_center+center_offset),
                                            sigma=dict(value=narrow_guess_sigma, min=minimum_sigma, max=maximum_sigma),
                                            amplitude=dict(value=narrow_amp, max=maximum_amp, min=minimum_amp)))
                single_result = single_g1.fit(line_fluxes, single_params, x=line_wavelengths)
                single_redchi = single_result.redchi
                
                # Double gaussian
                double_g1_2 = GaussianModel(prefix="g1_")
                double_g2_2 = GaussianModel(prefix="g2_")
                double_params_2 = double_g1_2.guess(line_fluxes, x=line_wavelengths)
                double_params_2.update(double_g1_2.make_params(center=dict(value=line_center, min=line_center-center_cutoff, max=line_center+center_cutoff),
                                            sigma=dict(value=narrow_guess_sigma, min=minimum_sigma, max=maximum_sigma),
                                            amplitude=dict(value=narrow_amp, max=maximum_amp*100, min=minimum_amp*0.01)))
                double_params_2.update(double_g2_2.make_params(center=dict(value=line_center, min=line_center-center_cutoff, max=line_center+center_cutoff),
                                            sigma=dict(value=narrow_guess_sigma, min=minimum_sigma, max=maximum_sigma),
                                            amplitude=dict(value=narrow_amp, max=maximum_amp*100, min=minimum_amp*0.01)))
                double_model_2 = double_g1_2 + double_g2_2
                double_result = double_model_2.fit(line_fluxes, double_params_2, x=line_wavelengths)
                
                
                # Record reduced chi squared values
                try:
                    single_redchi = single_result.redchi
                except:
                    single_redchi = 1e7
                try:
                    double_redchi = double_result.redchi
                except:
                    double_redchi = 1e7
                
                # Select the model with the lowest chi square
                if single_redchi <= double_redchi:
                    best_result = single_result
                    best_redchi = single_redchi    
                else:
                    best_result = double_result
                    best_redchi = double_redchi  
                
                amplitudes, centers, sigmas = self.validate_snr(best_result, 
                                                                spaxel_line, 
                                                                background_flux, 
                                                                line_wavelengths)
                
                best_idx = 2
                if len(amplitudes) == 0:
                    for idx in range(2,13):
                        best_fit_values[idx].append(np.nan)
                    continue
                for idx, amp in enumerate(amplitudes):
                    best_fit_values[best_idx + idx * 3].append(amp)
                    best_fit_values[best_idx + idx * 3 + 1].append(centers[idx])
                    best_fit_values[best_idx + idx * 3 + 2].append(sigmas[idx])
                last_idx = best_idx + idx * 3 + 2
                for idx in range(last_idx+1,11):
                    best_fit_values[idx].append(np.nan)
                """for prefix in prefixes:
                    for pname in names:
                        if prefix+pname in best_result.best_values:
                            best_fit_values[best_idx].append(best_result.best_values[prefix+pname])
                            best_idx += 1
                        else:
                            best_fit_values[best_idx].append(np.nan)
                            best_idx += 1"""
                best_fit_values[11].append(best_redchi)
                best_fit_values[12].append(len(amplitudes))
                
                
                
        if self.mode == "cube":
            #for listy in best_fit_values:
            #    print(len(listy))
            self.fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA", "REDCHI", "NCOMP"))
            self.fitparams.write(self.output + "twogaussian_raw.dat", format="ipac", overwrite=True)
            self.label = "twogaussian" '''

    def _double_model_preferred(self, single_result, double_result, alpha=0.05):
        """Decides whether the double-Gaussian fit is a statistically
        justified improvement over the single-Gaussian fit.

        A double Gaussian has strictly more free parameters than a single
        Gaussian, so it will nearly always achieve an equal-or-lower raw
        chi-square even when the extra component is just fitting noise.
        Comparing `redchi` directly is therefore biased toward the more
        complex model. Instead, this runs an F-test on the two *nested*
        models: it asks whether the reduction in chi-square from adding
        the second component is larger than what the extra degrees of
        freedom would explain by chance alone, at the given significance
        level `alpha`.

        Args:
            single_result (lmfit.model.ModelResult): The single-Gaussian fit
                result.
            double_result (lmfit.model.ModelResult): The double-Gaussian fit
                result.
            alpha (float): The significance level required to prefer the
                double-Gaussian model (default 0.05, i.e. the improvement must
                be significant at >95% confidence).

        Returns:
            bool: ``True`` if the double-Gaussian model is preferred.
        """
        return _double_model_preferred_test(single_result, double_result, alpha)

    def _build_fit_config(self, snr_threshold=3.0, alpha=0.05, max_nfev=2000,
                          noise_floor=True):
        """Precomputes the per-line fitting constants shared by every
        spaxel (they only depend on the line, not the spaxel data).
        Returned as a plain dict so it can be shipped to worker processes.

        Args:
            snr_threshold (float): The amplitude/uncertainty ratio a fitted
                component must exceed to be kept (see ``validate_snr``).
            alpha (float): The F-test significance level required to prefer the
                double-Gaussian model (see ``_double_model_preferred``).
            max_nfev (int): Cap on optimizer function evaluations per fit.
                Healthy fits use well under 200; the cap only truncates
                degenerate chi-square valley crawls.
            noise_floor (bool): Rescale each spaxel's errors by the excess
                scatter measured in its line-free continuum bands before
                fitting (see ``_empirical_noise_scale``). When enabled, lmfit's
                post-hoc covariance rescaling (``scale_covar``) is turned off,
                so parameter uncertainties come from the continuum-calibrated
                weights rather than from the line-fit residuals (which conflate
                model mismatch with noise).
        """
        line_params = self.line_dict[self.name]
        line_center = line_params[4]
        resolving_power = self.mrs_resolving_power(line_center)
        minimum_sigma = line_center / resolving_power / 2
        cutoff_sigma = self.vel_to_sigma(line_params[8], line_center).value
        maximum_sigma = self.vel_to_sigma(line_params[10], line_center).value
        return {
            "wv_unit": str(self._wv_unit),
            "lower_continuum": (line_params[2][0], line_params[3][0]),
            "upper_continuum": (line_params[2][1], line_params[3][1]),
            "line_region": (line_params[2][0], line_params[2][1]),
            "line_center": line_center,
            "center_offset": self.vel_to_wv(line_params[5], line_center).value,
            "center_cutoff": self.vel_to_wv(line_params[6], line_center).value,
            "minimum_sigma": minimum_sigma,
            "cutoff_sigma": cutoff_sigma,
            "maximum_sigma": maximum_sigma,
            "narrow_guess_sigma": (minimum_sigma + cutoff_sigma) / 2,
            "broad_guess_sigma": (maximum_sigma + cutoff_sigma) / 2,
            "snr_threshold": snr_threshold,
            "alpha": alpha,
            "max_nfev": max_nfev,
            "noise_floor": noise_floor,
            "scale_covar": not noise_floor,
        }

    def two_gaussian_fit(self, diagnose=False, n_workers=None, noise_floor=True):
        """Fits each spaxel with a single- and a double-Gaussian model and
        keeps the model preferred by an F-test on the error-weighted fits.

        The double-Gaussian model is fit from several qualitatively distinct
        starting guesses to avoid local minima (see
        ``_fit_double_gaussian_ladder``).

        Args:
            diagnose (bool): Displays diagnostic fit plots for each spaxel.
                Forces serial fitting.
            n_workers (int, optional): Number of worker processes to fit
                spaxels in parallel. Defaults to the number of CPUs; pass 1 to
                disable multiprocessing. Only ``cube`` mode without
                ``diagnose`` runs in parallel.
            noise_floor (bool): Rescale each spaxel's errors by the excess
                scatter measured in its line-free continuum bands, calibrating
                the chi-square scale and parameter uncertainties against the
                actual noise (default ``True``; see
                ``_empirical_noise_scale``).
        """
        cfg = self._build_fit_config(noise_floor=noise_floor)
        wvs = self.wvs.value

        # Spaxel coordinates in the original output-table order
        # (x fastest, y slowest).
        coords = [(x_pix, y_pix)
                  for y_pix in range(self.im_shape[1])
                  for x_pix in range(self.im_shape[0])]
        if self.mode == "spaxels":
            fit_coords = [coord for coord in coords if coord in self.test_spaxels]
        elif self.mode == "spaxel":
            fit_coords = [coord for coord in coords if coord == self.test_spaxel]
        else:
            fit_coords = coords

        def build_task(x_pix, y_pix):
            flux, flux_err, dq = self.spaxel_values(x_pix, y_pix)
            dq_flag = len(dq) > 0 and 513 in np.unique(dq[0])
            if "flux_unit" not in cfg:
                cfg["flux_unit"] = str(flux.unit)
            return (x_pix, y_pix, np.asarray(flux.value, dtype=float),
                    np.asarray(flux_err.value, dtype=float), dq_flag, wvs, cfg)

        results = {}
        use_parallel = (self.mode == "cube" and not diagnose
                        and (n_workers is None or n_workers > 1))
        if use_parallel:
            tasks = [build_task(x_pix, y_pix) for x_pix, y_pix in fit_coords]
            worker_count = n_workers or os.cpu_count() or 1
            chunksize = max(1, len(tasks) // (worker_count * 8))
            try:
                with ProcessPoolExecutor(max_workers=worker_count) as executor:
                    for x_pix, y_pix, row in tqdm(
                            executor.map(_fit_spaxel_task, tasks, chunksize=chunksize),
                            total=len(tasks)):
                        results[(x_pix, y_pix)] = row
            except BrokenProcessPool:
                print("Parallel fitting failed -- this usually means the calling "
                      "script has no `if __name__ == '__main__':` guard, which "
                      "multiprocessing requires. Falling back to serial fitting.")
                results = {}
                use_parallel = False
        if not use_parallel:
            for x_pix, y_pix in tqdm(fit_coords):
                _, _, flux, flux_err, dq_flag, _, _ = build_task(x_pix, y_pix)
                if diagnose:
                    row, details = _fit_spaxel_arrays(flux, flux_err, dq_flag, wvs,
                                                      cfg, return_details=True)
                    if details is not None:
                        single_result = details["single_result"]
                        double_result = details["double_result"]
                        if double_result is not None:
                            print(single_result.redchi, double_result.redchi)
                            self.display_fitted_models(details["spaxel_line"],
                                                       details["background_flux"],
                                                       details["line_wavelengths"],
                                                       single_result, double_result)
                else:
                    row = _fit_spaxel_arrays(flux, flux_err, dq_flag, wvs, cfg)
                results[(x_pix, y_pix)] = row

        # Assemble the output columns in the original layout. Spaxels not
        # fit (the un-selected coordinates in `spaxel`/`spaxels` modes) get
        # an all-NaN row, so every mode yields a full-grid parameter table
        # (and hence full-grid FITS maps) rather than ragged columns.
        best_fit_values = [[] for _ in range(2 + _ROW_LENGTH)]
        for x_pix, y_pix in coords:
            best_fit_values[0].append(x_pix)
            best_fit_values[1].append(y_pix)
            row = results.get((x_pix, y_pix))
            if row is None:
                row = [np.nan] * _ROW_LENGTH
            for col_idx, value in enumerate(row):
                best_fit_values[2 + col_idx].append(value)

        self.fitparams = Table(best_fit_values, names=("XPIX", "YPIX", "G1AMP", "G1CEN", "G1SIGMA", "G2AMP", "G2CEN", "G2SIGMA", "G3AMP", "G3CEN", "G3SIGMA", "REDCHI", "NCOMP"))
        self.label = "twogaussian"
        # Reorder components so the outflow (furthest from line center) is
        # the second component, consistently across the table, the maps,
        # and the FITS cube (see `_order_fitparams_by_outflow`).
        self._order_fitparams_by_outflow(cfg["line_center"])
        # The IPAC table is the full-cube data product; keep it cube-only.
        if self.mode == "cube":
            self.fitparams.write(self.output + "twogaussian_raw.dat", format="ipac", overwrite=True)
        self.write_fits(cfg=cfg)

    # Alternative fitting: q3dfit
    #
    # q3dfit (https://q3dfit.readthedocs.io) is a separate PSF-decomposition
    # and emission-line fitting package tailored to JWST NIRSpec/MIRI IFU
    # cubes. `q3dfit_fit` drives it natively -- it builds a `q3din` config
    # from this cube's own metadata and writes q3dfit's per-spaxel `.npy`
    # outputs (and `q3dcollect`ed science products) into `self.output`,
    # rather than translating into `self.fitparams`. Use q3dfit's own
    # `q3dout`/`q3dpro` tooling to inspect and map those products.

    # JWST/MIRI BAND keyword -> the sub-channel letter used in q3dfit's
    # dispersion-file naming (e.g. CHANNEL 3 + BAND LONG -> "ch3c").
    _MIRI_BAND_LETTERS = {"SHORT": "a", "MEDIUM": "b", "LONG": "c"}

    def _header_value(self, key):
        """Returns `key` from the general header, falling back to the
        science header, or None if absent from both."""
        for header in (self.general_header, self.science_header):
            if key in header:
                return header[key]
        return None

    def _resolve_q3d_line(self, tol=0.01):
        """Resolves the q3dfit line label for this cube's line by matching
        `self.line_center` against q3dfit's master linelist (nearest rest
        wavelength, vacuum microns -- the convention ``self.line_center``
        already follows).

        Args:
            tol (float): Maximum allowed separation, in microns, between the
                line center and the nearest catalog line.

        Returns:
            str: The q3dfit line label (e.g. "[SIV]10.51").

        Raises:
            ValueError: If no catalog line lies within ``tol``; pass an
                explicit ``q3d_line`` to ``q3dfit_fit`` in that case.
        """
        from q3dfit.linelist import linelist
        master = linelist(waveunit="micron", vacuum=True)
        waves = np.asarray(master["lines"], dtype=float)
        idx = int(np.argmin(np.abs(waves - self.line_center)))
        if abs(waves[idx] - self.line_center) > tol:
            raise ValueError(
                f"No q3dfit line within {tol} um of {self.name} "
                f"(center {self.line_center} um); pass q3d_line explicitly.")
        return str(master["name"][idx])

    def _q3d_spect_convol(self):
        """Builds q3dfit's spectral-convolution instruction from this cube's
        JWST/MIRI header (channel + band -> e.g. {'JWST_MIRI': ['ch3c']}),
        or returns None if the instrument/band can't be mapped (so the
        caller can skip convolution)."""
        instrument = str(self._header_value("INSTRUME") or "").upper()
        if instrument != "MIRI":
            return None
        channel = self._header_value("CHANNEL")
        band = str(self._header_value("BAND") or "").upper()
        band_letter = self._MIRI_BAND_LETTERS.get(band)
        if channel is None or band_letter is None:
            return None
        return {"JWST_MIRI": [f"ch{channel}{band_letter}"]}

    def q3dfit_fit(self, q3d_line=None, maxncomp=2, siglim_gas=None,
                   contfit_order=3, spect_convol=True, sigcut=3.0,
                   ncores=1, cols=None, rows=None, collect=True, quiet=True,
                   nocrash=True, fluxnorm=1e-18):
        """Fits this cube's emission line with q3dfit, as an alternative to
        the built-in `two_gaussian_fit`.

        The `q3din` configuration is derived from this cube's own metadata:
        the FITS path and redshift, the fit window from the line's continuum
        bounds, and (for JWST/MIRI) the spectral-resolution kernel from the
        header channel/band. q3dfit reads the cube through its own reader,
        so its science/error/DQ extensions (1/2/3, matching this cube) and
        the `error=True` flag are set for it. Results are written by q3dfit
        as per-spaxel `.npy` files (and, unless `collect=False`, combined by
        `q3dcollect`) under `self.output`; nothing is written to
        `self.fitparams`.

        Which spaxels are fit follows `self.mode`: `cube` fits all spaxels,
        `spaxel` fits `self.test_spaxel`, and `spaxels` loops over
        `self.test_spaxels`. Passing `cols`/`rows` (q3dfit's own unity-offset
        column/row selectors) overrides the mode selection.

        Args:
            q3d_line (str, optional): The q3dfit line label to fit (e.g.
                "[SIV]10.51"). Defaults to the catalog line nearest
                ``self.line_center`` (see ``_resolve_q3d_line``).
            maxncomp (int): Maximum number of velocity components per spaxel
                (default 2, to mirror ``two_gaussian_fit``).
            siglim_gas (array-like, optional): Lower/upper Gaussian sigma
                bounds in km/s. Defaults to [5, line's max sigma velocity] from
                the line parameters.
            contfit_order (int): Order of the q3dfit polynomial continuum fit
                (default 3).
            spect_convol (bool): Convolve models with the instrumental LSF
                derived from the header (default True; ignored if it can't be
                mapped).
            sigcut (float): Significance cut (in sigma) below which q3dfit's
                ``checkcomp`` discards a component and refits (default 3).
            ncores (int): Number of cores for q3dfit (default 1). With
                ncores > 1, q3dfit splits the spaxels across MPI ranks via
                ``mpiexec``; this method saves the config to a ``.npy`` file
                and hands q3dfit that path (MPI ranks can't receive the
                in-memory object). Requires ``mpiexec`` (openmpi) and
                ``mpi4py`` to be installed -- a clear error is raised if
                ``mpiexec`` is not on PATH. Multi-core helps ``cube`` mode; the
                single-region ``spaxel``/``spaxels`` modes see little benefit.
            cols (int or list, optional): q3dfit unity-offset column selector;
                overrides ``self.mode``.
            rows (int or list, optional): q3dfit unity-offset row selector;
                overrides ``self.mode``.
            collect (bool): Run ``q3dcollect`` to combine the per-spaxel
                results into science products (default True; skipped in
                ``spaxels`` mode, whose scattered coordinates don't form a
                rectangular region).
            quiet (bool): Suppress q3dfit's per-spaxel fit logging (default
                True).
            nocrash (bool): Skip (rather than abort on) spaxels whose fit
                raises, and continue to the next (default True). Needed for
                full-cube runs: FOV-edge spaxels that are almost entirely
                NaN/masked leave fewer data points than the model has free
                parameters (e.g. the order-N continuum poly), which otherwise
                raises and kills the whole run. Skipped spaxels write no output
                file and so are simply absent (NaN) from the q3dcollect maps.
                Set False when fitting a single diagnostic spaxel, where a
                crash is informative.
            fluxnorm (float): Flux normalization (in the internal
                erg/s/cm^2/micron units) that q3dfit divides the spectra by
                before fitting, purely for numerical conditioning (default
                1e-18, the value used in q3dfit's JWST/MIRI example). Without
                it the flux values are ~1e-12 and the fitted amplitudes stick
                at their bounds with undefined errors, so every component is
                rejected. The plotting methods multiply the fitted fluxes back
                by ``fluxnorm`` to recover physical units.

        Returns:
            q3din: The configured q3dfit input object, for downstream use with
            q3dfit's ``q3dout``/``q3dpro`` tooling.
        """
        from q3dfit.q3din import q3din
        from q3dfit.q3df import q3dfit as run_q3dfit

        line_params = self.line_dict[self.name]
        q3d_line = q3d_line or self._resolve_q3d_line()
        # Line names carry glob metacharacters (e.g. "[KVII]"); q3dcollect
        # globs the output filenames, so sanitize the label it builds them
        # from.
        label = re.sub(r"[^A-Za-z0-9_.+-]", "", self.name) + "_q3dfit"
        logfile = os.path.join(self.output, label + "-fitlog.txt")
        # Full fit window, in the OBSERVED frame: q3dfit fits the cube's own
        # (un-de-redshifted) wavelength axis, unlike the built-in fitter,
        # which works in the rest frame. Spans the blue edge of the lower
        # continuum band to the red edge of the upper (see
        # `_build_fit_config`), shifted to observed by (1 + redshift).
        obs = 1 + self.redshift
        fitrange = [line_params[2][0] * obs, line_params[3][1] * obs]

        q3di = q3din(self.filepath, label, outdir=self.output, logfile=logfile,
                     name=self.name, zsys_gas=self.redshift, fitrange=fitrange,
                     datext=1, varext=2, dqext=3)
        # This cube's extension 2 is an error (stddev) spectrum, not variance;
        # `fluxnorm` rescales the flux for numerical conditioning of the fit.
        q3di.argsreadcube = {"error": True, "fluxnorm": fluxnorm}
        if spect_convol:
            convol = self._q3d_spect_convol()
            if convol is not None:
                q3di.spect_convol["ws_instrum"] = convol

        # Load the cube up front so ncols/nrows are set before init_linefit
        # (which otherwise lazy-loads it through a broken code path in
        # q3dfit 1.2) and so q3dcollect has the grid dimensions.
        q3di.load_cube()

        if siglim_gas is None:
            siglim_gas = [5.0, float(line_params[10])]
        q3di.init_linefit([q3d_line], linetie=q3d_line, maxncomp=maxncomp,
                          siglim_gas=np.asarray(siglim_gas, dtype=float))
        # Automatically discard insignificant components (one at a time) and
        # refit; use the continuum-residual error for the significance test,
        # since the pipeline error underestimates the real per-channel noise.
        q3di.checkcomp = True
        q3di.argscheckcomp["sigcut"] = sigcut
        q3di.argscheckcomp["subone"] = True
        q3di.perror_useresid = True

        q3di.init_contfit("fitpoly")
        q3di.argscontfit["fitord"] = contfit_order

        # Always persist the config so it can be reloaded later (e.g. by
        # `load_q3dfit` for rendering in a fresh session). Multi-core runs
        # additionally require it: q3dfit spawns `mpiexec` and reloads the
        # config from this `.npy` path in each MPI rank (it can't ship the
        # in-memory object to the ranks). A single core fits the object.
        npy_path = os.path.join(self.output, label + ".npy")
        np.save(npy_path, q3di)
        initproc = q3di
        if ncores > 1:
            if shutil.which("mpiexec") is None:
                raise RuntimeError(
                    "q3dfit_fit(ncores > 1) needs MPI, but 'mpiexec' was not "
                    "found on PATH. Install openmpi and mpi4py (e.g. add them "
                    "to the pixi environment), or call with ncores=1.")
            initproc = npy_path

        collect_cols, collect_rows = self._run_q3dfit(q3di, initproc, run_q3dfit,
                                                       cols, rows, ncores, quiet,
                                                       nocrash)
        if collect and collect_cols is not None:
            from q3dfit.q3dcollect import q3dcollect
            q3dcollect(q3di, cols=collect_cols, rows=collect_rows,
                       compsortpar="sigma", compsortdir="down")
        # Remember the fit for the render_q3d_* methods (and its flux scale).
        self.q3di = q3di
        self._q3d_fluxnorm = fluxnorm
        return q3di

    def _run_q3dfit(self, q3di, initproc, run_q3dfit, cols, rows, ncores, quiet,
                    nocrash):
        """Runs q3dfit over the mode-selected spaxels and returns the
        (cols, rows) region for `q3dcollect`, or (None, None) when there is
        no rectangular region to collect.

        `initproc` is what q3dfit is invoked on: the in-memory `q3din` for a
        serial run, or the path to its saved `.npy` for a multi-core (MPI)
        run. `q3di` is still used for the grid dimensions. `nocrash` lets
        q3dfit skip and continue past spaxels whose fit raises.

        This cube indexes spaxels as science_data[:, x_pix, y_pix], so
        x_pix runs along the FITS second axis (q3dfit's row) and y_pix along
        the first (q3dfit's column); the +1 converts to q3dfit's unity
        offset.
        """
        if cols is not None or rows is not None:
            run_q3dfit(initproc, cols=cols, rows=rows, ncores=ncores,
                       quiet=quiet, nocrash=nocrash)
            return cols, rows
        if self.mode == "spaxel":
            x_pix, y_pix = self.test_spaxel
            col, row = y_pix + 1, x_pix + 1
            run_q3dfit(initproc, cols=col, rows=row, ncores=ncores,
                       quiet=quiet, nocrash=nocrash)
            return col, row
        if self.mode == "spaxels":
            for x_pix, y_pix in self.test_spaxels:
                run_q3dfit(initproc, cols=y_pix + 1, rows=x_pix + 1,
                           ncores=ncores, quiet=quiet, nocrash=nocrash)
            return None, None
        # cube: fit every spaxel, then collect over the full grid.
        run_q3dfit(initproc, cols=None, rows=None, ncores=ncores,
                   quiet=quiet, nocrash=nocrash)
        return [1, q3di.ncols], [1, q3di.nrows]

    # FITS output

    # Structural/geometry keywords of the source cube that describe its
    # own array layout and 3D (spectral) WCS -- they do not apply to the
    # 2D parameter maps, so they are dropped when copying provenance.
    _NONPROVENANCE_KEYWORDS = {
        "SIMPLE", "BITPIX", "EXTEND", "XTENSION", "PCOUNT", "GCOUNT",
        "EXTNAME", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "WCSAXES",
        "CRPIX1", "CRPIX2", "CRPIX3", "CRVAL1", "CRVAL2", "CRVAL3",
        "CDELT1", "CDELT2", "CDELT3", "CTYPE1", "CTYPE2", "CTYPE3",
        "CUNIT1", "CUNIT2", "CUNIT3", "PC1_1", "PC1_2", "PC1_3",
        "PC2_1", "PC2_2", "PC2_3", "PC3_1", "PC3_2", "PC3_3",
        "PC01_01", "PC02_02", "PC03_03", "BUNIT",
    }

    def _provenance_header(self):
        """Collects the provenance keywords of the original datacube
        (target, instrument, program, calibration, etc.) as a Header, so
        each output map can be traced back to the cube it was fitted from.

        The cube's own structural and 3D-WCS keywords are dropped (see
        `_NONPROVENANCE_KEYWORDS`); the correct 2D spatial WCS and per-map
        units are attached separately in `_map_extension_header`.
        """
        prov = fits.Header()
        for header in (self.general_header, self.science_header):
            for card in header.cards:
                key = card.keyword
                if key in ("", "COMMENT", "HISTORY"):
                    continue
                if key in self._NONPROVENANCE_KEYWORDS:
                    continue
                if key in prov:
                    continue
                prov[key] = (card.value, card.comment)
        return prov

    def _add_seed_keywords(self, header, cfg):
        """Records the key values used to seed the fit into `header`, so
        the initial conditions (line center, allowed center/width ranges,
        acceptance thresholds) travel with the results."""
        wv_unit = cfg["wv_unit"]
        header["LINENAME"] = (self.name, "Fitted emission line identifier")
        header["LINECEN"] = (cfg["line_center"], f"Seed line center [{wv_unit}]")
        header["CENOFF"] = (cfg["center_offset"], f"Max primary center offset [{wv_unit}]")
        header["CENCUT"] = (cfg["center_cutoff"], f"Max secondary center offset [{wv_unit}]")
        header["SIGMIN"] = (cfg["minimum_sigma"], f"Minimum Gaussian sigma [{wv_unit}]")
        header["SIGMAX"] = (cfg["maximum_sigma"], f"Maximum Gaussian sigma [{wv_unit}]")
        header["SIGNAR"] = (cfg["narrow_guess_sigma"], f"Narrow sigma seed [{wv_unit}]")
        header["SIGBRD"] = (cfg["broad_guess_sigma"], f"Broad sigma seed [{wv_unit}]")
        header["SNRTHR"] = (cfg["snr_threshold"], "Component amplitude/error SNR threshold")
        header["FALPHA"] = (cfg["alpha"], "F-test significance to prefer double model")
        header["WVUNIT"] = (wv_unit, "Wavelength unit of center/sigma/seed values")

    def _map_extension_header(self, cfg, prov, bunit, moment, comment, has_components):
        """Builds the header for one parameter-map extension: 2D spatial
        WCS from the cube + provenance + this array's unit and layout +
        the fit seed values."""
        header = self.wcs.to_header()
        header.extend(prov, update=False)
        header["BUNIT"] = (bunit, "Physical unit of the data array")
        header["MOMENT"] = (moment, comment)
        if has_components:
            # Axis 3 stacks the (narrow-to-broad ordered) components; it
            # carries no WCS, only the plane count.
            header["NCOMPMAX"] = (3, "Number of component planes along axis 3")
        header["EXTNAME"] = moment
        self._add_seed_keywords(header, cfg)
        return header

    def _total_flux_map(self):
        """Builds the per-spaxel total line flux (summed over all fitted
        components) in W/m^2; spaxels with no fit are NaN.

        The lmfit Gaussian ``amplitude`` is the integrated area under a
        component in (wavelength, surface brightness) space, i.e.
        MJy/sr * micron. The frequency-integrated surface brightness is
        therefore amplitude * c / lambda^2, and multiplying by the spaxel
        solid angle gives the component's line flux.
        """
        base_array = np.full(self.im_shape, np.nan)
        sb_area_unit = u.MJy / u.steradian * u.micron
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            if np.isnan(self.fitparams["G1AMP"][idx]):
                continue
            flux_val = 0.0
            for comp in range(1, 4):
                amp = self.fitparams[f"G{comp}AMP"][idx]
                cen = self.fitparams[f"G{comp}CEN"][idx]
                if np.isnan(amp):
                    break
                comp_flux = (amp * sb_area_unit * const.c / (cen * u.micron) ** 2
                             * self.area_sr).to(u.watt / u.meter ** 2)
                flux_val += comp_flux.value
            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = max(flux_val, 0.0)
        return base_array

    def _order_fitparams_by_outflow(self, line_center):
        """Reorders each spaxel's fitted components in `self.fitparams` by
        increasing distance of the fitted center from the line center, so
        the systemic component occupies the G1 columns and the outflow --
        by convention the component offset furthest from the line center
        -- the G2 columns.

        This replaces, for the output products, the narrow-to-broad
        (sigma) ordering imposed by `_validate_components`: velocity
        offset, not width, is the outflow discriminant we want the
        G1/G2/G3 columns to encode. Applied once to the assembled table so
        the IPAC table, the rendered maps, and the FITS cube all share the
        same component identity. Ties keep their incoming (sigma) order;
        rows with fewer than two fitted components are left unchanged.
        """
        amp_cols = ["G1AMP", "G2AMP", "G3AMP"]
        cen_cols = ["G1CEN", "G2CEN", "G3CEN"]
        sig_cols = ["G1SIGMA", "G2SIGMA", "G3SIGMA"]
        for row_idx in range(len(self.fitparams)):
            centers = np.array([self.fitparams[col][row_idx] for col in cen_cols])
            valid = np.where(np.isfinite(centers))[0]
            if len(valid) < 2:
                continue
            order = valid[np.argsort(np.abs(centers[valid] - line_center),
                                     kind="stable")]
            if np.array_equal(order, valid):
                continue
            amps = [self.fitparams[amp_cols[i]][row_idx] for i in order]
            cens = [self.fitparams[cen_cols[i]][row_idx] for i in order]
            sigs = [self.fitparams[sig_cols[i]][row_idx] for i in order]
            for slot, dest in enumerate(valid):
                self.fitparams[amp_cols[dest]][row_idx] = amps[slot]
                self.fitparams[cen_cols[dest]][row_idx] = cens[slot]
                self.fitparams[sig_cols[dest]][row_idx] = sigs[slot]

    def write_fits(self, filepath=None, cfg=None):
        """Writes the fitted parameter maps to a multi-extension FITS file.

        Five data arrays are stored, one per fitted quantity:

        * ``AMPLITUDE`` -- the Gaussian amplitude (integrated line area,
          i.e. the 0th moment) of each component;
        * ``CENTER``    -- the line center (1st moment);
        * ``SIGMA``     -- the velocity dispersion (2nd moment);
        * ``TOTFLUX``   -- the total line flux summed over components;
        * ``CHISQR``    -- the reduced chi-square of the preferred fit.

        The three moment arrays have shape ``(3, XPIX, YPIX)``: plane k
        holds the k-th component's map, ordered by increasing distance
        from the line center (G1 systemic, G2 outflow; see
        `_order_fitparams_by_outflow`), NaN where that component was not
        fitted. ``TOTFLUX`` and ``CHISQR``
        are single ``(XPIX, YPIX)`` maps. Every extension header inherits the source
        cube's provenance keywords and 2D spatial WCS, carries the
        physical unit of its array (``BUNIT``) and its layout, and records
        the values used to seed the fit (line center, allowed center
        offsets, sigma bounds; see `_add_seed_keywords`).

        Args:
            filepath (str, optional): Output path. Defaults to
                ``<output>twogaussian_fit.fits``.
            cfg (dict, optional): The fit config (as built by
                ``_build_fit_config``) whose seed values are recorded in the
                headers. Rebuilt if not supplied.
        """
        if filepath is None:
            filepath = self.output + "twogaussian_fit.fits"
        if cfg is None:
            cfg = self._build_fit_config()

        nx, ny = self.im_shape
        amplitude = np.full((3, nx, ny), np.nan)
        center = np.full((3, nx, ny), np.nan)
        sigma = np.full((3, nx, ny), np.nan)
        redchi = np.full((nx, ny), np.nan)

        x_pixels = np.asarray(self.fitparams["XPIX"], dtype=int)
        y_pixels = np.asarray(self.fitparams["YPIX"], dtype=int)
        for row_idx in range(len(x_pixels)):
            x_pix, y_pix = x_pixels[row_idx], y_pixels[row_idx]
            for comp in range(3):
                amplitude[comp, x_pix, y_pix] = self.fitparams[f"G{comp+1}AMP"][row_idx]
                center[comp, x_pix, y_pix] = self.fitparams[f"G{comp+1}CEN"][row_idx]
                sigma[comp, x_pix, y_pix] = self.fitparams[f"G{comp+1}SIGMA"][row_idx]
            redchi[x_pix, y_pix] = self.fitparams["REDCHI"][row_idx]
        total_flux = self._total_flux_map()

        amp_unit = (self._sb_unit * self._wv_unit).to_string("fits")
        wv_unit = self._wv_unit.to_string("fits")
        flux_unit = (u.watt / u.meter ** 2).to_string("fits")

        prov = self._provenance_header()

        primary = fits.PrimaryHDU(header=prov)
        primary.header["LINENAME"] = (self.name, "Fitted emission line identifier")
        primary.header.add_history("Multicomponent Gaussian spaxel fit (irspec.Spaxelcube).")
        primary.header.add_history("Maps indexed [component, XPIX, YPIX]; components "
                                   "ordered by distance from line center (G1 systemic, "
                                   "G2 outflow).")
        primary.header.add_history(f"Source cube: {os.path.basename(str(self.filepath))}")

        hdus = [primary]
        for data, moment, bunit, comment, has_components in (
            (amplitude, "AMPLITUDE", amp_unit,
             "Gaussian amplitude / line area (0th moment)", True),
            (center, "CENTER", wv_unit, "Line center (1st moment)", True),
            (sigma, "SIGMA", wv_unit, "Velocity dispersion sigma (2nd moment)", True),
            (total_flux, "TOTFLUX", flux_unit, "Total line flux summed over components", False),
            (redchi, "CHISQR", "", "Reduced chi-square of the preferred fit", False),
        ):
            header = self._map_extension_header(cfg, prov, bunit, moment,
                                                comment, has_components)
            hdus.append(fits.ImageHDU(data=data.astype(np.float32),
                                      header=header, name=moment))

        fits.HDUList(hdus).writeto(filepath, overwrite=True)
        return filepath


    def validate_snr(self, best_result):
        """Ensures that fitted components pass the SNR threshold and returns
        the array of fitted values.

        SNR here is each component's fitted line flux (the lmfit Gaussian
        `amplitude` parameter, i.e. the integrated area under the curve)
        divided by its fit uncertainty. This is only meaningful because
        the fit itself is weighted by the real per-pixel flux errors
        (see `two_gaussian_fit`) -- lmfit then derives `amplitude.stderr`
        from the resulting covariance matrix, rather than us having to
        estimate detection significance from the noiseless best-fit curve.

        Surviving components are sorted narrow-to-broad (ascending sigma)
        so the G1/G2/G3 output columns keep a consistent physical identity
        across spaxels.
        """
        return _validate_components(best_result)
    
    
    def reconstruct_gaussian_fluxes(self, best_result, spaxel_line, background_flux, line_wavelengths):
        """Regenerates the continuum flux, the total flux, and all of the
        gaussian component fluxes for an lmfit gaussian decomposition model."""
        
        prefixes = ["g1_", "g2_", "g3_"]
        model_components = []
        total_flux = np.copy(background_flux)
        for prefix in prefixes:
            if prefix+"amplitude" in best_result.best_values:
                params = Parameters()
                params.add(name="amplitude", 
                           value=best_result.best_values[prefix + "amplitude"])
                params.add(name="center", 
                           value=best_result.best_values[prefix + "center"])
                params.add(name="sigma", 
                           value=best_result.best_values[prefix + "sigma"])
                model_component = GaussianModel()
                component_flux = model_component.eval(params=params, x=line_wavelengths)
                total_flux += component_flux
                model_components.append(background_flux + component_flux)
        return total_flux, model_components
    
    
    def reconstruct_spaxel_fluxes(self, x_pix, y_pix):
        """Regenerates the continuum flux, the total flux, and all of the
        gaussian component fluxes for a given spaxel using the completed 
        fitting parameter dictionary."""
        (flux, flux_err, dq) = self.spaxel_values(y_pix, x_pix)
        spaxel = Spectrum1D(flux, self.wvs)
        spaxel_continuum, spaxel_line = self._mask_spaxel(spaxel)
        line_wavelengths, line_fluxes, background_flux, _ = self.background_subtraction(spaxel_continuum, spaxel_line)

        # Generate component fluxes and total flux
        model_components = []
        total_flux = np.copy(background_flux)
        spaxel_idx = x_pix * self.im_shape[0] + y_pix
        for idx in range(1, int(self.fitparams["NCOMP"][spaxel_idx]) + 1):
            params = Parameters()
            params.add(name="amplitude", value=self.fitparams[f"G{idx}AMP"][spaxel_idx])
            params.add(name="center", value=self.fitparams[f"G{idx}CEN"][spaxel_idx])
            params.add(name="sigma", value=self.fitparams[f"G{idx}SIGMA"][spaxel_idx])
            model_component = GaussianModel()
            component_flux = model_component.eval(params=params, x=line_wavelengths)
            total_flux += component_flux
            model_components.append(background_flux + component_flux)
        return spaxel_line, background_flux, total_flux, model_components
    
    
    def render_spaxel_fit(self, x_pix, y_pix, file_stem=None, savefig=False):
        """Renders the multicomponent gaussian fit for an individual spaxel."""
        spaxel_line, background_flux, total_flux, model_components = self.reconstruct_spaxel_fluxes(x_pix, y_pix)

        # Plot
        self.pltparams.apply()
        # The data and continuum are drawn in the palette's foreground: a
        # hardcoded black spectrum is invisible on the dark background.
        data_color = self.pltparams.foreground
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)
        ax.scatter(spaxel_line.spectral_axis.value, spaxel_line.flux.value, label="Data", c=data_color)
        ax.plot(spaxel_line.spectral_axis.value, background_flux, label="Continuum", c=data_color, ls="dashed")
        ax.plot(spaxel_line.spectral_axis.value, total_flux, label="Model", c="gray")
        #for idx, component in enumerate(model_components):
        #    ax.plot(spaxel_line.spectral_axis.value, component, label=f"Component {idx+1}")
        ax.set_title(f"Spaxel: ({x_pix}, {y_pix})", loc="right")
        ax.set_title(f"Gaussian Decomposition", loc="left")
        ax.set_xlabel("Wavelength [μm]")
        ax.set_ylabel("Surface Brightness [MJy/sr]")
        ax.legend()
        if savefig:
            filename = (f"{self.name}_spaxel_{x_pix}_{y_pix}"
                        f"{self.pltparams.file_suffix}.png")
            plt.savefig(self._plot_path(file_stem, filename), dpi=300)
            plt.close()
        else:
            plt.show()

    # Diagnostic value display methods
    
    def display_masked_spaxel(self, spaxel_continuum, spaxel_line):
        
        fig = plt.figure()
        fig.set_size_inches(10, 8)
        ax = plt.subplot()
        ax.plot(spaxel_continuum.spectral_axis.value, spaxel_continuum.flux.value, label="Continuum")
        ax.plot(spaxel_line.spectral_axis.value, spaxel_line.flux.value, label="Line")
        ax.set_xlabel("Wavelength [μm]")
        ax.set_ylabel("Flux [MJy/sr]")
        ax.legend()
        plt.show()
    
    def display_background_subtraction(self, line_wavelengths, line_fluxes, background_flux):
        fig = plt.figure()
        fig.set_size_inches(10, 8)
        ax = plt.subplot()
        ax.plot(line_wavelengths, line_fluxes, label="Line Flux")
        ax.plot(line_wavelengths, background_flux, label="Background Flux")
        ax.set_xlabel("Wavelength [μm]")
        ax.set_ylabel("Flux [MJy/sr]")
        ax.legend()
        plt.show()
    
    def display_fitted_models(self, spaxel_line, background_flux, line_wavelengths, single_result, double_result):
        
        single_total_flux, single_model_components = self.reconstruct_gaussian_fluxes(single_result, 
                                                                        spaxel_line, 
                                                                        background_flux, 
                                                                        line_wavelengths)
        double_total_flux, double_model_components = self.reconstruct_gaussian_fluxes(double_result, 
                                                                        spaxel_line, 
                                                                        background_flux, 
                                                                        line_wavelengths)
        self.pltparams.apply()
        fig = plt.figure()
        fig.set_size_inches(10, 8)
        ax = plt.subplot()
        vels = self.wv_to_vel(line_wavelengths, self.line_center)
        ax.plot(vels, spaxel_line.flux.value, color=self.pltparams.foreground,
                lw=2, label="Data")
        ax.plot(vels, single_total_flux, color="blue", ls="dotted", lw=2, label="Single Gaussian Model")
        ax.plot(vels, double_total_flux, color="red", ls="dotted", lw=2, label="Double Gaussian Model")
        ax.set_xlabel("Velocity Offset [km/s]")
        ax.set_ylabel("Flux [MJy/sr]")
        ax.legend()
        plt.show()
    
    # Fitted value display methods
    
    # Map rendering helpers

    def _line_title(self):
        """Returns the display title for the fitted line."""
        titles = {"[H_2_S_1]": r"H$_{2}$ S(1)",
                  "[H_2_S_2]": r"H$_{2}$ S(2)",
                  "[H_2_S_3]": r"H$_{2}$ S(3)",
                  "[H_2_S_5]": r"H$_{2}$ S(5)",
                  "[NeV]_14": "[NeV]" + r"$_{14}$"}
        return titles.get(self.name, self.name)

    def _setup_map_axes(self, use_wcs):
        """Creates a spaxel-map figure and axes: a WCS sky projection with
        a physical scalebar (from `self.distance`), or plain pixel axes.

        Note on the pixel axes: map arrays are indexed [XPIX][YPIX], and
        imshow displays the first (row) index vertically, so XPIX is the
        *vertical* axis. This matches the orientation of the WCS maps
        (XPIX runs along declination), keeping every map in the same
        frame regardless of projection.

        This object's palette is re-asserted first: rcParams are global, so
        without it the figure inherits whatever palette another module (or
        another Spaxelcube) applied most recently.
        """
        self.pltparams.apply()
        fig = plt.figure()
        fig.set_size_inches(10, 8)
        ax = (plt.subplot(projection=self.wcs) if use_wcs else plt.subplot())
        self._decorate_map_axes(ax, use_wcs)
        return fig, ax

    def _decorate_map_axes(self, ax, use_wcs, scalebar=True, axis_labels=True):
        """Applies this module's standard map furniture to an existing axis:
        coordinate labels, the physical scalebar and compass, the framing
        margin, and the footprint or contour overlay.

        Split out of `_setup_map_axes` so a multi-panel figure can dress each
        of its panels the same way a standalone map is dressed, rather than
        keeping a second, drifting copy of the same decoration. Panels that
        share a row or column suppress the scalebar and the axis labels, which
        would otherwise be repeated six times over.
        """
        if use_wcs:
            if scalebar:
                scalebar_angle = (1 * u.kpc / self.distance).to(
                    u.deg, equivalencies=u.dimensionless_angles())
                # Follow the active palette and scaling: a hardcoded black
                # scalebar is invisible on the dark background style, and a
                # hardcoded size leaves 24pt text on a 14pt paper figure.
                fontprops = fm.FontProperties(size=self.pltparams.annotation_size,
                                              family=self.pltparams.font)
                add_scalebar(ax, scalebar_angle, label="1 kpc",
                             color=self.pltparams.foreground,
                             fontproperties=fontprops)
                self._add_compass(ax)
            if axis_labels:
                ax.set_xlabel("Right Ascension")
                ax.set_ylabel("Declination")
            else:
                # WCSAxes labels its coordinate axes from the FITS ctype
                # unless told otherwise, so simply not setting a label leaves
                # "pos.eq.ra"/"pos.eq.dec" on the panel. Blank them.
                ax.coords[0].set_axislabel("")
                ax.coords[1].set_axislabel("")
        elif axis_labels:
            ax.set_xlabel("YPIX")
            ax.set_ylabel("XPIX")
        self._apply_map_limits(ax)
        # The footprint outline and the imaging contours are alternatives,
        # not additions: both are reference geometry drawn over the data, and
        # together they clutter the map. The contours need the sky projection
        # to register against, so pixel axes always get the outline.
        if self.contour_image is not None and use_wcs:
            self._draw_image_contours(ax)
        else:
            self._draw_footprint(ax)
        return ax

    # Blank margin around the map, in spaxels.
    _MAP_MARGIN = 2

    def _apply_map_limits(self, ax):
        """Frames the map with a margin of empty sky around the cube grid.

        Most of these cubes are cut off by their own spatial grid -- the
        footprint runs to the first or last row or column -- so without a
        margin the outermost measured spaxels sit flush against the axes
        spine, which covers them and their share of the footprint outline.

        Set here rather than after the image because `set_xlim` turns
        autoscaling off, so the later `imshow` inherits these limits instead
        of snapping back to the array bounds.
        """
        rows, cols = self.im_shape
        margin = self._MAP_MARGIN
        ax.set_xlim(-0.5 - margin, cols - 0.5 + margin)
        ax.set_ylim(-0.5 - margin, rows - 0.5 + margin)

    def footprint_mask(self):
        """Boolean map of the spaxels the cube actually measured.

        This is the instrument footprint, not the fitted region: it covers
        every spaxel that carries a spectrum, including those where no fit
        converged. The distinction is the point of drawing it -- a map of
        fitted spaxels alone cannot tell "the line is absent here" from
        "nothing was observed here".

        A single finite channel is enough to count. That is deliberately the
        loosest possible criterion, and it is what makes the outline honest:
        fitting a spaxel requires finite data in the line window, so every
        fitted spaxel necessarily lands inside, with no threshold to tune. A
        stricter cut (half the channels finite, say) looks tidier but cuts
        off the ragged edge of the dither pattern, leaving a scatter of
        fitted spaxels stranded outside their own footprint.

        The criterion is the whole cube, not the line's fit window, so every
        map of a given cube carries the same outline.
        """
        return np.any(np.isfinite(self.science_data), axis=0)

    def _draw_footprint(self, ax):
        """Outlines the cube footprint on a map.

        Called from `_setup_map_axes`, before the map image exists. That is
        safe because the contour carries a higher zorder than an image (which
        defaults to 0), so it still draws on top -- and it keeps the call in
        one place rather than repeated in every render_* method.
        """
        segments = self._footprint_edges(self.footprint_mask())
        if len(segments) == 0:
            return
        # autolim=False: the view is already fixed by `_apply_map_limits`, and
        # the outline should not be able to argue with it.
        ax.add_collection(LineCollection(
            segments, colors=self.pltparams.foreground,
            linewidths=self.pltparams.annotation_size / 12, zorder=3),
            autolim=False)

    # How much wider than the cube grid the contour cutout is taken, so
    # contours run to the edge of the map instead of stopping inside it.
    _CONTOUR_CUTOUT_MARGIN = 1.4

    def _contour_cutout(self):
        """The patch of the contour image covering this cube, and its WCS.

        Cut down from the full mosaic because the mosaic is tens of
        arcseconds across and a cube covers a few: contouring it whole would
        set the levels from sky that is never shown, and spend its time on
        vertices outside the map.
        """
        _, path = self.contour_image
        image, image_wcs = load_contour_image(path)

        rows, cols = self.im_shape
        corners = self.wcs.pixel_to_world([-0.5, cols - 0.5, cols - 0.5, -0.5],
                                          [-0.5, -0.5, rows - 0.5, rows - 0.5])
        x, y = image_wcs.world_to_pixel(corners)
        size = self._CONTOUR_CUTOUT_MARGIN * max(x.max() - x.min(),
                                                 y.max() - y.min())
        centre = (0.5 * (x.min() + x.max()), 0.5 * (y.min() + y.max()))
        # mode="partial" so a cube near the edge of the mosaic still yields a
        # cutout, padded with NaN, rather than raising.
        cutout = Cutout2D(image, centre, int(round(size)), wcs=image_wcs,
                          mode="partial", copy=True)
        return cutout.data, cutout.wcs

    def _draw_image_contours(self, ax):
        """Overlays the contour image's intensity contours on a map.

        Drawn through `ax.get_transform`, which converts from the image's own
        pixel grid to the map's. That keeps the imaging at full resolution --
        NIRCam's pixels are four to ten times finer than an MRS spaxel, so
        resampling it onto the cube grid first would throw away exactly the
        detail the overlay is there to show.
        """
        data, image_wcs = self._contour_cutout()
        # `contour_field` closes the level sets over the saturated nucleus
        # and caps the levels below it; without that the innermost contours
        # ring the hole the saturation leaves and appear to place the source
        # ~0.2" from where it is.
        field, levels = contour_field(data)
        if field is None:
            return
        # A single foreground-coloured set rather than a colormapped one: the
        # map underneath already spends a full colormap on the fitted
        # quantity, and a second one competing with it reads as data.
        ax.contour(field, levels=levels,
                   colors=self.pltparams.foreground,
                   linewidths=self.pltparams.annotation_size / 20,
                   alpha=0.9, transform=ax.get_transform(image_wcs), zorder=3)

    @staticmethod
    def _footprint_edges(mask):
        """The outline of `mask`, as segments running along spaxel edges.

        `contour` interpolates between cell *centres*, so on a staircase
        boundary -- which is what a rotated footprint on a square grid always
        is -- it cuts the corners, drawing a smooth diagonal straight through
        the spaxels it is meant to be enclosing. Walking the cell edges puts
        the outline exactly where the measurements stop.

        A spaxel contributes an edge wherever its neighbour is unmeasured,
        counting everything outside the array as unmeasured. That closes the
        outline along the edge of the cube grid, where most of these cubes cut
        their footprint off, and around any gap inside it -- both for free,
        with no padding to arrange.
        """
        if not mask.any():
            return []
        padded = np.pad(mask, 1, constant_values=False)
        # Each entry: the neighbours in one direction, and the two corners of
        # the shared edge, as offsets from the spaxel's centre.
        sides = (
            (padded[1:-1, :-2], (-0.5, -0.5), (-0.5, +0.5)),   # left
            (padded[1:-1, 2:],  (+0.5, -0.5), (+0.5, +0.5)),   # right
            (padded[:-2, 1:-1], (-0.5, -0.5), (+0.5, -0.5)),   # below
            (padded[2:, 1:-1],  (-0.5, +0.5), (+0.5, +0.5)),   # above
        )
        segments = []
        for neighbour, (dx0, dy0), (dx1, dy1) in sides:
            rows, cols = np.nonzero(mask & ~neighbour)
            if len(rows) == 0:
                continue
            segments.extend(np.stack(
                [np.column_stack((cols + dx0, rows + dy0)),
                 np.column_stack((cols + dx1, rows + dy1))], axis=1))
        return segments

    # Compass anchor, in axes fractions, and arrow length as a fraction of
    # the shorter pixel axis. Top left: the MRS and NIRSpec footprints are
    # rotated rectangles, so their corners are empty sky, and the scalebar
    # already holds the bottom right.
    #
    # The anchor is inset far enough that an arrow plus its letter stays
    # inside the frame whichever way the WCS rotation points them: the reach
    # is length x _COMPASS_LABEL_OFFSET in any direction, so the anchor needs
    # at least that much clearance on both sides.
    _COMPASS_ANCHOR = (0.18, 0.80)
    _COMPASS_LENGTH = 0.09
    _COMPASS_LABEL_OFFSET = 1.3

    def _add_compass(self, ax):
        """Draws a north/east compass on a WCS map.

        The directions are read off the WCS rather than assumed: a cube is
        rotated on the sky by the position angle of its observation, so north
        is not up and east is not left. Each arrow is built by stepping a
        small angular distance from the anchor along a position angle -- 0 for
        north, 90 for east, since position angle runs east of north -- and
        projecting the result back into pixel coordinates.

        East is the direction right ascension increases, which on a sky image
        points opposite the way it would in a plain cartesian plot; drawing it
        from the WCS gets the handedness right for free.
        """
        rows, cols = self.im_shape
        x_anchor = self._COMPASS_ANCHOR[0] * cols
        y_anchor = self._COMPASS_ANCHOR[1] * rows
        anchor = self.wcs.pixel_to_world(x_anchor, y_anchor)

        length = self._COMPASS_LENGTH * min(rows, cols)
        color = self.pltparams.foreground
        # A halo in the background colour. The MRS footprints leave their
        # corners empty, but the NIRSpec cube fills its whole grid, so there
        # the compass has no choice but to sit on data; the halo keeps it
        # readable over a colormap without a box hiding the spaxels under it.
        halo = [pe.withStroke(linewidth=3, foreground=self.pltparams.background)]
        for position_angle, label in ((0 * u.deg, "N"), (90 * u.deg, "E")):
            # The step only has to fix a direction, not a distance: it is
            # normalised away below, so it just needs to be small enough to
            # stay in the local tangent plane.
            tip = self.wcs.world_to_pixel(
                anchor.directional_offset_by(position_angle, 1 * u.arcsec))
            dx, dy = float(tip[0]) - x_anchor, float(tip[1]) - y_anchor
            norm = np.hypot(dx, dy)
            if norm == 0:
                continue
            dx, dy = dx / norm * length, dy / norm * length
            ax.annotate("", xy=(x_anchor + dx, y_anchor + dy),
                        xytext=(x_anchor, y_anchor),
                        arrowprops=dict(arrowstyle="-|>", color=color,
                                        linewidth=1.5, shrinkA=0, shrinkB=0,
                                        path_effects=halo),
                        zorder=4)
            ax.text(x_anchor + dx * self._COMPASS_LABEL_OFFSET,
                    y_anchor + dy * self._COMPASS_LABEL_OFFSET, label,
                    color=color, fontsize=self.pltparams.annotation_size,
                    ha="center", va="center", path_effects=halo, zorder=4)

    def _plot_path(self, file_stem, filename):
        """Resolves where a rendered figure is written.

        `file_stem` overrides the destination for a single call; when it is
        None the figure goes to `self.plot_output` (the per-run directory
        under `plots/line_maps/<line>/`). The directory is created on demand,
        so callers never have to pre-make it.
        """
        directory = self.plot_output if file_stem is None else file_stem
        os.makedirs(directory, exist_ok=True)
        return os.path.join(directory, filename)

    def _finish_map(self, suffix, file_stem, savefig):
        """Saves the current map figure to disk or displays it.

        The palette tag (``_light``/``_dark``) closes the filename, so the
        two palettes of the same map are separate files in the run directory
        instead of one overwriting the other. A contour overlay adds its own
        tag ahead of that, for the same reason: the overlaid variation of a
        map and the footprint version are different figures of the same fit.
        """
        if savefig:
            overlay = ("" if self.contour_image is None
                       else f"_{self.contour_image[0]}")
            filename = (f"{self.name}_{suffix}{overlay}"
                        f"{self.pltparams.file_suffix}.png")
            plt.savefig(self._plot_path(file_stem, filename), dpi=300)
            plt.close()
        else:
            plt.show()

    # Fraction of a fit bound at or above which a value counts as pinned to
    # that bound rather than measured.
    _PINNED_FRACTION = 0.98

    def _unpinned(self, values, bound=None):
        """The finite entries of `values` that are not sitting on a fit bound.

        lmfit returns the bound itself when a parameter walks into its limit,
        so those spaxels record where the fit stopped, not what was measured.
        They stay on the map -- saturated at the end of the color bar -- but
        they must not set the color limits: a percentile taken over them
        stretches the scale across the fitter's entire allowed range and
        flattens the real structure into one flat tone.

        Falls back to the full finite set when no bound is known (the q3dfit
        maps, which carry their own bounds) or when every spaxel is pinned.
        """
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if bound is None or len(finite) == 0:
            return finite
        unpinned = finite[np.abs(finite) < self._PINNED_FRACTION * bound]
        return unpinned if len(unpinned) > 0 else finite

    @staticmethod
    def _robust_spread(values):
        """Half the 16th-84th percentile range: the 1-sigma width of the bulk
        of a distribution, which a heavy tail cannot inflate the way it
        inflates a standard deviation."""
        return 0.5 * (np.percentile(values, 84) - np.percentile(values, 16))

    def _draw_component_map(self, ax, base_array):
        """Draws a component-count map: one discrete color per integer count.

        `viridis` rather than a grayscale, because a grayscale's endpoints are
        black and white -- exactly the two palettes' background colors -- so
        the highest count would read as a hole in the map under the light
        palette, and the lowest under the dark one.
        """
        finite = base_array[np.isfinite(base_array)]
        if len(finite) == 0:
            low, high = 0, 1
        else:
            low, high = int(np.min(finite)), int(np.max(finite))
        cmap = self.pltparams.colormap("viridis", high - low + 1)
        image = ax.imshow(base_array, cmap=cmap, vmin=low - 0.5, vmax=high + 0.5,
                          origin="lower")
        plt.colorbar(image, ticks=np.arange(low, high + 1))
        ax.set_title("# of Components", loc="right")
        ax.set_title(self._line_title(), loc="left")

    def render_multicomponent_plot(self, file_stem=None, use_wcs=True, savefig=False):
        """Renders a spaxel map of the number of fitted components.

        Unfitted spaxels stay NaN so the palette paints them as empty sky,
        the way every other map treats them. Filling them with 0 instead
        invents a "fitted, found nothing" value the fit table never contains
        -- and, being the bottom of the color scale, buries the scalebar in a
        block of solid color.
        """
        base_array = np.full(self.im_shape, np.nan)
        ncomp = np.asarray(self.fitparams["NCOMP"], dtype=float)
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            if not np.isnan(ncomp[idx]):
                base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = ncomp[idx]

        fig, ax = self._setup_map_axes(use_wcs)
        self._draw_component_map(ax, base_array)
        self._finish_map("components", file_stem, savefig)

    def render_totflux_plot(self, file_stem=None, use_wcs=True, savefig=False):
        """Renders a spaxel map of the total fitted line flux in W/m^2
        (summed over components; see `_total_flux_map`)."""
        base_array = self._total_flux_map()

        fig, ax = self._setup_map_axes(use_wcs)
        # Percentile-clipped log stretch: an auto-ranged LogNorm lets a
        # few near-zero spaxels stretch the scale by decades and wash out
        # the bright structure.
        positive = base_array[np.isfinite(base_array) & (base_array > 0)]
        if len(positive) > 0 and np.percentile(positive, 2) < np.percentile(positive, 99.5):
            norm = LogNorm(vmin=np.percentile(positive, 2),
                           vmax=np.percentile(positive, 99.5))
        else:
            norm = None
        image = ax.imshow(base_array, cmap=self.pltparams.colormap("plasma"),
                          norm=norm, origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, r"[W/m$^2$]")
        ax.set_title("Total Flux", loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map("total_flux", file_stem, savefig)

    def render_rel_vel_plot(self, param="G2CEN", file_stem=None, use_wcs=True, savefig=False):
        """Renders a spaxel map of a fitted component's velocity offset
        from the systemic line center."""
        base_array = np.full(self.im_shape, np.nan)
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            rel_vel = self.wv_to_vel(self.fitparams[param][idx],
                                              self.line_dict[self.name][4])
            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = rel_vel.value

        # Symmetric color limits, so the diverging colormap stays centered on
        # systemic and neither sign can saturate the other. The 98th
        # percentile of |v| alone is not enough: spaxels pinned at
        # `center_cutoff` are dropped first (see `_unpinned`), and the result
        # is capped at three robust sigma. A narrow component typically has a
        # body of a few tens of km/s over a tail reaching the cutoff, and the
        # percentile hands almost the whole color range to that tail --
        # rendering the rotation field as a flat sheet. The cap can only
        # tighten the stretch, so a genuinely broad component keeps its
        # percentile limits.
        vels = self._unpinned(base_array, self.line_dict[self.name][6])
        vel_limit = 1.0
        if len(vels) > 0:
            vel_limit = min(np.percentile(np.abs(vels), 98),
                            3 * self._robust_spread(vels))
        if not vel_limit > 0:
            vel_limit = 1.0

        comp_name = r"$v_{" + param[1] + r"}$"

        fig, ax = self._setup_map_axes(use_wcs)
        image = ax.imshow(base_array, vmin=-vel_limit, vmax=vel_limit,
                          cmap=self.pltparams.colormap("bwr"), origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, "[km/s]")
        ax.set_title(comp_name, loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map(param, file_stem, savefig)

    def render_sigma_plot(self, param="G2SIGMA", file_stem=None, use_wcs=True,
                          savefig=False, subtract_lsf=True):
        """Renders a spaxel map of a fitted component's velocity dispersion.

        By default the instrumental line width (sigma = lambda / R /
        2.3548, from the MRS resolving power) is subtracted in quadrature
        so the map shows the intrinsic gas dispersion rather than the
        instrument; spaxels at or below the instrumental width map to 0.
        Pass `subtract_lsf=False` for the observed dispersion.
        """
        line_center = self.line_dict[self.name][4]
        sigma_inst_vel = 0.0
        if subtract_lsf:
            resolving_power = self.mrs_resolving_power(line_center)
            sigma_inst_um = line_center / resolving_power / 2.3548
            sigma_inst_vel = ((const.c * sigma_inst_um / line_center)
                              .to(u.kilometer / u.second)).value

        base_array = np.full(self.im_shape, np.nan)
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            vel_disp = ((const.c * np.abs(self.fitparams[param][idx]) / line_center)
                        .to(u.kilometer / u.second)).value
            if subtract_lsf and np.isfinite(vel_disp):
                vel_disp = np.sqrt(max(vel_disp ** 2 - sigma_inst_vel ** 2, 0.0))
            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = vel_disp

        # Color limits from the plotted component only: pooling every
        # component's range compresses the narrow component into the
        # bottom of the color scale.
        #
        # Spaxels pinned at `maximum_sigma` are dropped before the
        # percentiles, and the upper limit is capped at three robust sigma
        # above the median, for the reason given in `render_rel_vel_plot`:
        # the narrow component's dispersion piles up near the instrumental
        # width, so the handful of spaxels sitting on the bound would
        # otherwise own most of the color range. The bound is LSF-corrected
        # to match the map when `subtract_lsf` is on.
        max_sigma = self.line_dict[self.name][10]
        if subtract_lsf:
            max_sigma = np.sqrt(max(max_sigma ** 2 - sigma_inst_vel ** 2, 0.0))
        disps = self._unpinned(base_array, max_sigma)
        if len(disps) > 0:
            min_disp = np.percentile(disps, 1)
            max_disp = min(np.percentile(disps, 99),
                           np.median(disps) + 3 * self._robust_spread(disps))
            if max_disp <= min_disp:
                max_disp = min_disp + 1.0
        else:
            min_disp, max_disp = None, None

        comp_name = r"$\sigma_{" + param[1] + r"}$"

        fig, ax = self._setup_map_axes(use_wcs)
        image = ax.imshow(base_array, vmin=min_disp, vmax=max_disp,
                          cmap=self.pltparams.colormap("gist_heat"), origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, "[km/s]")
        ax.set_title(comp_name, loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map(param, file_stem, savefig)

    # The three moments of a line profile, in the sense used by spectral-cube
    # (https://spectral-cube.readthedocs.io/en/latest/moments.html):
    #
    #     M0 = INT I dv                      integrated intensity
    #     M1 = INT v I dv / M0               intensity-weighted velocity
    #     M2 = INT (v - M1)^2 I dv / M0      variance, i.e. velocity width^2
    #
    # Note M2 is the variance, so its unit is the square of the spectral axis
    # unit; spectral-cube's `linewidth_sigma` is its square root.
    _MOMENT_TITLES = (r"$M_0$", r"$M_1$", r"$M_2$")
    _MOMENT_UNITS = (r"[MJy sr$^{-1}$ km s$^{-1}$]", "[km/s]",
                     r"[km$^2$ s$^{-2}$]")
    # Component roles, top row first, as drawn by `render_moment_panels`.
    _MOMENT_ROWS = ("broad", "narrow")

    def _component_moment_maps(self, subtract_lsf=False):
        """Builds the 0th, 1st and 2nd moment maps of the narrow and the broad
        component of every spaxel's fit.

        Each fitted component is a single Gaussian, so the integrals above have
        closed forms and are evaluated exactly rather than summed over a
        sampled profile. For a component of lmfit amplitude ``A`` (the area
        under the profile in MJy/sr * micron), center ``lambda_c`` and width
        ``sigma_lambda``, converting the integration variable from wavelength
        to velocity (``dv = c dlambda / lambda_0``) gives

            M0 = A c / lambda_0,   M1 = c (lambda_c - lambda_0) / lambda_0,
            M2 = (c sigma_lambda / lambda_0)^2

        Narrow and broad are decided per spaxel by width, taking the smallest
        and largest ``|sigma|`` of the components that survived validation.
        They are *not* the G1/G2 columns: `_order_fitparams_by_outflow` sorts
        those by velocity offset, so G2 is the outflow component, which is
        frequently but not always the broader one. A spaxel fitted with a
        single component contributes to the narrow maps only, leaving the
        broad row NaN there -- calling a lone component "broad" would assert a
        decomposition the fit never made.

        Args:
            subtract_lsf (bool, optional): Remove the instrumental width from
                M2 in quadrature. Off by default, so M2 is the moment as
                defined above -- the observed variance. Note this is the
                opposite default from `render_sigma_plot`, which reports the
                intrinsic dispersion.

        Returns:
            dict: ``{'narrow': [M0, M1, M2], 'broad': [M0, M1, M2]}``, each an
            ``im_shape`` array that is NaN where that component is absent.
        """
        line_center = self.line_dict[self.name][4]
        sigma_inst_vel = 0.0
        if subtract_lsf:
            sigma_inst_um = (line_center / self.mrs_resolving_power(line_center)
                             / 2.3548)
            sigma_inst_vel = ((const.c * sigma_inst_um / line_center)
                              .to(u.kilometer / u.second)).value
        # dv/dlambda at the line center, which turns an area integrated over
        # wavelength into one integrated over velocity.
        vel_per_um = (const.c / (line_center * u.micron)).to(
            u.kilometer / u.second / u.micron).value

        maps = {role: [np.full(self.im_shape, np.nan) for _ in range(3)]
                for role in self._MOMENT_ROWS}

        for idx, _ in enumerate(self.fitparams["XPIX"]):
            components = []
            for comp in range(1, 4):
                amp = self.fitparams[f"G{comp}AMP"][idx]
                cen = self.fitparams[f"G{comp}CEN"][idx]
                sig = self.fitparams[f"G{comp}SIGMA"][idx]
                if not np.all(np.isfinite([amp, cen, sig])):
                    continue
                components.append((abs(sig), amp, cen))
            if not components:
                continue
            components.sort(key=lambda component: component[0])
            chosen = {"narrow": components[0]}
            if len(components) > 1:
                chosen["broad"] = components[-1]

            x_pix = self.fitparams["XPIX"][idx]
            y_pix = self.fitparams["YPIX"][idx]
            for role, (sig, amp, cen) in chosen.items():
                sigma_vel = ((const.c * sig / line_center)
                             .to(u.kilometer / u.second)).value
                if subtract_lsf:
                    sigma_vel = np.sqrt(
                        max(sigma_vel ** 2 - sigma_inst_vel ** 2, 0.0))
                maps[role][0][x_pix][y_pix] = amp * vel_per_um
                maps[role][1][x_pix][y_pix] = self.wv_to_vel(
                    cen, line_center).value
                maps[role][2][x_pix][y_pix] = sigma_vel ** 2
        return maps

    def _draw_moment_panel(self, ax, data, moment):
        """Draws one moment map and returns its image.

        Each panel scales to its own data rather than sharing limits down a
        column. Pooling the narrow and broad component onto one stretch is the
        failure mode `render_sigma_plot` documents: the broad component's range
        pushes the narrow one into the bottom of the color bar, and the
        structure that the narrow map exists to show goes flat.
        """
        if moment == 0:
            # Percentile-clipped log stretch, as in `render_totflux_plot`: a
            # few near-zero spaxels would otherwise own decades of the scale.
            positive = data[np.isfinite(data) & (data > 0)]
            norm = None
            if len(positive) > 0 and np.percentile(positive, 2) < np.percentile(positive, 99.5):
                norm = LogNorm(vmin=np.percentile(positive, 2),
                               vmax=np.percentile(positive, 99.5))
            return ax.imshow(data, cmap=self.pltparams.colormap("plasma"),
                             norm=norm, origin="lower")

        if moment == 1:
            # Symmetric limits so the diverging map stays centered on
            # systemic; bound-pinned spaxels are excluded and the result
            # capped at three robust sigma (see `render_rel_vel_plot`).
            vels = self._unpinned(data, self.line_dict[self.name][6])
            limit = 1.0
            if len(vels) > 0:
                limit = min(np.percentile(np.abs(vels), 98),
                            3 * self._robust_spread(vels))
            if not limit > 0:
                limit = 1.0
            return ax.imshow(data, vmin=-limit, vmax=limit,
                             cmap=self.pltparams.colormap("bwr"),
                             origin="lower")

        # M2 is a variance, so the bound it can pin against is the sigma
        # bound squared.
        variances = self._unpinned(data, self.line_dict[self.name][10] ** 2)
        low, high = None, None
        if len(variances) > 0:
            low = np.percentile(variances, 1)
            high = min(np.percentile(variances, 99),
                       np.median(variances) + 3 * self._robust_spread(variances))
            if high <= low:
                high = low + 1.0
        return ax.imshow(data, vmin=low, vmax=high,
                         cmap=self.pltparams.colormap("gist_heat"),
                         origin="lower")

    def render_moment_panels(self, file_stem=None, use_wcs=True, savefig=False,
                             subtract_lsf=False):
        """Renders the six-panel moment diagram of the fit.

        Two rows -- the broad component on top, the narrow one below -- by
        three columns, the 0th, 1st and 2nd moment from left to right, as
        defined in `_component_moment_maps`. Spaxels where the component was
        not fitted are NaN and paint as empty sky, so an all-blank top row
        means no spaxel needed two components.

        Only the bottom-left panel carries the scalebar, compass and axis
        labels; repeating them on all six crowds out the maps.
        """
        maps = self._component_moment_maps(subtract_lsf=subtract_lsf)

        self.pltparams.apply()
        fig = plt.figure(figsize=(24, 15))
        for row, role in enumerate(self._MOMENT_ROWS):
            for moment in range(3):
                position = row * 3 + moment + 1
                ax = fig.add_subplot(2, 3, position,
                                     projection=self.wcs if use_wcs else None)
                corner = (row == len(self._MOMENT_ROWS) - 1 and moment == 0)
                self._decorate_map_axes(ax, use_wcs, scalebar=corner,
                                        axis_labels=corner)
                image = self._draw_moment_panel(ax, maps[role][moment], moment)
                # fraction/pad keep the bar against its own panel; the
                # matplotlib default reserves enough width that it drifts into
                # the gap between columns and reads as belonging to neither.
                cax = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
                self.pltparams.label_colorbar(cax, self._MOMENT_UNITS[moment])
                ax.set_title(f"{role.capitalize()} "
                             f"{self._MOMENT_TITLES[moment]}", loc="left")

        # suptitle defaults to `figure.titlesize`, which the scaling does not
        # set, so it comes out at matplotlib's default while everything around
        # it is at poster size.
        fig.suptitle(self._line_title(),
                     fontsize=plt.rcParams["axes.titlesize"])
        fig.tight_layout()
        self._finish_map("moments", file_stem, savefig)

    def render_chisqr_plot(self, file_stem=None, use_wcs=True, savefig=False):
        """Renders a spaxel map of the reduced chi-square of the preferred
        fit, the goodness-of-fit companion to the parameter maps.

        A well-calibrated fit scatters around 1; values >> 1 flag spaxels
        the model fits poorly (unmodelled line structure, bad continuum,
        underestimated noise) and << 1 flag over-fit / over-inflated
        errors. Only fitted spaxels carry a value; the rest stay NaN.
        """
        base_array = np.full(self.im_shape, np.nan)
        redchi = np.asarray(self.fitparams["REDCHI"], dtype=float)
        for idx, _ in enumerate(self.fitparams["XPIX"]):
            base_array[self.fitparams["XPIX"][idx]][self.fitparams["YPIX"][idx]] = redchi[idx]

        # Percentile-clipped limits so a handful of badly-fit spaxels can't
        # flatten the stretch across the well-fit majority.
        finite = base_array[np.isfinite(base_array)]
        if len(finite) > 0:
            min_chi = np.percentile(finite, 1)
            max_chi = np.percentile(finite, 99)
        else:
            min_chi, max_chi = None, None

        fig, ax = self._setup_map_axes(use_wcs)
        image = ax.imshow(base_array, vmin=min_chi, vmax=max_chi,
                          cmap=self.pltparams.colormap("cividis"), origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, r"$\chi^2_\nu$")
        ax.set_title(r"Reduced $\chi^2$", loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map("chisqr", file_stem, savefig)

    # q3dfit output rendering
    #
    # These mirror the built-in render_* maps above, but read the parameter
    # maps produced by `q3dfit_fit`/`q3dcollect` (via q3dfit's `LineData`)
    # instead of `self.fitparams`. Components are ordered as q3dfit sorted
    # them at collect time (by sigma, descending: c1 is the broadest), so
    # `comp=1` is the broad component. Requires a completed q3dfit fit --
    # call `q3dfit_fit` first, or `load_q3dfit` to reload a saved one.

    def load_q3dfit(self, filepath=None):
        """Reloads a saved q3dfit config so the render_q3d_* methods can run
        in a fresh session (the analogue of `load_fit` for the q3dfit
        backend). `q3dfit_fit` saves this `.npy` automatically.

        Args:
            filepath (str, optional): Path to the saved ``q3din`` ``.npy``.
                Defaults to ``<output><sanitized name>_q3dfit.npy``.
        """
        if filepath is None:
            label = re.sub(r"[^A-Za-z0-9_.+-]", "", self.name) + "_q3dfit"
            filepath = os.path.join(self.output, label + ".npy")
        self.q3di = np.load(filepath, allow_pickle=True).item()
        self._q3d_fluxnorm = self.q3di.argsreadcube.get("fluxnorm", 1.0)
        self._q3d_linedata_cache = None
        return self.q3di

    def _q3d_linedata(self):
        """Returns q3dfit's `LineData` for this cube's fit (cached), raising
        a clear error if no fit is loaded."""
        if getattr(self, "q3di", None) is None:
            raise RuntimeError(
                "No q3dfit results available; run q3dfit_fit(...) or "
                "load_q3dfit(...) before rendering q3dfit maps.")
        if getattr(self, "_q3d_linedata_cache", None) is None:
            from q3dfit.q3dpro import LineData
            self._q3d_linedata_cache = LineData(self.q3di)
        return self._q3d_linedata_cache

    def _q3d_line_label(self):
        """The q3dfit line label of the loaded fit (e.g. "[SIV]10.51")."""
        return self.q3di.lines[0]

    def _q3d_to_map(self, array):
        """Reorients a q3dfit ``(ncols, nrows)`` array to this cube's map
        layout. q3dfit indexes [col=NAXIS1, row=NAXIS2]; the render maps are
        indexed [XPIX=NAXIS2][YPIX=NAXIS1], so the array is transposed."""
        return np.asarray(array, dtype=float).T

    def render_q3d_multicomponent_plot(self, file_stem=None, use_wcs=True,
                                       savefig=False):
        """q3dfit analogue of `render_multicomponent_plot`: a map of the
        number of components q3dfit kept per spaxel."""
        linedata = self._q3d_linedata()
        base_array = self._q3d_to_map(linedata.get_ncomp(self._q3d_line_label()))

        fig, ax = self._setup_map_axes(use_wcs)
        self._draw_component_map(ax, base_array)
        self._finish_map("q3d_components", file_stem, savefig)

    def render_q3d_totflux_plot(self, file_stem=None, use_wcs=True, savefig=False):
        """q3dfit analogue of `render_totflux_plot`: a map of the total line
        flux (summed over components), in erg/s/cm^2.

        q3dfit reports the flux in units of `fluxnorm`; it is multiplied back
        here to recover physical erg/s/cm^2."""
        linedata = self._q3d_linedata()
        flux = linedata.get_flux(self._q3d_line_label(), FLUXSEL="ftot")["flux"]
        base_array = self._q3d_to_map(flux) * self._q3d_fluxnorm

        fig, ax = self._setup_map_axes(use_wcs)
        # Percentile-clipped log stretch, as in render_totflux_plot.
        positive = base_array[np.isfinite(base_array) & (base_array > 0)]
        if len(positive) > 0 and np.percentile(positive, 2) < np.percentile(positive, 99.5):
            norm = LogNorm(vmin=np.percentile(positive, 2),
                           vmax=np.percentile(positive, 99.5))
        else:
            norm = None
        image = ax.imshow(base_array, cmap=self.pltparams.colormap("plasma"),
                          norm=norm, origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, r"[erg s$^{-1}$ cm$^{-2}$]")
        ax.set_title("Total Flux", loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map("q3d_total_flux", file_stem, savefig)

    def render_q3d_rel_vel_plot(self, comp=1, file_stem=None, use_wcs=True,
                                savefig=False):
        """q3dfit analogue of `render_rel_vel_plot`: a map of a component's
        velocity offset from the systemic line center.

        The offset is computed from q3dfit's fitted (observed-frame) central
        wavelength relative to the systemic observed wavelength
        (rest line center x (1 + redshift))."""
        linedata = self._q3d_linedata()
        wav = linedata.get_wave(self._q3d_line_label(), COMPSEL=comp)["wav"]
        systemic = self.line_center * (1 + self.redshift)
        rel_vel = (const.c * (self._q3d_to_map(wav) - systemic) / systemic
                   ).to(u.kilometer / u.second).value

        # Symmetric, robustly capped color limits, as in render_rel_vel_plot.
        # No bound is passed: q3dfit's limits are its own, not this module's
        # `center_cutoff`, so the cap does the work here on its own.
        vels = self._unpinned(rel_vel)
        vel_limit = 1.0
        if len(vels) > 0:
            vel_limit = min(np.percentile(np.abs(vels), 98),
                            3 * self._robust_spread(vels))
        if not vel_limit > 0:
            vel_limit = 1.0

        fig, ax = self._setup_map_axes(use_wcs)
        image = ax.imshow(rel_vel, vmin=-vel_limit, vmax=vel_limit,
                          cmap=self.pltparams.colormap("bwr"), origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, "[km/s]")
        ax.set_title(r"$v_{" + str(comp) + r"}$", loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map(f"q3d_vel_c{comp}", file_stem, savefig)

    def render_q3d_sigma_plot(self, comp=1, file_stem=None, use_wcs=True,
                              savefig=False):
        """q3dfit analogue of `render_sigma_plot`: a map of a component's
        velocity dispersion, in km/s.

        No instrumental-width subtraction is applied here (unlike
        `render_sigma_plot`): q3dfit convolves its models with the
        instrumental LSF during the fit (see `spect_convol`), so the fitted
        sigma is already the intrinsic dispersion."""
        linedata = self._q3d_linedata()
        sigma = linedata.get_sigma(self._q3d_line_label(), COMPSEL=comp)["sig"]
        base_array = self._q3d_to_map(sigma)

        # Robustly capped upper limit, as in render_sigma_plot.
        disps = self._unpinned(base_array)
        if len(disps) > 0:
            min_disp = np.percentile(disps, 1)
            max_disp = min(np.percentile(disps, 99),
                           np.median(disps) + 3 * self._robust_spread(disps))
            if max_disp <= min_disp:
                max_disp = min_disp + 1.0
        else:
            min_disp, max_disp = None, None

        fig, ax = self._setup_map_axes(use_wcs)
        image = ax.imshow(base_array, vmin=min_disp, vmax=max_disp,
                          cmap=self.pltparams.colormap("gist_heat"), origin="lower")
        cax = plt.colorbar(image)
        self.pltparams.label_colorbar(cax, "[km/s]")
        ax.set_title(r"$\sigma_{" + str(comp) + r"}$", loc="right")
        ax.set_title(self._line_title(), loc="left")
        self._finish_map(f"q3d_sigma_c{comp}", file_stem, savefig)

    def render_q3d_spaxel_fit(self, x_pix, y_pix, file_stem=None, savefig=False):
        """q3dfit analogue of `render_spaxel_fit`: plots the fitted line
        model over the data for a single spaxel, via q3dfit's own `q3dout`
        line plot.

        The (x_pix, y_pix) indexing matches `render_spaxel_fit`; it is mapped
        to q3dfit's unity-offset (col, row)."""
        if getattr(self, "q3di", None) is None:
            raise RuntimeError(
                "No q3dfit results available; run q3dfit_fit(...) or "
                "load_q3dfit(...) before rendering q3dfit spaxels.")
        from q3dfit.q3dout import load_q3dout
        q3do = load_q3dout(self.q3di, y_pix + 1, x_pix + 1)
        plotargs = {"nx": 1, "ny": 1, "line": [self._q3d_line_label()],
                    "figsize": [8, 5]}
        outfile = None
        if savefig:
            self.pltparams.apply()
            outfile = self._plot_path(
                file_stem, f"{self.name}_q3d_spaxel_{x_pix}_{y_pix}"
                           f"{self.pltparams.file_suffix}")
        q3do.plot_line(self.q3di, savefig=savefig, outfile=outfile,
                       plotargs=plotargs)

    def render_q3d_maps(self, file_stem=None, use_wcs=True, savefig=False,
                        maxncomp=2):
        """Renders the standard set of q3dfit spaxel maps (component count,
        total flux, and per-component velocity offset and dispersion),
        mirroring the built-in `render_*` set driven by
        `generate_spaxel_map.render_maps`."""
        self.render_q3d_multicomponent_plot(file_stem, use_wcs, savefig)
        self.render_q3d_totflux_plot(file_stem, use_wcs, savefig)
        for comp in range(1, maxncomp + 1):
            self.render_q3d_rel_vel_plot(comp, file_stem, use_wcs, savefig)
            self.render_q3d_sigma_plot(comp, file_stem, use_wcs, savefig)

    def render_master_plot(self):
        """This renders a combination of the component, flux, relative
        velocity, and velocity dispersion plots."""
    
    def gaussian_integral(amplitude, center, sigma, x_1, x_2):
        return amplitude * (sp.erf((center - x_1)) / (np.sqrt(2) * sigma) - sp.erf((center - x_2)) / (np.sqrt(2) * sigma)) / 2
    
    def snr_cut(self, snr_threshold):
        """ 
        Removes and flattens fitted component array based on an SNR threshold.
        """
        
        for i in range(1, 4):
            pass
        
        pass
    
    
    def calc_line_flux(self):
        pass
    
    #def remove_component(self):
    #    
    
    