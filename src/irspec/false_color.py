"""
false_color.py
==============

Generate 3-color (RGB) false-color images from JWST/HST FITS imaging.

The imaging under ``irspec/src/irspec/image_data/{nircam,miri,hubble}`` all covers
the same sky field (the LIRG IR23128) but at different pixel scales, array shapes,
and WCS orientations, so the frames must be reprojected onto a common grid before
they can be composited. This module:

1. loads the science plane and celestial WCS from each FITS file,
2. reprojects all three channels onto an optimal common WCS
   (``reproject.mosaicking.find_optimal_celestial_wcs``), which covers the union
   of the input footprints,
3. stretches each channel and composites them into an RGB image using one of two
   selectable methods, and
4. writes the result to ``images/false_color/``.

The input file list is interpreted in ``[R, G, B]`` order: the first file becomes
the red channel, the second green, the third blue. Conventionally the longest
wavelength band is placed in red.

Three stretch methods are available via the ``method`` toggle:

* ``"lupton"``     -- ``astropy.visualization.make_lupton_rgb`` (Lupton asinh).
* ``"percentile"`` -- per-channel ``PercentileInterval`` + ``AsinhStretch``.
* ``"dark"``       -- per-channel black point pinned just above the robust sky
  level, white point at a high percentile, then ``AsinhStretch``. Use this when
  the background has to render black while the galaxy keeps internal structure;
  ``"percentile"`` puts its black point far down in the noise, which greys the
  sky out.

Because the frames carry very different units (JWST in MJy/sr, HST in
ELECTRONS/S), all methods normalize each channel independently before combining.

Two output grids are available:

* the default *union* grid, built by ``find_optimal_celestial_wcs`` at the finest
  input pixel scale over the union of the input footprints (peak memory scales
  with that -- a full-resolution NIRCam frame produces a large output grid), and
* a *cutout* grid (``center`` + ``size``), a north-up TAN grid of fixed angular
  size centred on a sky position. Use this to frame the merger itself instead of
  the full mosaic footprint; it also keeps the output grid small when one of the
  channels is a full-resolution NIRCam mosaic.

Run the built-in example from the repo root with the pixi environment::

    pixi run python irspec/src/irspec/false_color.py

or drive it from the command line::

    pixi run python irspec/src/irspec/false_color.py \\
        --r path/to/red.fits --g path/to/green.fits --b path/to/blue.fits \\
        --out my_rgb --method percentile

    pixi run python irspec/src/irspec/false_color.py \\
        --center 348.944566 -59.053651 --size 24 --method dark --out merger

The bundled example carries two named framings of IR23128-5919, selected with
``--view``: ``merger`` (20 arcsec, on the two nuclei) and ``tails`` (64 arcsec,
opening up the faint end in NIRCam to show the tidal tails). The three bands are
chosen with ``--bands`` from the ``BANDS`` registry, in ``[R, G, B]`` order::

    pixi run python irspec/src/irspec/false_color.py --view tails
    pixi run python irspec/src/irspec/false_color.py --bands f770 f356 f150

Stretches are tuned per (band trio, view) in ``STRETCHES``; an untuned
combination falls back to the default trio's stretch and says so. Override any
of it with ``--black-sigma`` / ``--white-sigma`` / ``--asinh-a`` / ``--smooth``.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.visualization import (
    make_lupton_rgb,
    PercentileInterval,
    AsinhStretch,
)
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

from irspec import paths

# The JWST/HST imaging and the output location for rendered false-color
# plots both live outside the package; they are resolved lazily through
# irspec.paths (IRSPEC_DATA_ROOT / IRSPEC_OUTPUT_ROOT) so importing this
# module needs no configuration.


def _default_data_dir():
    return paths.image_data_dir()


def _default_out_dir():
    return paths.output_root() / "images" / "false_color"

# Midpoint of the two IR23128-5919 nuclei, measured from the F770W peaks (they
# sit ~4.9 arcsec apart, N and S).
MERGER_CENTER = SkyCoord(348.944566, -59.053651, unit="deg")

# Named framings of the merger, each with a stretch tuned by eye against the
# F770W/F560W/F200W cutout. Both put the black point well above the sky so the
# background renders black, and both set per-channel white points: F560W carries
# a far larger peak-to-sky ratio than the other two bands (the southern nucleus
# dominates it), so a common cut floods the frame with green MIRI diffraction
# spikes.
#
# "tails" additionally drops the black point on the blue channel alone. The tidal
# tails are stellar, not dusty -- F770W has essentially no signal past ~25 arcsec
# -- so the faint end only needs opening up in NIRCam. Lowering all three
# together instead pulls up MIRI 1/f striping across the whole frame and buys no
# extra tail.
#
# 64 arcsec is the useful limit for "tails": the northern arc closes inside it,
# and beyond that the frame is mostly empty sky. Note that a north-up square of
# that size runs past the corners of the (rotated) MIRI mosaic, so ~26% of the
# frame has no MIRI coverage. That costs nothing here because those corners are
# blank sky in every band, but it does mean colour out there is NIRCam-only.
#
# "tails" also coarsens the output grid. The default pixel scale is the finest
# input (NIRCam's 0.031"/pix), which is right for the merger view but oversamples
# a low-surface-brightness view: at native sampling the NIRCam 1/f striping sits
# right at the level of the tails and is plainly visible across them. Binning to
# roughly the MIRI scale averages the stripes down while costing no real detail
# at these surface brightnesses.
VIEWS = {
    "merger": dict(size=20.0, pixscale=None),
    "tails": dict(size=64.0, pixscale=0.10),
}

# The bundled IR23128 imaging, by short band name. The NIRCam frames are the
# CAL_VER 3.0.0 "image3" reduction; the older IR23128_*.fits products they
# replaced are gone.
BANDS = {
    "f560": ("miri", "jw03368-o047_t031_miri_f560w-brightsky_i2d.fits"),
    "f770": ("miri", "jw03368-o047_t031_miri_f770w-brightsky_i2d.fits"),
    "f1500": ("miri", "jw03368-o047_t031_miri_f1500w-brightsky_i2d.fits"),
    "f150": ("nircam", "image3_f150w_i2d.fits"),
    "f200": ("nircam", "image3_f200w_i2d.fits"),
    "f277": ("nircam", "image3_f277w_i2d.fits"),
    "f356": ("nircam", "image3_f356w_i2d.fits"),
}

DEFAULT_BANDS = ("f770", "f560", "f200")

# Stretches, keyed by (band trio, view), each optionally carrying its own
# per-channel ``smooth``. They have to be per trio, not per view alone: a white
# point in sky sigma is comparable across reductions of the *same* band (which is
# why these survived the NIRCam re-reduction untouched), but not across different
# bands, whose peak-to-sky ratios differ.
#
# The "tails" entries drop the black point on whichever channel actually carries
# the tidal tails and hold the others high, and smooth that channel to trade
# resolution for surface brightness. Which channel that is depends on the trio,
# and it is worth measuring rather than assuming -- counting pixels above 3 sigma
# in radial annuli, the tails show up in:
#
#   band    10-15"  15-20"  20-25"  25-30"
#   F356W     39%     19%     12%      3%
#   F770W     14%      4%      3%      0%
#   F150W      8%      6%      2%      1%
#
# so f770/f560/f200 opens up blue (F200W) while f770/f356/f150 opens up *green*
# (F356W), which is the strongest tail tracer of the lot. Lowering all three
# together instead just pulls up MIRI 1/f striping across the frame.
#
# f770/f560/f200 merger -- F560W's white point is pushed high because the
# southern nucleus dominates that band; a common cut floods the frame with green
# MIRI diffraction spikes.
#
# f770/f356/f150 merger -- both short bands are NIRCam stellar tracers, so they
# co-vary and the body washes out to pale blue under the settings above. This
# trio wants a harder black point and a more linear stretch (larger ``a``).
STRETCHES = {
    (("f770", "f560", "f200"), "merger"): dict(
        black_sigma=3.0, white_sigma=(700.0, 1200.0, 300.0), a=0.02
    ),
    (("f770", "f560", "f200"), "tails"): dict(
        black_sigma=(3.0, 3.0, 0.5),
        white_sigma=(800.0, 1400.0, 250.0),
        a=0.01,
        smooth=(0.0, 0.0, 3.0),
    ),
    (("f770", "f356", "f150"), "merger"): dict(
        black_sigma=5.0, white_sigma=(500.0, 900.0, 200.0), a=0.07
    ),
    (("f770", "f356", "f150"), "tails"): dict(
        black_sigma=(3.0, 0.4, 0.9),
        white_sigma=(500.0, 900.0, 200.0),
        a=0.025,
        # Both short bands need it here, not just blue: F356W carries the tails
        # and F150W is the SW frame that carries the striping.
        smooth=(0.0, 2.0, 3.0),
    ),
}


def resolve_bands(bands, data_dir):
    """Map three short band names to FITS paths, in [R, G, B] order."""
    bands = tuple(b.lower() for b in bands)
    if len(bands) != 3:
        raise ValueError(f"Expected 3 bands in [R, G, B] order, got {len(bands)}.")
    unknown = [b for b in bands if b not in BANDS]
    if unknown:
        raise ValueError(
            f"Unknown band(s) {unknown}; choose from {sorted(BANDS)}."
        )
    data_dir = Path(data_dir)
    return bands, [data_dir / BANDS[b][0] / BANDS[b][1] for b in bands]


def resolve_stretch(bands, view):
    """Look up the tuned stretch for a band trio and view.

    Falls back to the default trio's stretch for that view, which is a starting
    point rather than a tuned result -- the caller is told so it can say as much.

    The returned dict holds ``scale_dark`` keyword arguments plus, optionally,
    a ``smooth`` entry for ``make_false_color``; the caller pops it off.

    Returns
    -------
    stretch : dict
    tuned : bool
        False when the returned stretch is the fallback.
    """
    key = (tuple(bands), view)
    if key in STRETCHES:
        return dict(STRETCHES[key]), True
    return dict(STRETCHES[(DEFAULT_BANDS, view)]), False


def repair_holes(data, *, kernel_width=3.0, max_passes=5):
    """Interpolate over bad pixels that are fully enclosed by valid data.

    "Bad" means either exactly ``0.0`` or non-finite, because the two product
    generations flag pixels differently: the MIRI i2d frames write ``0.0`` for
    pixels rejected during resampling, while the CAL_VER 3.0.0 NIRCam frames
    write NaN. Both conventions also appear as the off-array fill, so the two
    cases are separated by connectivity rather than by value -- a bad region
    touching the frame border is off-array and is left alone, and only enclosed
    regions are interpolated across.

    This matters at the saturated cores of the merger: F770W has a zero patch at
    the southern nucleus and F356W a NaN patch at the same place. Left alone,
    each reprojects to a hole in its own channel and the brightest point of the
    merger takes on the complementary colour (cyan without red, magenta without
    green).

    Border fill is normalised to ``0.0`` on the way out. That matches what
    ``reproject_to_common`` / ``reproject_to_cutout`` produce anyway, keeps it
    out of the interpolation, and is the value ``sky_level`` knows to ignore.

    Arguments
    ---------
    data : numpy.ndarray
        Science image, modified out of place.
    kernel_width : float
        Gaussian kernel width, in pixels, used to fill the holes.
    max_passes : int
        Interpolation passes. A hole wider than the kernel is not filled in one
        pass, so passes repeat until no enclosed bad pixels remain.

    Returns
    -------
    numpy.ndarray
        A copy of ``data`` with enclosed holes filled and border fill set to 0.
    """
    bad = ~np.isfinite(data) | (data == 0.0)
    if not bad.any():
        return data

    labels, _ = ndimage.label(bad)
    border = np.concatenate(
        [labels[0, :], labels[-1, :], labels[:, 0], labels[:, -1]]
    )
    edge_labels = set(np.unique(border)) - {0}
    off_array = bad & np.isin(labels, list(edge_labels))
    interior = bad & ~off_array
    if not interior.any():
        filled = data.copy()
        filled[off_array] = 0.0
        return filled

    filled = data.copy()
    filled[interior] = np.nan
    # Take the off-array fill out of play: it is NaN in the newer products, and
    # interpolate_replace_nans would otherwise bleed the frame edge into it.
    filled[off_array] = 0.0

    kernel = Gaussian2DKernel(kernel_width)
    for _ in range(max_passes):
        filled = interpolate_replace_nans(filled, kernel)
        if np.isfinite(filled).all():
            break
    filled[off_array] = 0.0
    return filled


def load_science(filepath, *, repair=False):
    """Load the science plane and celestial WCS from a FITS file.

    Reads the ``SCI`` extension when present (standard for JWST i2d and HST drz
    products), otherwise falls back to HDU 1.

    Arguments
    ---------
    filepath : str or Path
        Path to the FITS file.
    repair : bool
        If True, interpolate over enclosed bad-pixel holes (see
        ``repair_holes``).

    Returns
    -------
    data : numpy.ndarray
        The science image as ``float64``.
    wcs : astropy.wcs.WCS
        The celestial (2-D) WCS for the science plane.
    """
    with fits.open(filepath) as hdul:
        hdu = hdul["SCI"] if "SCI" in hdul else hdul[1]
        data = np.asarray(hdu.data, dtype=np.float64)
        # WCS construction on JWST/HST headers emits noisy datfix/obsfix
        # FITSFixedWarnings; they do not affect the celestial solution we use.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(hdu.header).celestial
    if repair:
        data = repair_holes(data)
    return data, wcs


def reproject_to_common(entries):
    """Reproject a set of (data, wcs) channels onto an optimal common grid.

    The output WCS covers the union of the input footprints at the finest input
    pixel scale, via ``find_optimal_celestial_wcs``. Pixels that fall outside a
    given channel's footprint (returned as NaN by the reprojection) are set to 0.

    Arguments
    ---------
    entries : list of (numpy.ndarray, astropy.wcs.WCS)
        One ``(data, wcs)`` tuple per channel.

    Returns
    -------
    arrays : list of numpy.ndarray
        The reprojected channels, all sharing ``shape_out``.
    wcs_out : astropy.wcs.WCS
        The common output WCS.
    shape_out : tuple of int
        The common output array shape ``(ny, nx)``.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs_out, shape_out = find_optimal_celestial_wcs(
            [(data, wcs) for data, wcs in entries]
        )

    arrays = _reproject_all(entries, wcs_out, shape_out)
    return arrays, wcs_out, shape_out


def _reproject_all(entries, wcs_out, shape_out):
    """Reproject every ``(data, wcs)`` entry onto a given output grid."""
    arrays = []
    for data, wcs in entries:
        arr, _ = reproject_interp((data, wcs), wcs_out, shape_out=shape_out)
        arrays.append(np.nan_to_num(arr, nan=0.0))
    return arrays


def pixel_scale_arcsec(wcs):
    """Return the mean absolute pixel scale of a celestial WCS, in arcsec."""
    return float(np.mean(np.abs(proj_plane_pixel_scales(wcs))) * 3600.0)


def make_cutout_wcs(center, size, pixscale):
    """Build a north-up, east-left TAN grid centred on a sky position.

    Unlike ``find_optimal_celestial_wcs``, which spans the union of the input
    footprints, this fixes the field of view up front. That both frames the
    target and bounds the output grid, which matters when one channel is a
    full-resolution NIRCam mosaic.

    Arguments
    ---------
    center : astropy.coordinates.SkyCoord
        Centre of the output field.
    size : float
        Field of view on a side, in arcsec.
    pixscale : float
        Output pixel scale, in arcsec/pixel.

    Returns
    -------
    wcs_out : astropy.wcs.WCS
        The cutout WCS.
    shape_out : tuple of int
        The output array shape ``(ny, nx)``.
    """
    npix = max(int(round(size / pixscale)), 1)
    wcs_out = WCS(naxis=2)
    wcs_out.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    # FITS pixel centres are 1-based, so the grid centre sits at (npix + 1) / 2.
    wcs_out.wcs.crpix = [(npix + 1) / 2.0, (npix + 1) / 2.0]
    wcs_out.wcs.crval = [center.ra.deg, center.dec.deg]
    wcs_out.wcs.cdelt = [-pixscale / 3600.0, pixscale / 3600.0]
    return wcs_out, (npix, npix)


def reproject_to_cutout(entries, center, size, pixscale=None):
    """Reproject channels onto a fixed-size cutout grid centred on a position.

    Arguments
    ---------
    entries : list of (numpy.ndarray, astropy.wcs.WCS)
        One ``(data, wcs)`` tuple per channel.
    center : astropy.coordinates.SkyCoord
        Centre of the output field.
    size : float
        Field of view on a side, in arcsec.
    pixscale : float or None
        Output pixel scale in arcsec/pixel. Defaults to the finest input scale,
        so the sharpest channel is not degraded.

    Returns
    -------
    arrays : list of numpy.ndarray
    wcs_out : astropy.wcs.WCS
    shape_out : tuple of int
    """
    if pixscale is None:
        pixscale = min(pixel_scale_arcsec(wcs) for _, wcs in entries)
    wcs_out, shape_out = make_cutout_wcs(center, size, pixscale)
    return _reproject_all(entries, wcs_out, shape_out), wcs_out, shape_out


def scale_lupton(r, g, b, *, Q=8, stretch=5, percentile=99.5):
    """Composite three channels into RGB using the Lupton asinh method.

    Each channel is first clipped to a common percentile interval so the very
    different native units (MJy/sr vs ELECTRONS/S) sit on a comparable footing,
    then combined with ``make_lupton_rgb``.

    Arguments
    ---------
    r, g, b : numpy.ndarray
        Reprojected red, green and blue channels (same shape).
    Q : float
        Asinh softening parameter passed to ``make_lupton_rgb``.
    stretch : float
        Linear stretch parameter passed to ``make_lupton_rgb``.
    percentile : float
        Central percentile kept when normalizing each channel.

    Returns
    -------
    numpy.ndarray
        ``(ny, nx, 3)`` ``uint8`` RGB image.
    """
    interval = PercentileInterval(percentile)
    r, g, b = (interval(chan) for chan in (r, g, b))
    return make_lupton_rgb(r, g, b, Q=Q, stretch=stretch)


def scale_percentile(r, g, b, *, percentile=99.5, a=0.1):
    """Composite three channels into RGB using per-channel asinh normalization.

    Each channel is independently clipped to a percentile interval and passed
    through an ``AsinhStretch`` before stacking, giving per-channel control over
    contrast.

    Arguments
    ---------
    r, g, b : numpy.ndarray
        Reprojected red, green and blue channels (same shape).
    percentile : float
        Central percentile kept when normalizing each channel.
    a : float
        Asinh stretch parameter (smaller -> stronger faint-end boost).

    Returns
    -------
    numpy.ndarray
        ``(ny, nx, 3)`` float RGB image in ``[0, 1]``.
    """
    interval = PercentileInterval(percentile)
    stretch = AsinhStretch(a)
    channels = []
    for chan in (r, g, b):
        scaled = stretch(interval(chan))
        channels.append(np.clip(np.nan_to_num(scaled, nan=0.0), 0.0, 1.0))
    return np.dstack(channels)


def smooth_channels(arrays, smooth):
    """Gaussian-smooth each channel by a per-channel sigma, in output pixels.

    Trading resolution for surface brightness. This is worth doing on a wide,
    faint view: NIRCam's 1/f striping lands at the same level as the tidal
    tails, and smoothing suppresses it while the tails, being extended, survive.
    It also lowers the measured sky sigma, so a black point expressed in sigma
    follows it down and fainter structure clears the cut.

    Note the smoothing runs after reprojection, so it bleeds across the edge of a
    channel's footprint (where the fill value is 0). That is harmless as long as
    the footprint edge falls on blank sky, which is the case for these frames.

    Arguments
    ---------
    arrays : sequence of numpy.ndarray
        The three reprojected channels.
    smooth : float or sequence of 3 floats
        Gaussian sigma per channel, in output pixels. Zero leaves a channel
        untouched.

    Returns
    -------
    list of numpy.ndarray
    """
    return [
        ndimage.gaussian_filter(chan, sig) if sig else chan
        for chan, sig in zip(arrays, _per_channel(smooth))
    ]


def _per_channel(value):
    """Broadcast a scalar to a 3-tuple, or pass a 3-sequence through."""
    if np.isscalar(value):
        return (value,) * 3
    values = tuple(value)
    if len(values) != 3:
        raise ValueError(f"Expected a scalar or 3 values, got {len(values)}.")
    return values


def sky_level(chan, *, ignore_zero=True, sigma=3.0, maxiters=10):
    """Robust (median, sigma) of the sky in one channel.

    Sigma-clipped rather than a plain median, so the estimate does not depend on
    how much of the frame the source fills: on a tight cutout the galaxy covers
    enough pixels to drag a plain median up, which would raise the black point
    and quietly suppress the faint envelope as the field of view is narrowed.

    Exactly-zero pixels are the fill value written by ``reproject_to_common`` /
    ``reproject_to_cutout`` outside a channel's footprint, and are excluded by
    default so an unfilled corner cannot pull the sky estimate to zero.

    Returns
    -------
    median : float
        Clipped sky median.
    sigma : float
        Clipped sky standard deviation.
    """
    finite = chan[np.isfinite(chan)]
    if ignore_zero:
        finite = finite[finite != 0.0]
    if finite.size == 0:
        return 0.0, 0.0
    _, median, std = sigma_clipped_stats(finite, sigma=sigma, maxiters=maxiters)
    return float(median), float(std)


def scale_dark(
    r, g, b, *, black_sigma=1.0, white_sigma=None, hi_percentile=99.8, a=0.1
):
    """Composite three channels with the sky pinned to black.

    For each channel the black point is set to ``sky_median + black_sigma *
    sky_sigma``, a white point is chosen, the channel is normalized between them
    and passed through an ``AsinhStretch``. Anchoring black to the measured sky
    rather than to a low percentile is what keeps the background dark: a
    percentile black point sits well inside the noise distribution, so the sky
    lands at mid-grey once the asinh stretch boosts the faint end.

    The white point comes from ``white_sigma`` when given, otherwise from
    ``hi_percentile``. Prefer ``white_sigma`` for colour balance: a percentile
    white point depends on how much of the frame the source fills and on each
    band's own dynamic range, so the three channels end up cut at unrelated
    physical levels and the composite takes on a spurious colour cast. Placing
    white at the same multiple of each channel's sky noise ties the channels to
    a common reference, so colour then reflects the band-to-band contrast of the
    source.

    Arguments
    ---------
    r, g, b : numpy.ndarray
        Reprojected red, green and blue channels (same shape).
    black_sigma : float or sequence of 3 floats
        Black point, in robust sigma above the sky median. Raise it for a darker
        background, lower it (or go negative) to reveal more faint emission.
    white_sigma : float or sequence of 3 floats or None
        White point, in robust sigma above the sky median. Takes precedence over
        ``hi_percentile`` when not None. Raise it to keep bright cores from
        saturating, lower it to brighten the source.
    hi_percentile : float or sequence of 3 floats
        Percentile of the frame mapped to white, used only when ``white_sigma``
        is None.
    a : float or sequence of 3 floats
        Asinh stretch parameter (smaller -> stronger faint-end boost).

    Returns
    -------
    numpy.ndarray
        ``(ny, nx, 3)`` float RGB image in ``[0, 1]``.
    """
    use_sigma = white_sigma is not None
    whites = _per_channel(white_sigma if use_sigma else hi_percentile)

    channels = []
    for chan, sig, white_spec, aa in zip(
        (r, g, b), _per_channel(black_sigma), whites, _per_channel(a)
    ):
        median, sigma = sky_level(chan)
        black = median + sig * sigma
        white = (
            median + white_spec * sigma
            if use_sigma
            else float(np.nanpercentile(chan, white_spec))
        )
        if not white > black:
            # Degenerate channel (flat, or clipped away by the black point);
            # fall back to a unit span so the division stays finite.
            white = black + 1.0
        scaled = AsinhStretch(aa)(np.clip((chan - black) / (white - black), 0.0, 1.0))
        channels.append(np.clip(np.nan_to_num(scaled, nan=0.0), 0.0, 1.0))
    return np.dstack(channels)


# Available stretch/composite methods, keyed by the ``method`` toggle.
_SCALERS = {
    "lupton": scale_lupton,
    "percentile": scale_percentile,
    "dark": scale_dark,
}


def make_false_color(
    files,
    output_name,
    *,
    method="lupton",
    show_wcs=False,
    center=None,
    size=None,
    pixscale=None,
    outdir=None,
    repair=True,
    smooth=0.0,
    **kwargs,
):
    """Build and save a 3-color false-color image from three FITS files.

    Arguments
    ---------
    files : sequence of (str or Path)
        Exactly three FITS files in ``[R, G, B]`` order.
    output_name : str
        Base name (without extension) for the output PNG.
    method : {"lupton", "percentile", "dark"}
        Stretch/composite method to use.
    show_wcs : bool
        If True, render with RA/Dec coordinate axes; otherwise a clean,
        axis-free image.
    center : astropy.coordinates.SkyCoord or None
        Centre of a fixed-size cutout. If None, the union grid is used.
    size : float or None
        Cutout field of view on a side, in arcsec. Required with ``center``.
    pixscale : float or None
        Output pixel scale for the cutout, in arcsec/pixel. Defaults to the
        finest input scale.
    outdir : str or Path or None
        Directory for the output PNG. Defaults to images/false_color under
        the irspec output root.
    repair : bool
        Interpolate over enclosed bad-pixel holes in each input (see
        ``repair_holes``).
    smooth : float or sequence of 3 floats
        Per-channel Gaussian sigma, in output pixels, applied after reprojection
        (see ``smooth_channels``).
    **kwargs
        Extra keyword arguments forwarded to the selected scaler
        (e.g. ``Q``, ``stretch``, ``percentile``, ``black_sigma``, ``a``).

    Returns
    -------
    pathlib.Path
        Path to the written PNG.
    """
    files = list(files)
    if len(files) != 3:
        raise ValueError(
            f"Expected exactly 3 files in [R, G, B] order, got {len(files)}."
        )
    if method not in _SCALERS:
        raise ValueError(
            f"Unknown method {method!r}; choose from {sorted(_SCALERS)}."
        )
    if (center is None) != (size is None):
        raise ValueError("center and size must be given together.")

    entries = [load_science(f, repair=repair) for f in files]
    if center is None:
        (r, g, b), wcs_out, _ = reproject_to_common(entries)
    else:
        (r, g, b), wcs_out, _ = reproject_to_cutout(
            entries, center, size, pixscale=pixscale
        )

    if np.any(_per_channel(smooth)):
        r, g, b = smooth_channels((r, g, b), smooth)

    rgb = _SCALERS[method](r, g, b, **kwargs)

    if show_wcs:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection=wcs_out)
        ax.imshow(rgb, origin="lower")
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
        ax.coords.grid(color="white", alpha=0.2)
    else:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(rgb, origin="lower")
        ax.axis("off")

    out_dir = Path(outdir) if outdir is not None else _default_out_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{output_name}.png"
    fig.savefig(out_path, dpi=400, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    """Command-line entry point / bundled example.

    With no ``--r/--g/--b`` arguments it renders a default NIRCam+MIRI tri-band
    from the bundled imaging.
    """
    parser = argparse.ArgumentParser(
        description="Generate a 3-color false-color image from three FITS files."
    )
    parser.add_argument("--r", help="FITS file for the red channel (longest wavelength).")
    parser.add_argument("--g", help="FITS file for the green channel.")
    parser.add_argument("--b", help="FITS file for the blue channel (shortest wavelength).")
    parser.add_argument("--out", default="false_color", help="Output PNG base name.")
    parser.add_argument(
        "--method",
        default="lupton",
        choices=sorted(_SCALERS),
        help="Stretch/composite method.",
    )
    parser.add_argument(
        "--show-wcs",
        action="store_true",
        help="Render with RA/Dec coordinate axes.",
    )
    parser.add_argument(
        "--center",
        nargs=2,
        type=float,
        metavar=("RA", "DEC"),
        help="Cutout centre in decimal degrees. Requires --size.",
    )
    parser.add_argument(
        "--size",
        type=float,
        help="Cutout field of view on a side, in arcsec. Requires --center.",
    )
    parser.add_argument(
        "--pixscale",
        type=float,
        help="Cutout pixel scale in arcsec/pixel (default: finest input scale).",
    )
    parser.add_argument(
        "--view",
        default="merger",
        choices=sorted(VIEWS),
        help="Named framing + stretch preset for the bundled example.",
    )
    parser.add_argument(
        "--black-sigma",
        nargs="+",
        type=float,
        metavar="S",
        help="Override the preset black point, in sky sigma. One value, or three "
        "for R G B.",
    )
    parser.add_argument(
        "--white-sigma",
        nargs="+",
        type=float,
        metavar="S",
        help="Override the preset white point, in sky sigma. One value, or three "
        "for R G B.",
    )
    parser.add_argument(
        "--asinh-a",
        type=float,
        help="Override the preset asinh stretch parameter (smaller -> stronger "
        "faint-end boost).",
    )
    parser.add_argument(
        "--smooth",
        nargs="+",
        type=float,
        metavar="SIGMA",
        help="Override the preset Gaussian smoothing, in output pixels. One "
        "value, or three for R G B.",
    )
    parser.add_argument(
        "--bands",
        nargs=3,
        default=list(DEFAULT_BANDS),
        metavar=("R", "G", "B"),
        help=f"Three bundled bands in [R, G, B] order. Choose from "
        f"{sorted(BANDS)}. Default: {' '.join(DEFAULT_BANDS)}.",
    )
    parser.add_argument(
        "--data-dir",
        default=None,
        help="Directory holding the survey imaging (default: image_data under "
        "the irspec data root).",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Directory for the output PNG (default: images/false_color under "
        "the irspec output root).",
    )
    args = parser.parse_args()

    center = SkyCoord(args.center[0], args.center[1], unit="deg") if args.center else None

    def _one_or_three(values, flag):
        """Accept a single value or an explicit R G B triple from the CLI."""
        if values is None:
            return None
        if len(values) == 1:
            return values[0]
        if len(values) == 3:
            return tuple(values)
        raise SystemExit(f"{flag} takes 1 or 3 values, got {len(values)}.")

    overrides = {
        "black_sigma": _one_or_three(args.black_sigma, "--black-sigma"),
        "white_sigma": _one_or_three(args.white_sigma, "--white-sigma"),
        "a": args.asinh_a,
    }
    overrides = {k: v for k, v in overrides.items() if v is not None}

    if args.r and args.g and args.b:
        files = [args.r, args.g, args.b]
        out_path = make_false_color(
            files,
            args.out,
            method=args.method,
            show_wcs=args.show_wcs,
            center=center,
            size=args.size,
            pixscale=args.pixscale,
            outdir=args.outdir,
            smooth=_one_or_three(args.smooth, "--smooth") or 0.0,
            **(overrides if args.method == "dark" else {}),
        )
        print(f"Wrote {out_path}")
        return

    # Bundled example: the IR23128-5919 merger in the selected bands (longest
    # wavelength -> red), framed by the selected view with the sky pinned to
    # black.
    try:
        bands, files = resolve_bands(args.bands, args.data_dir or _default_data_dir())
    except ValueError as exc:
        raise SystemExit(str(exc))
    view = VIEWS[args.view]
    stretch, tuned = resolve_stretch(bands, args.view)
    preset_smooth = stretch.pop("smooth", 0.0)
    if not tuned and not overrides:
        print(
            f"Note: no stretch tuned for {'/'.join(bands)} + {args.view}; "
            f"falling back to {'/'.join(DEFAULT_BANDS)}. Expect to adjust it "
            "with --black-sigma / --white-sigma / --asinh-a."
        )

    method = args.method if args.method != "lupton" else "dark"
    default_name = f"ir23128_{'_'.join(bands)}_{args.view}"
    out_path = make_false_color(
        files,
        args.out if args.out != "false_color" else default_name,
        method=method,
        show_wcs=args.show_wcs,
        center=center if center is not None else MERGER_CENTER,
        size=args.size if args.size is not None else view["size"],
        pixscale=args.pixscale if args.pixscale is not None else view["pixscale"],
        smooth=_one_or_three(args.smooth, "--smooth") if args.smooth else preset_smooth,
        outdir=args.outdir,
        **({**stretch, **overrides} if method == "dark" else {}),
    )
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
