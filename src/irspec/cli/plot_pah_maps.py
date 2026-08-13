"""Render spatial maps of the PAH features fitted by a CAFE grid run.

A grid fit (``run_cafe.py`` on a grid param file) leaves a per-spaxel parameter
cube at ``cafe_output/grid/{NAME}_grid/{NAME}_parcube.fits``. Every PAH in it is
stored as a Drude triplet -- ``d<n>_PAH<label>_<wave0>_{Wave,Gamma,Peak}`` --
which describes a profile, not a flux. This script reconstructs CAFE's model at
each spaxel, integrates every PAH profile, and maps the result.

Observed vs intrinsic power
---------------------------
CAFE fits PAHs behind an extinction screen, so a Drude triplet on its own gives
the *intrinsic* (extinction-corrected) feature strength. That quantity is only
as trustworthy as the extinction, and at low signal-to-noise the fit is free to
pair a huge intrinsic amplitude with a huge opacity: in this cube ``PAH_TAU``
runs to >1000 in the faint outskirts, where the intrinsic 6.2 um flux then
exceeds the nucleus by 3x while the observed flux is 100x fainter.

The default maps are therefore the *observed* (extinguished) power -- CAFE's own
``pah_power_ext``, the Drude profile multiplied through by ``extPAH`` -- which is
what the data actually constrain. ``--power intrinsic`` renders the
extinction-corrected version instead, and ``--power both`` does each; read the
intrinsic maps alongside the ``PAH_TAU`` map that is always written.

Reconstructing ``extPAH`` needs the same continuum profiles and input parameter
file the fit used. The profiles depend only on the wavelength grid, so they are
built once for the cube rather than per spaxel, and the fit's own settings are
read back out of the run manifest so they cannot silently disagree.

Spaxels CAFE never fitted are NaN in the parameter cube -- ``make_parcube``
initialises it that way -- and stay masked, reading as empty sky rather than as
a measured zero.

Figures go to ``{plot-root}/pah_maps/{NAME}/``: one map per feature, an overview
panel of every feature on a shared colour scale, and the ``PAH_TAU`` map. The
palette tags each filename (``_light``/``_dark``), so both can live side by side.

Run from anywhere; the script anchors its own paths::

    python plot_pah_maps.py --name pahmap_ch12
    python plot_pah_maps.py --name pahmap_ch12 --power both
    python plot_pah_maps.py --name pahmap_ch12 --palette dark --scaling paper
"""

import argparse
import os
import re
import warnings

import numpy as np
import astropy.io.fits as fits
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import add_scalebar
from specutils import Spectrum1D

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import LogNorm

import CAFE.cafe as cafe
from CAFE.cafe_helper import CAFE_prof_generator, parcube2parobj
from CAFE.cafe_lib import get_model_fluxes, get_feat_pars, mask_spec
from CAFE.component_model import drude_int_fluxes

from irspec.plotparams import PlotParams
from irspec.cli import runner_utils as ru
from irspec import paths

# Luminosity distance to IR 23128, matching the other map scripts.
DISTANCE = 194.99 * u.Mpc
# CAFE/CRETA/plot roots resolve lazily through irspec.paths.

# CAFE names each Drude parameter d<n>_PAH<label>_<wave0>_<attr>, where <wave0>
# is the rest wavelength in units of 1e-4 um (e.g. 62123 -> 6.2123 um).
_DRUDE = re.compile(r"^d\d+_(?P<label>PAH\d+)_(?P<w0>\d+)_(?P<attr>Wave|Gamma|Peak)$")
_WAVE0_SCALE = 1e4


# ---------------------------------------------------------------- run metadata

def read_manifest_entry(root, name, run_type="grid"):
    """Returns the ``runs.txt`` block for one run as a dict, or ``{}``.

    The fit's mode and redshift have to be reproduced exactly to rebuild its
    extinction, so they are read back from the manifest the runner wrote rather
    than re-specified on this command line and left to drift.
    """
    manifest = os.path.join(ru.group_dir(root, run_type), ru.MANIFEST_NAME)
    if not os.path.exists(manifest):
        return {}
    header = f"[{ru.run_dirname(name, run_type)}]"
    entry, inside = {}, False
    with open(manifest) as handle:
        for line in handle:
            if line.strip() == header:
                inside = True
                continue
            if inside:
                if line.strip() and set(line.strip()) == {"-"}:
                    break
                key, _, value = line.partition("=")
                entry[key.strip()] = value.strip()
    return entry


def repair_ctypes(header):
    """Fixes the malformed celestial CTYPEs CRETA writes into its grid cubes.

    A FITS CTYPE is 8 characters, so the declination axis of a tangent
    projection is ``DEC--TAN``. CRETA's cube writer pads it to ``DEC---TAN``,
    which CAFE copies into the parameter cube; wcslib then truncates to
    ``DEC---TA`` and rejects ``-TA`` as a projection code. Normalising the
    separator dashes here keeps the repair in one place instead of every caller
    having to know about it.
    """
    for axis in (1, 2):
        key = f"CTYPE{axis}"
        value = header.get(key)
        if value is None or "-" not in value:
            continue
        stem, _, projection = value.partition("-")
        projection = projection.strip("-")
        if projection:
            header[key] = f"{stem:-<{8 - len(projection)}}{projection}"
    return header


# ------------------------------------------------------------------ model side

def feature_catalogue(parnames):
    """Every PAH feature in the parameter cube, as ``(key, label, wave0)``.

    ``key`` is how ``get_feat_pars`` names the feature in ``drude[3]``: the
    parameter stem with the leading ``d`` dropped (``d3_PAH62_62123_Wave`` ->
    ``3_PAH62_62123``). Matching on it rather than on position matters because
    a spaxel only carries the features its own data covered -- ``make_parobj``
    leaves the rest out and ``parcube2parobj`` skips their NaNs -- so the Drude
    arrays are shorter, and differently ordered, at some spaxels.
    """
    catalogue, seen = [], set()
    for parname in parnames:
        match = _DRUDE.match(str(parname))
        if match is None:
            continue
        key = str(parname)[1:-len(f"_{match['attr']}")]
        if key in seen:
            continue
        seen.add(key)
        catalogue.append((key, match["label"],
                          int(match["w0"]) / _WAVE0_SCALE))
    return catalogue


def load_run(name, cafe_root, creta_root, mode, redshift, extract):
    """Loads a grid run and returns everything the maps need.

    Returns a dict with the celestial ``wcs``, the ``features`` catalogue, the
    per-feature ``observed``/``intrinsic`` power maps (W/m^2), and the
    ``pah_tau`` map.
    """
    cube_dir = os.path.join(creta_root, "grid", f"{name}_grid")
    par_dir = os.path.join(cafe_root, "grid", f"{name}_grid")
    cube_fn = f"{name}_cube.fits"
    parcube_fn = f"{name}_parcube.fits"
    for directory, filename in ((cube_dir, cube_fn), (par_dir, parcube_fn)):
        if not os.path.exists(os.path.join(directory, filename)):
            raise FileNotFoundError(
                f"Missing {os.path.join(directory, filename)}. Run the CRETA "
                f"extraction and the CAFE grid fit for '{name}' first.")

    cafe_dir = os.path.dirname(cafe.__file__) + "/"
    inppar = cafe_dir + f"inp_parfiles/inpars_jwst_nirspec-miri_{mode}.ini"
    optpar = cafe_dir + "opt_parfiles/default_opt.cafe"

    model = cafe.cubemod(cafe_dir)
    model.read_cube(cube_fn, file_dir=cube_dir + "/", extract=extract,
                    trim=True, z=redshift)
    model.read_parcube_file(parcube_fn, file_dir=par_dir + "/")

    with fits.open(os.path.join(par_dir, parcube_fn)) as hdul:
        parnames = [str(row[1]) for row in hdul["PARNAME"].data]
        wcs = WCS(repair_ctypes(hdul["VALUE"].header)).celestial
    features = feature_catalogue(parnames)
    if not features:
        raise ValueError(f"No PAH features in {parcube_fn}.")

    # The continuum profiles depend only on the wavelength grid, not on the
    # fitted values, so build them once off the highest-SNR spaxel and reuse
    # them for every spaxel (this is what CAFE's own make_profcube does).
    snr_image = np.nansum(model.fluxes[10:20, :, :], axis=0)
    y_ref, x_ref = np.unravel_index(np.nanargmax(snr_image), snr_image.shape)
    wave, flux, flux_unc, _, _ = mask_spec(model, x=x_ref, y=y_ref)
    spec = Spectrum1D(spectral_axis=wave * u.micron, flux=flux * u.Jy,
                      uncertainty=StdDevUncertainty(flux_unc), redshift=model.z)
    cont_profs = CAFE_prof_generator(spec, inppar, optpar, None,
                                     cafe_path=cafe_dir).make_cont_profs()

    ny, nx = model.ny, model.nx
    shape = (len(features), ny, nx)
    powers = {"observed": np.full(shape, np.nan),
              "intrinsic": np.full(shape, np.nan)}
    pah_tau = np.full((ny, nx), np.nan)

    index_of = {key: index for index, (key, _, _) in enumerate(features)}
    # Rest-wavelength range each spaxel actually covers. A grid that spans more
    # than one channel's footprint has spaxels carrying only the redder channel,
    # and CAFE still returns a Drude for features outside their coverage --
    # unconstrained, so the fit is free to pair an enormous amplitude with an
    # enormous opacity. Those are extrapolations, not measurements, and are left
    # NaN below. CRETA writes 0.0 rather than NaN for apertures off the edge of
    # a band's array, so zeros are "no data" here too.
    observed = np.isfinite(model.fluxes) & (model.fluxes != 0)
    wave_lo = np.full((ny, nx), np.nan)
    wave_hi = np.full((ny, nx), np.nan)
    for y in range(ny):
        for x in range(nx):
            covered = model.waves[observed[:, y, x]]
            if covered.size:
                wave_lo[y, x] = covered.min()
                wave_hi[y, x] = covered.max()

    fitted = np.isfinite(model.parcube["VALUE"].data[0])
    for y in range(ny):
        for x in range(nx):
            if not fitted[y, x]:
                continue
            params = parcube2parobj(model.parcube, x=x, y=y)
            comp_fluxes, _, ext_comps, _, _, _ = get_model_fluxes(
                params, wave, cont_profs, comps=True)
            _, drude, _ = get_feat_pars(params, apply_vgrad2waves=True)
            if len(drude[0]) == 0:
                continue
            intrinsic = drude_int_fluxes(comp_fluxes["wave"], drude).value
            observed_power = drude_int_fluxes(comp_fluxes["wave"], drude,
                                              ext=ext_comps["extPAH"]).value
            # Scatter into the catalogue by name: this spaxel may carry only a
            # subset of the features, so position in `drude` means nothing.
            # Anything it lacks stays NaN, i.e. "not measured here".
            for position, key in enumerate(drude[3]):
                plane = index_of.get(str(key))
                if plane is None:
                    continue
                if not (wave_lo[y, x] <= features[plane][2] <= wave_hi[y, x]):
                    continue
                powers["intrinsic"][plane, y, x] = intrinsic[position]
                powers["observed"][plane, y, x] = observed_power[position]
            pah_tau[y, x] = params["PAH_TAU"].value

    return {"wcs": wcs, "features": features, "powers": powers,
            "pah_tau": pah_tau, "fitted": fitted}


# ----------------------------------------------------------------- figure side

def feature_title(wave0):
    """Display name for a feature, e.g. ``PAH 6.2 um``."""
    return rf"PAH {wave0:.1f} $\mu$m"


def log_norm(*maps):
    """A log normalisation spanning every supplied map, or ``None`` if they
    hold no positive data.

    The floor is the 1st percentile of the positive values, which keeps a
    handful of near-zero edge spaxels from flattening the rest of the stretch.
    A feature can be in the parameter cube yet fitted at no spaxel at all --
    likelier the wider the catalogue and the smaller the grid -- and taking a
    percentile of that empty set raises, so hand matplotlib its own autoscale
    instead.
    """
    stacked = np.concatenate(
        [m[np.isfinite(m) & (m > 0)].ravel() for m in maps])
    if stacked.size == 0:
        return None
    return LogNorm(vmin=np.percentile(stacked, 1), vmax=stacked.max())


def tidy_coords(ax, show_x, show_y):
    """Thins and labels a panel's celestial axes.

    Three ticks per axis: a WCS axis left to choose its own on a panel a couple
    of arcseconds across writes enough RA labels that they overlap illegibly.
    """
    ax.coords[0].set_ticks(number=3)
    ax.coords[1].set_ticks(number=4)
    ax.coords[0].set_ticklabel(exclude_overlapping=True)
    ax.coords[1].set_ticklabel(exclude_overlapping=True)
    ax.coords[0].set_ticklabel_visible(show_x)
    ax.coords[1].set_ticklabel_visible(show_y)
    if show_x:
        ax.set_xlabel("Right Ascension")
    if show_y:
        ax.set_ylabel("Declination")


def add_kpc_scalebar(ax, pltparams):
    scalebar_angle = (1 * u.kpc / DISTANCE).to(
        u.deg, equivalencies=u.dimensionless_angles())
    fontprops = fm.FontProperties(size=pltparams.annotation_size,
                                  family=pltparams.font)
    add_scalebar(ax, scalebar_angle, label="1 kpc",
                 color=pltparams.foreground, fontproperties=fontprops)


def render_single(data, wcs, pltparams, out_path, title, bar_label,
                  norm=None, cmap="plasma"):
    """Writes one full-size map. Returns the path.

    What the map shows goes in the colorbar label rather than a second title:
    at poster scaling a left/right title pair overruns the panel and the two
    collide in the middle.
    """
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(projection=wcs)
    image = ax.imshow(np.ma.masked_invalid(data), origin="lower",
                      cmap=pltparams.colormap(cmap),
                      norm=norm if norm is not None else log_norm(data))
    pltparams.label_colorbar(fig.colorbar(image, ax=ax, pad=0.02), bar_label)
    tidy_coords(ax, show_x=True, show_y=True)
    add_kpc_scalebar(ax, pltparams)
    ax.set_title(title, loc="left")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_overview(features, maps, wcs, pltparams, out_path, title):
    """Writes the all-features panel on a shared colour scale. Returns the path."""
    ncols = min(4, len(features))
    nrows = int(np.ceil(len(features) / ncols))
    norm = log_norm(*maps)

    fig = plt.figure(figsize=(5.0 * ncols, 4.8 * nrows))
    for position, ((_, _, wave0), data) in enumerate(zip(features, maps)):
        ax = fig.add_subplot(nrows, ncols, position + 1, projection=wcs)
        image = ax.imshow(np.ma.masked_invalid(data), origin="lower",
                          cmap=pltparams.colormap("plasma"), norm=norm)
        ax.set_title(feature_title(wave0), loc="left")
        # Coordinate labels only on the outside of the grid: a panel with a
        # neighbour below it is not the bottom of its column, whatever row it
        # sits in, and the last row is usually short.
        tidy_coords(ax,
                    show_x=position + ncols >= len(features),
                    show_y=position % ncols == 0)

    bar = fig.colorbar(image, ax=fig.get_axes(), pad=0.02, fraction=0.03)
    pltparams.label_colorbar(bar, r"[W m$^{-2}$]")
    fig.suptitle(title, y=0.995)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


POWER_TITLES = {"observed": "Observed", "intrinsic": "Extinction-corrected"}


def render_power(kind, run, pltparams, out_dir, name):
    """Renders every feature map plus the overview for one power definition."""
    maps = [run["powers"][kind][index] for index in range(len(run["features"]))]
    label = POWER_TITLES[kind]
    drawn = []
    for (_, feature_label, wave0), data in zip(run["features"], maps):
        detected = np.isfinite(data) & (data > 0)
        if not detected.any():
            # Catalogued but never fitted anywhere on this grid. Reported
            # rather than skipped silently: an absent map file should not be
            # something the reader has to notice for themselves.
            print(f"  {feature_title(wave0):>18}    0 spaxels  -- not fitted "
                  f"anywhere; no map written")
            continue
        out = os.path.join(
            out_dir, f"{name}_{feature_label}_{kind}{pltparams.file_suffix}.png")
        render_single(data, run["wcs"], pltparams, out,
                      feature_title(wave0), f"{label} flux [W m$^{{-2}}$]")
        drawn.append(((None, feature_label, wave0), data))
        print(f"  {feature_title(wave0):>18}  {int(detected.sum()):3d} spaxels  "
              f"max {np.nanmax(data):.2e} W/m^2  -> {os.path.basename(out)}")

    overview = os.path.join(
        out_dir, f"{name}_overview_{kind}{pltparams.file_suffix}.png")
    render_overview([f for f, _ in drawn], [d for _, d in drawn], run["wcs"],
                    pltparams, overview, f"{name} -- {label.lower()} PAH flux")
    print(f"[plot_pah_maps] overview -> {os.path.basename(overview)}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Render PAH feature maps from a CAFE grid parameter cube.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--name", required=True,
                        help="Grid run name, as passed to run_cafe.py.")
    parser.add_argument("--power", default="observed",
                        choices=["observed", "intrinsic", "both"],
                        help="Map the extinguished power, the extinction-"
                             "corrected power, or both.")
    parser.add_argument("--cafe-root", default=None,
                        help="Root holding the grid/ group with the fit "
                             "(default: cafe_output/ under the irspec output root).")
    parser.add_argument("--creta-root", default=None,
                        help="Root holding the grid/ group with the extraction "
                             "(default: creta_extractions/ under the irspec output root).")
    parser.add_argument("--plot-root", default=None,
                        help="Root under which pah_maps/{NAME}/ is created "
                             "(default: plots/ under the irspec output root).")
    parser.add_argument("--mode", default=None, choices=["AGN", "SB", "SB_and_AGN"],
                        help="CAFE mode the fit used; read from the run "
                             "manifest when omitted.")
    parser.add_argument("--redshift", type=float, default=None,
                        help="Redshift the fit used; read from the run "
                             "manifest when omitted.")
    parser.add_argument("--extract", default=None,
                        help="CRETA flux column the fit used; read from the "
                             "run manifest when omitted.")
    parser.add_argument("--palette", choices=list(PlotParams.PALATTES),
                        default=PlotParams.DEFAULT_PALATTE,
                        help="Figure palette.")
    parser.add_argument("--scaling", choices=list(PlotParams.SCALINGS),
                        default=PlotParams.DEFAULT_SCALING,
                        help="Text scaling.")
    return parser


def main():
    args = build_parser().parse_args()
    if args.cafe_root is None:
        args.cafe_root = str(paths.cafe_dir())
    if args.creta_root is None:
        args.creta_root = str(paths.creta_dir())
    if args.plot_root is None:
        args.plot_root = str(paths.plots_dir())
    warnings.filterwarnings("ignore")

    entry = read_manifest_entry(args.cafe_root, args.name)
    mode = args.mode or entry.get("mode", "AGN")
    redshift = args.redshift if args.redshift is not None else float(
        entry.get("redshift", ru.DEFAULT_REDSHIFT))
    extract = args.extract or entry.get("extract", "Flux_st")
    print(f"[plot_pah_maps] {args.name}: mode={mode} z={redshift} "
          f"extract={extract}")

    run = load_run(args.name, args.cafe_root, args.creta_root, mode, redshift,
                   extract)

    out_dir = os.path.join(args.plot_root, "pah_maps", args.name)
    os.makedirs(out_dir, exist_ok=True)
    pltparams = PlotParams(palatte=args.palette, scaling=args.scaling)

    fitted = run["fitted"]
    tau = run["pah_tau"]
    print(f"[plot_pah_maps] {len(run['features'])} PAH features, "
          f"{int(fitted.sum())} of {fitted.size} spaxels fitted; "
          f"PAH_TAU median {np.nanmedian(tau):.1f}, max {np.nanmax(tau):.0f}")

    for kind in (["observed", "intrinsic"] if args.power == "both"
                 else [args.power]):
        render_power(kind, run, pltparams, out_dir, args.name)

    tau_out = os.path.join(out_dir, f"{args.name}_PAH_TAU{pltparams.file_suffix}.png")
    render_single(tau, run["wcs"], pltparams, tau_out, r"PAH opacity",
                  r"Fitted $\tau_{\rm PAH}$", cmap="viridis")
    print(f"[plot_pah_maps] opacity -> {os.path.basename(tau_out)}")
    print(f"[plot_pah_maps] figures in {out_dir}")


if __name__ == "__main__":
    main()
