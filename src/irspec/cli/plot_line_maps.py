"""
Moment maps and P-V cuts for any of the ALMA lines observed toward IR 23128.

This generalises what plot_co_kinematics.py did for CO(3-2). Pick a line with
`--line` and it writes, into <plot-root>/radio_maps/:

    alma_<line>_moment0_<palette>.png   integrated intensity
    alma_<line>_moment1_<palette>.png   intensity-weighted velocity
    alma_<line>_moment2_<palette>.png   intensity-weighted velocity dispersion
    alma_<line>_pv_<palette>.png        a P-V cut through each nucleus

All four come out of one pass, because the expensive step is the mask, not
the arithmetic.

Choosing a line
---------------
LINES below records, for each, the cube that holds it and its rest frequency.
That pairing is not guessable: ALMA spectral-window numbers do not run in
frequency order (in one of these datasets spw 23/25/27/29 sit at 318.6, 330.6,
328.7 and 316.9 GHz), and the archive's frequency_support field is sorted by
frequency rather than by window, so the two cannot be matched from metadata
alone. Each entry here was fixed by reading the cube's own header.

Rest frequencies are constants here rather than header values on purpose. The
header RESTFRQ of the CO(3-2) cube is 345.5 GHz, which is not CO(3-2) at
345.7960 GHz but the nominal tuning value that redshifts to the window centre;
trusting it puts every velocity 257 km/s out.

Why the mask is built on the cube
---------------------------------
Moment 0 tolerates being computed first and blanked afterwards, because noise
averages down in a sum. Moment 1 is a ratio and does not: a spaxel whose
denominator lands near zero returns an arbitrary velocity, and negative noise
in the numerator throws velocities outside the window entirely.

The mask is the usual dilated kind -- smooth, cut a secure core at `--core`
sigma, then grow that core through everything above `--wing` sigma, so faint
line wings survive without a low threshold dragging in noise. Noise is
measured per channel, since window edges are markedly noisier than the middle,
and on the flat-noise cube (pbcor * pb) so the mask does not follow the
primary-beam correction outward.

Read moment 2 as an upper limit: it cannot separate turbulence from a velocity
gradient the beam fails to resolve.
"""
import argparse
import os
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Ellipse

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, AsinhStretch
from astropy.visualization.wcsaxes import add_scalebar
from astropy.nddata import Cutout2D
from reproject import reproject_interp
from scipy import ndimage

from irspec.plotparams import PlotParams
from irspec import paths

warnings.filterwarnings("ignore", category=FITSFixedWarning)

REDSHIFT = 0.044601                  # IR 23128
LUMINOSITY_DISTANCE = 194.99 * u.Mpc
ANGULAR_DIAMETER_DISTANCE = LUMINOSITY_DISTANCE / (1 + REDSHIFT) ** 2

# Data and output locations resolve lazily through irspec.paths.
DEFAULT_SUBDIR = "radio_maps"


def _alma_dir():
    return str(paths.image_data_dir() / "alma")


def _nircam_file():
    return str(paths.image_data_dir() / "nircam" / "IR23128_f200w.fits")

LINES = {
    "co32": {
        "stem": ("member.uid___A001_X2ca_X36.ari_l.ESO_148-IG002_sci"
                 ".spw1_330633MHz.12m.cube.I"),
        "rest": 345.7960, "label": "CO(3-2)", "project": "2015.1.00102.S"},
    "co10": {
        "stem": "member.uid___A001_X879_X7a0.ESO148-ig002_sci.spw25.cube.I",
        "rest": 115.2712, "label": "CO(1-0)", "project": "2016.1.00140.S"},
    "co10-ext": {
        "stem": ("member.uid___A001_X1284_X223e.IRAS_23128-5919_sci"
                 ".spw25.cube.I"),
        "rest": 115.2712, "label": "CO(1-0)", "project": "2017.1.01398.S"},
    "hcn32": {
        "stem": "member.uid___A001_X1289_X55c.ESO148-IG002_sci.spw27.cube.I",
        "rest": 265.8864, "label": "HCN(3-2)", "project": "2017.1.00759.S"},
    "hcop32": {
        "stem": "member.uid___A001_X1289_X55c.ESO148-IG002_sci.spw25.cube.I",
        "rest": 267.5576, "label": "HCO+(3-2)", "project": "2017.1.00759.S"},
    "hcn43": {
        "stem": ("member.uid___A001_X5a3_X132.ari_l.ESO_148-IG002_sci"
                 ".spw1_339564MHz.12m.cube.I"),
        "rest": 354.5055, "label": "HCN(4-3)", "project": "2015.1.00102.S"},
    "hcop43": {
        "stem": ("member.uid___A001_X5a3_X132.ari_l.ESO_148-IG002_sci"
                 ".spw2_341325MHz.12m.cube.I"),
        "rest": 356.7342, "label": "HCO+(4-3)", "project": "2015.1.00102.S"},
    "cs76": {
        "stem": ("member.uid___A001_X2ca_X36.ari_l.ESO_148-IG002_sci"
                 ".spw2_328676MHz.12m.cube.I"),
        "rest": 342.8829, "label": "CS(7-6)", "project": "2015.1.00102.S"},
}

LIGHT_SPEED = 299792.458 * u.km / u.s

WINDOW = 400.0 * u.km / u.s          # velocity range considered
BASELINE = 500.0 * u.km / u.s        # beyond this defines the continuum
MIN_BASELINE_CHANNELS = 10

CROP_RADIUS = 8.0 * u.arcsec

# Mask smoothing, in beams and channels rather than pixels: these cubes run
# from 0.12" to 0.92" resolution over pixel scales that differ by a factor of
# several, so a fixed pixel count would smooth one cube by a fraction of its
# beam and another by many beams.
SMOOTH_BEAMS = 1.0
SMOOTH_SPECTRAL = 1.5                # channels, sigma

DEFAULT_CORE_SIGMA = 4.0
DEFAULT_WING_SIGMA = 2.0
MIN_CHANNELS = 2
DEFAULT_MOMENT0_SNR = 5.0

PEAK_EXCLUSION = 2.0 * u.arcsec
SEARCH_RADIUS = 6.0 * u.arcsec

# One systemic velocity for every line, so the velocity fields are directly
# comparable and a faint line does not set its own zero point. Flux-weighting
# each line separately looks principled but is not: HCN(3-2) is 27 times
# fainter than CO(3-2) in integrated intensity, and weighting its noisy
# outskirts puts its "systemic" at +87 km/s even though its line peaks within
# 10 km/s of CO's. This value is the flux-weighted mean of CO(3-2), which is
# the best detected line here. Override with --systemic.
SYSTEMIC = 2.9

GRADIENT_RADIUS = 1.0 * u.arcsec
PV_LENGTH = 3.0 * u.arcsec

SCALEBAR_CHOICES = (0.2, 0.5, 1, 2, 5, 10) * u.kpc
MOMENT0_CMAP = "magma"
MOMENT1_CMAP = "RdBu_r"
MOMENT2_CMAP = "viridis"
PV_CMAP = "magma"
CONTOUR_COLOR = "white"
ASINH_SOFTENING = 0.15
CONTOUR_SIGMAS = (3, 10, 30, 100, 300)

# Optional stellar-morphology underlay. NIRCam F200W spans several decades
# across these nuclei, so its contours are log-spaced between percentiles
# rather than evenly; even spacing puts every level on the two nuclei and
# none on the tidal structure between them.

# Fixed absolute surface brightnesses, in the image's own MJy/sr, rather than
# percentiles of whatever happens to fall in each field. The line cubes have
# different pixel scales and field sizes, so percentile levels came out
# different in every figure and the underlays could not be compared with each
# other. These half-decade steps run from the outer stellar envelope (1, near
# the 80th percentile over a 20" box) to the nuclei (300, the 99.98th).
NIRCAM_LEVELS = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0)


def velocity_axis(header, rest):
    """Velocity of each channel relative to the redshifted line, in km/s."""
    channels = np.arange(header["NAXIS3"])
    frequency = ((channels - (header["CRPIX3"] - 1)) * header["CDELT3"]
                 + header["CRVAL3"]) * u.Hz
    observed = rest * u.GHz / (1 + REDSHIFT)
    return (LIGHT_SPEED * (observed - frequency) / observed).to(u.km / u.s)


def load_cropped_cube(spec):
    """Continuum-subtracted cube and flat-noise cube around the phase centre."""
    stem = os.path.join(_alma_dir(), spec["stem"])
    hdu = fits.open(stem + ".pbcor.fits", memmap=True)[0]
    header = hdu.header
    velocity = velocity_axis(header, spec["rest"])

    scale = abs(header["CDELT1"]) * 3600.0
    half = int(round(CROP_RADIUS.value / scale))
    half = min(half, header["NAXIS1"] // 2, header["NAXIS2"] // 2)
    cx, cy = header["NAXIS1"] // 2, header["NAXIS2"] // 2
    xs, ys = slice(cx - half, cx + half), slice(cy - half, cy + half)

    cube = np.squeeze(hdu.data)[:, ys, xs].astype(np.float32)

    with fits.open(stem + ".pb.fits.gz") as pbhdu:
        pb = np.squeeze(pbhdu[0].data)[header["NAXIS3"] // 2, ys, xs]
    pb = pb.astype(np.float32)

    # Widen the baseline definition if this window is too narrow to spare
    # MIN_BASELINE_CHANNELS at the nominal cut.
    baseline = BASELINE
    while (np.abs(velocity) >= baseline).sum() < MIN_BASELINE_CHANNELS:
        baseline *= 0.8
        if baseline < WINDOW:
            raise SystemExit("window too narrow to define a continuum baseline")
    off = np.abs(velocity) >= baseline

    # A first-order baseline, not a constant. These windows are not flat: the
    # HCN(3-2) window rises by about 6 mJy/beam from one edge to the other,
    # which is larger than the line peak, so subtracting a single mean leaves
    # a ramp that drags the intensity-weighted velocity redward by ~90 km/s.
    # Fitting a slope over the line-free channels removes it.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        design = np.vstack([velocity.value[off], np.ones(off.sum())]).T
        flat_pixels = np.nan_to_num(cube[off]).reshape(off.sum(), -1)
        coefficients, *_ = np.linalg.lstsq(design, flat_pixels, rcond=None)
        model = (np.outer(velocity.value, coefficients[0])
                 + coefficients[1][None, :]).reshape(cube.shape)
        cube -= model.astype(np.float32)

    width = abs(velocity[1] - velocity[0])
    print(f"  {cube.shape[2]}x{cube.shape[1]} pixels at {scale:.3f}\"/pix, "
          f"{header['NAXIS3']} channels of {width:.1f}")
    print(f"  beam {header['BMAJ'] * 3600:.3f}\" x {header['BMIN'] * 3600:.3f}\", "
          f"baseline beyond {baseline:.0f} ({off.sum()} channels)")

    celestial = WCS(header).celestial.to_header()
    for key in ("BMAJ", "BMIN", "BPA"):
        celestial[key] = header[key]
    return cube, cube * pb, velocity, WCS(celestial)[ys, xs], celestial


def build_mask(flat, velocity, header, scale, core_sigma, wing_sigma):
    """Dilated mask: a core cut at core_sigma grown through wing_sigma."""
    inside = np.abs(velocity) <= WINDOW
    beam_pixels = header["BMAJ"] * 3600.0 / scale
    spatial = max(1.0, SMOOTH_BEAMS * beam_pixels / 2.355)
    smoothed = ndimage.gaussian_filter(
        np.nan_to_num(flat), (SMOOTH_SPECTRAL, spatial, spatial))

    sigma = np.array([sigma_clipped_stats(plane, sigma=3.0, maxiters=5)[2]
                      for plane in smoothed], dtype=np.float32)
    sigma[~np.isfinite(sigma) | (sigma <= 0)] = np.inf
    print(f"  smoothing by {spatial:.1f} px ({SMOOTH_BEAMS:g} beam); smoothed "
          f"RMS {np.nanmin(sigma) * 1e3:.2f}-{np.nanmax(sigma[np.isfinite(sigma)]) * 1e3:.2f} mJy/beam")

    sigma = sigma[:, None, None]
    core = (smoothed > core_sigma * sigma) & inside[:, None, None]
    wing = (smoothed > wing_sigma * sigma) & inside[:, None, None]
    mask = ndimage.binary_propagation(core, mask=wing)
    mask &= (mask.sum(axis=0) >= MIN_CHANNELS)[None, :, :]
    print(f"  mask holds {mask.sum()} voxels over "
          f"{(mask.sum(axis=0) > 0).sum()} spaxels")
    return mask


def moments(cube, flat, mask, velocity, channel_rms, snr):
    """Masked moments 0, 1 and 2, gated on the flat-noise signal-to-noise."""
    width = abs(velocity[1] - velocity[0]).value
    values = velocity.value[:, None, None]
    weights = np.where(mask, cube, 0.0)

    total = weights.sum(axis=0)
    moment0 = total * width
    counts = mask.sum(axis=0)

    # Gate on the uncorrected moment: the primary-beam correction inflates
    # map and noise together, so comparing a corrected moment against an
    # uncorrected sigma lets the field edges through.
    gate = np.where(mask, flat, 0.0).sum(axis=0) * width

    with np.errstate(invalid="ignore", divide="ignore"):
        moment1 = np.where(total > 0, (weights * values).sum(axis=0) / total,
                           np.nan)
        spread = np.where(total > 0,
                          (weights * (values - moment1) ** 2).sum(axis=0) / total,
                          np.nan)
        moment2 = np.sqrt(np.where(spread > 0, spread, np.nan))

    noise = channel_rms * width * np.sqrt(np.maximum(counts, 1))
    good = (counts >= MIN_CHANNELS) & (gate > snr * noise)
    print(f"  {good.sum()} spaxels above {snr:g} sigma in moment 0")
    blank = lambda a: np.where(good, a, np.nan)
    return blank(moment0), blank(moment1), blank(moment2), good.sum()


def find_nuclei(moment0, wcs):
    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    rows, cols = np.mgrid[:moment0.shape[0], :moment0.shape[1]]
    centre = ((moment0.shape[0] - 1) / 2, (moment0.shape[1] - 1) / 2)
    near = np.hypot(rows - centre[0], cols - centre[1]) * scale <= SEARCH_RADIUS.value
    searchable = np.where(near, moment0, np.nan)
    if not np.isfinite(searchable).any():
        return {}
    first = np.unravel_index(np.nanargmax(searchable), searchable.shape)
    radius = np.hypot(rows - first[0], cols - first[1]) * scale
    rest = np.where(radius < PEAK_EXCLUSION.value, np.nan, searchable)
    if not np.isfinite(rest).any():
        return {"primary": first}
    second = np.unravel_index(np.nanargmax(rest), rest.shape)
    order = (first, second) if first[0] < second[0] else (second, first)
    return {"south": order[0], "north": order[1]}


def kinematic_angle(moment1, wcs, pixel):
    """Position angle of the steepest velocity gradient around a position."""
    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    rows, cols = np.mgrid[:moment1.shape[0], :moment1.shape[1]]
    east = -(cols - pixel[1]) * scale      # east is -x on a north-up image
    north = (rows - pixel[0]) * scale

    near = (np.hypot(east, north) <= GRADIENT_RADIUS.value) & np.isfinite(moment1)
    if near.sum() < 20:
        return None, None
    design = np.column_stack([east[near], north[near], np.ones(near.sum())])
    (grad_east, grad_north, _), *_ = np.linalg.lstsq(design, moment1[near],
                                                     rcond=None)
    return (np.degrees(np.arctan2(grad_east, grad_north)) % 360.0,
            np.hypot(grad_east, grad_north))


def extract_pv(cube, wcs, header, pixel, angle):
    """Position-velocity slice along `angle` through `pixel`."""
    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    width = header["BMAJ"] * 3600.0     # slit one beam wide
    offsets = np.arange(-PV_LENGTH.value / 2, PV_LENGTH.value / 2 + scale, scale)
    across = np.arange(-width / 2, width / 2 + scale, scale)

    theta = np.radians(angle)
    east = offsets[:, None] * np.sin(theta) + across[None, :] * np.cos(theta)
    north = offsets[:, None] * np.cos(theta) - across[None, :] * np.sin(theta)
    coordinates = np.array([(pixel[0] + north / scale).ravel(),
                            (pixel[1] - east / scale).ravel()])

    slices = np.empty((cube.shape[0], offsets.size), dtype=np.float32)
    for index, plane in enumerate(cube):
        sampled = ndimage.map_coordinates(np.nan_to_num(plane), coordinates,
                                          order=1, mode="constant")
        slices[index] = sampled.reshape(east.shape).mean(axis=1)
    return offsets, slices


def choose_scalebar(frame_width):
    half = (frame_width.to(u.rad).value / 2
            * ANGULAR_DIAMETER_DISTANCE).to(u.kpc)
    fits_in = [c for c in SCALEBAR_CHOICES if c <= half]
    return max(fits_in) if fits_in else min(SCALEBAR_CHOICES)


def nircam_contours(ax, wcs, shape, pltparams):
    """Log-spaced NIRCam F200W contours reprojected onto the moment grid.

    The mosaic is cut down before reprojecting: it is 4815x4830 pixels at
    0.031", and reprojecting all of it onto a few hundred pixels would read
    most of a gigabyte to use a corner of it.
    """
    hdu = fits.open(_nircam_file(), memmap=True)[1]
    source = WCS(hdu.header).celestial
    centre = wcs.pixel_to_world((shape[1] - 1) / 2, (shape[0] - 1) / 2)
    span = shape[1] * abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    try:
        cut = Cutout2D(hdu.data, centre, 1.5 * span * u.arcsec, wcs=source)
    except Exception as error:
        print(f"  NIRCam underlay skipped: {error}")
        return
    data, _ = reproject_interp((cut.data, cut.wcs), wcs, shape_out=shape)

    if not np.isfinite(data).any():
        print("  NIRCam underlay skipped: no overlap")
        return
    # Drawn above the moment map, not below it: under a filled colormap the
    # contours were only visible in the blanked gaps, which is where they
    # matter least.
    ax.contour(data, levels=NIRCAM_LEVELS, colors=pltparams.foreground,
               linewidths=0.6, alpha=0.8, zorder=3)


def decorate(ax, field, header, wcs, pltparams, title):
    ax.coords[0].set_axislabel("RA (J2000)")
    ax.coords[1].set_axislabel("Dec (J2000)")
    ax.coords[0].set_major_formatter("hh:mm:ss.s")
    ax.coords[1].set_major_formatter("dd:mm:ss")
    ax.set_title(title, loc="right")

    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.deg).value)
    inset = 0.06 * field.shape[0]
    ax.add_patch(Ellipse((inset, inset), width=header["BMIN"] / scale,
                         height=header["BMAJ"] / scale,
                         angle=header["BPA"] + 90,
                         facecolor=pltparams.foreground,
                         edgecolor=pltparams.foreground, alpha=0.8))
    length = choose_scalebar(field.shape[1] * scale * 3600 * u.arcsec)
    add_scalebar(ax, (length / ANGULAR_DIAMETER_DISTANCE).to(
                     u.deg, equivalencies=u.dimensionless_angles()),
                 label=f"{length:latex_inline}", color=pltparams.foreground,
                 fontproperties=fm.FontProperties(
                     size=pltparams.annotation_size, family=pltparams.font))


def render_map(field, moment0, wcs, header, pltparams, plot_root, subdir,
               name, title, label, cmap, vmin, vmax, contours=True,
               underlay=False):
    fig = plt.figure(figsize=(9, 8))
    ax = plt.subplot(projection=wcs)
    ax.set_facecolor(pltparams.background)

    colormap = pltparams.colormap(cmap)
    if underlay:
        # Blanked spaxels normally take the figure background, which is
        # opaque and would hide the underlay; make them transparent so the
        # stellar contours show through where there is no line emission.
        colormap = colormap.copy()
        colormap.set_bad((0, 0, 0, 0))
        nircam_contours(ax, wcs, field.shape, pltparams)

    image = ax.imshow(field, origin="lower", cmap=colormap,
                      vmin=vmin, vmax=vmax, zorder=2)
    decorate(ax, field, header, wcs, pltparams, title)

    if contours and np.isfinite(moment0).any():
        levels = np.nanmax(moment0) * np.array([0.05, 0.2, 0.5, 0.8])
        if np.all(np.diff(levels) > 0):
            ax.contour(moment0, levels=levels, colors=pltparams.foreground,
                       linewidths=0.5, alpha=0.6)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    colorbar.ax.tick_params(labelsize=pltparams.tick_size)
    pltparams.label_colorbar(colorbar, label)

    os.makedirs(os.path.join(plot_root, subdir), exist_ok=True)
    out = os.path.join(plot_root, subdir,
                       f"{name}{pltparams.file_suffix}.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def render_pv(cuts, velocity, systemic, rms, label, pltparams, plot_root,
              subdir, name):
    # Frequency rises with channel, so velocity falls with it; imshow with
    # origin="lower" would read the rows as ascending and mirror the axis.
    relative = velocity.value - systemic
    order = np.argsort(relative)
    relative = relative[order]

    fig, axes = plt.subplots(1, len(cuts), figsize=(6 * len(cuts), 6),
                             sharey=True, squeeze=False)
    for ax, (nucleus, offsets, slices, angle) in zip(axes[0], cuts):
        slices = slices[order]
        extent = [offsets[0], offsets[-1], relative[0], relative[-1]]
        ax.imshow(slices * 1e3, origin="lower", aspect="auto", extent=extent,
                  cmap=pltparams.colormap(PV_CMAP),
                  vmin=-2 * rms * 1e3, vmax=max(np.nanmax(slices) * 1e3, 1e-6))
        ax.contour(slices * 1e3, levels=np.array([3, 6, 12, 24]) * rms * 1e3,
                   colors="white", linewidths=0.6, alpha=0.6, extent=extent,
                   origin="lower")
        ax.axhline(0.0, color="white", lw=0.6, ls=":", alpha=0.7)
        ax.axvline(0.0, color="white", lw=0.6, ls=":", alpha=0.7)
        ax.set_ylim(-WINDOW.value, WINDOW.value)
        ax.set_xlabel("offset (arcsec)")
        ax.set_title(f"{label}, {nucleus} nucleus, PA {angle:.0f}", loc="right")
    axes[0][0].set_ylabel(r"$v - v_{\rm sys}$ (km s$^{-1}$)")
    fig.tight_layout()

    os.makedirs(os.path.join(plot_root, subdir), exist_ok=True)
    out = os.path.join(plot_root, subdir, f"{name}{pltparams.file_suffix}.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def build_parser():
    parser = argparse.ArgumentParser(
        description="Moment maps and P-V cuts for an ALMA line of IR 23128.")
    parser.add_argument("--line", choices=sorted(LINES), default="co32")
    parser.add_argument("--palette",
                        choices=list(PlotParams.PALATTES) + ["both"],
                        default=PlotParams.DEFAULT_PALATTE)
    parser.add_argument("--scaling", choices=list(PlotParams.SCALINGS),
                        default="paper")
    parser.add_argument("--core", type=float, default=DEFAULT_CORE_SIGMA)
    parser.add_argument("--wing", type=float, default=DEFAULT_WING_SIGMA)
    parser.add_argument("--snr", type=float, default=DEFAULT_MOMENT0_SNR)
    parser.add_argument("--systemic", type=float, default=None)
    parser.add_argument("--underlay", choices=["none", "nircam", "both"],
                        default="none",
                        help="Contour NIRCam F200W beneath the moment maps "
                             "and tag the filenames (default: none)")
    parser.add_argument("--plot-root", default=None)
    parser.add_argument("--subdir", default=DEFAULT_SUBDIR)
    return parser


def main():
    args = build_parser().parse_args()
    if args.plot_root is None:
        args.plot_root = str(paths.plots_dir())
    spec = LINES[args.line]
    label = spec["label"]
    print(f"{label} from {spec['project']}")

    cube, flat, velocity, wcs, header = load_cropped_cube(spec)
    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    mask = build_mask(flat, velocity, header, scale, args.core, args.wing)

    channel_rms = float(np.median(
        [sigma_clipped_stats(p, sigma=3.0, maxiters=5)[2] for p in flat]))
    moment0, moment1, moment2, kept = moments(cube, flat, mask, velocity,
                                              channel_rms, args.snr)
    if kept == 0:
        print(f"\n{label} is not detected above the threshold; no maps written")
        return

    systemic = SYSTEMIC if args.systemic is None else args.systemic
    print(f"  systemic {systemic:+.1f} km/s (shared), peak moment 0 "
          f"{np.nanmax(moment0):.2f} Jy/beam km/s")

    # "both" renders the plain and underlaid versions from this one masked
    # pass, which is the expensive part; rendering is not. The same applies to
    # the palettes, so a single run can produce all four combinations.
    modes = {"none": [("", False)], "nircam": [("_nircam", True)],
             "both": [("", False), ("_nircam", True)]}[args.underlay]
    palettes = (list(PlotParams.PALATTES) if args.palette == "both"
                else [args.palette])

    # Colour ranges come from the bright emission only. A percentile over
    # every surviving spaxel is set by the faint clumps at the edge of the
    # field, whose velocities are the least reliable in the map, and that
    # stretches the scale until the nuclei are all one colour.
    bright = moment0 > 0.2 * np.nanmax(moment0)
    field = moment1 - systemic
    scale_from = field[bright] if bright.sum() > 10 else field
    limit = float(np.nanpercentile(np.abs(scale_from), 98))
    dispersion_max = float(np.nanpercentile(
        moment2[bright] if bright.sum() > 10 else moment2, 98))

    # Kinematics are palette-independent, so they are worked out once.
    cuts = []
    for nucleus, pixel in find_nuclei(moment0, wcs).items():
        angle, gradient = kinematic_angle(moment1, wcs, pixel)
        if angle is None:
            print(f"  {nucleus}: too few valid spaxels to fit a gradient")
            continue
        coord = wcs.pixel_to_world(pixel[1], pixel[0])
        print(f"  {nucleus} at {coord.to_string('hmsdms', precision=2)}: "
              f"v = {moment1[pixel] - systemic:+.1f} km/s, "
              f"sigma = {moment2[pixel]:.0f} km/s, PA {angle:.1f}, "
              f"gradient {gradient:.0f} km/s/arcsec")
        offsets, slices = extract_pv(cube, wcs, header, pixel, angle)
        cuts.append((nucleus, offsets, slices, angle))

    written = []
    for palette in palettes:
        pltparams = PlotParams(palatte=palette, scaling=args.scaling)
        for tag, under in modes:
            written.append(render_map(
                moment0, moment0, wcs, header, pltparams, args.plot_root,
                args.subdir, f"alma_{args.line}_moment0{tag}",
                f"{label} integrated intensity",
                r"$\int S_\nu\,dv$  (Jy beam$^{-1}$ km s$^{-1}$)",
                MOMENT0_CMAP, 0.0, float(np.nanmax(moment0)),
                contours=False, underlay=under))
            written.append(render_map(
                field, moment0, wcs, header, pltparams, args.plot_root,
                args.subdir, f"alma_{args.line}_moment1{tag}",
                f"{label} velocity field",
                r"$v - v_{\rm sys}$ (km s$^{-1}$)", MOMENT1_CMAP,
                -limit, limit, underlay=under))
            written.append(render_map(
                moment2, moment0, wcs, header, pltparams, args.plot_root,
                args.subdir, f"alma_{args.line}_moment2{tag}",
                f"{label} velocity dispersion",
                r"$\sigma_v$ (km s$^{-1}$)", MOMENT2_CMAP, 0.0,
                dispersion_max, underlay=under))
        if cuts:
            written.append(render_pv(cuts, velocity, systemic, channel_rms,
                                     label, pltparams, args.plot_root,
                                     args.subdir, f"alma_{args.line}_pv"))
    for path in written:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
