"""
Maps the millimetre spectral index of the IR 23128 merger from two ALMA bands.

The two continuum maps are

    band 6   2017.1.00759.S   247.6 GHz   0.298" x 0.239"   MRS 3.95"
    band 7   2015.1.00102.S   323.8 GHz   0.294" x 0.254"   MRS 3.9"

which is an unusually clean pair to difference: the beams differ by a few
hundredths of an arcsecond and the maximum recoverable scales agree to within
a twentieth of one, so neither map is filtering out structure the other keeps.
A spectral index built from mismatched MRS is meaningless no matter how
carefully the beams are matched, and that is the usual reason these maps go
wrong.

The index is S propto nu^alpha, so for optically thin dust alpha = 2 + beta
and the dust emissivity index follows directly as beta = alpha - 2. Optically
thick dust sits at alpha = 2, optically thin free-free near -0.1 and
synchrotron near -0.7, so the map doubles as a check that these bands are
tracing dust at all rather than a mix.

The figure lands in <plot-root>/images/ as

    alma_spectral_index_<palette>.png

Usage examples (from anywhere; the script anchors its own paths):
    python plot_spectral_index.py
    python plot_spectral_index.py --palette dark --snr 8

Read the index with care
------------------------
Two things measured from these data say the index here is indicative, not
publishable, and both are worth knowing before the colours are believed.

First, the band 7 *aggregate* continuum is line contaminated. At this redshift
CO(3-2) falls at 331.03 GHz, inside the 330.6 GHz spectral window, and that
window is measurably more extended than the line-free ones (a factor 1.8
between r = 0.5" and 2", against 1.2-1.4 where no line falls). Differencing
against the aggregate raises alpha by roughly 1.4 at every aperture. Hence
the default here is the line-free 316.9 GHz window rather than the aggregate,
at the cost of a quarter of the bandwidth and so a noisier map.

Second, and not fixed by any choice of window, alpha depends on the aperture
it is measured in -- about +1.4 at r = 0.5" rising to +2.8 at r = 1.5" for the
line-free window. A real spectral index does not do that. It means the two
datasets are not recovering the same spatial scales despite their maximum
recoverable scales agreeing on paper, so the ratio picks up the difference in
what each observation filtered out. Fixing that needs the two re-imaged from
the visibilities on a common uv range, which is beyond what these delivered
products can support.

Matching the resolutions
------------------------
Rather than deconvolving either beam to a common one, each map is convolved
with the *other* map's beam. Both then carry an effective beam of B6 (*) B7 --
identical by construction, no deconvolution, no assumption that one beam
cleanly contains the other. It costs a little resolution (the common beam is
larger than either input) and that is the whole price.

The convolution is done on maps in Jy/pixel, not Jy/beam: a unit-area kernel
conserves total flux in the former and does not in the latter, and the ratio
of two Jy/beam maps is only the ratio of flux densities once their beams
already agree.
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
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.visualization.wcsaxes import add_scalebar
from reproject import reproject_interp

from irspec.plotparams import PlotParams
from irspec import paths

warnings.filterwarnings("ignore", category=FITSFixedWarning)

REDSHIFT = 0.044601                  # IR 23128
LUMINOSITY_DISTANCE = 194.99 * u.Mpc # IR 23128

# An angular scale wants the angular-diameter distance, not the luminosity
# distance; D_A = D_L/(1+z)^2 is the Etherington reciprocity relation, which
# holds in any FLRW metric, so deriving it keeps this tied to whatever
# cosmology produced the 194.99 Mpc.
ANGULAR_DIAMETER_DISTANCE = LUMINOSITY_DISTANCE / (1 + REDSHIFT) ** 2

# Data and output locations resolve lazily through irspec.paths.
def _alma_dir():
    return str(paths.image_data_dir() / "alma")
BAND6_STEM = "member.uid___A001_X1289_X55c.ESO148-IG002_sci.spw25_27_29_31.cont.I"

# Two choices of band 7 image. "linefree" is one spectral window clear of every
# bright line in the tuning, and is the default for the reason set out at the
# top; "aggregate" is the pipeline's four-window continuum, which is deeper but
# carries CO(3-2). The line-free window's beam also happens to sit within
# 0.02" of the band 6 beam, so it needs the least matching of the two.
BAND7_STEMS = {
# NOTE: "line-free" overstates it. 13CO(3-2) (rest 330.5880 GHz) lands at
# 316.47 GHz, inside this window's 315.93-317.80 GHz coverage. 13CO is far
# weaker than 12CO -- this window shows the *least* extended emission of the
# four, so the contamination is evidently small -- but the window is
# line-poor rather than line-free.
    "linefree": ("member.uid___A001_X2ca_X36.ari_l.ESO_148-IG002_sci"
                 ".spw3_316867MHz.12m.mfs.I"),
    "aggregate": ("member.uid___A001_X2ca_X36.ESO_148-IG002_sci"
                  ".spw23_25_27_29.cont.I"),
}
DEFAULT_BAND7 = "linefree"

SCALEBAR_LENGTH = 2 * u.kpc

# Frame width as a multiple of the nuclear separation. Tighter than the single
# band map: outside the nuclei there is nothing detected in both bands, so a
# wider frame is entirely masked.
DEFAULT_FRAME_FACTOR = 0.8

# Signal-to-noise required in *both* convolved maps before an index is
# computed for a pixel. A spectral index is a ratio, so its error runs away
# as either denominator approaches the noise; 5 keeps the formal uncertainty
# near +/-0.3 and blanks the rest rather than drawing colour where there is
# no measurement.
DEFAULT_SNR = 5.0

# Displayed range. 2 is optically thick dust, 4 is optically thin dust with
# beta = 2, so this brackets the physically interesting span; values outside
# are clipped by the colormap rather than rescaling everything to an outlier.
ALPHA_RANGE = (2.0, 4.5)

DISPLAY_PB_FLOOR = 0.3
PEAK_EXCLUSION = 2.0 * u.arcsec
ON_IMAGE_COLOR = "black"
COLORMAP = "Spectral_r"


def load_band(stem, label):
    """Load one band's continuum map, its beam and its noise level.

    The pipeline images carry four axes (RA, Dec, FREQ, STOKES) and only the
    primary-beam-corrected one ships, so the degenerate axes are squeezed, the
    WCS is reduced to its celestial part, and the flat-noise map whose RMS is
    quoted is rebuilt as pbcor * pb.
    """
    stem = os.path.join(_alma_dir(), stem)
    hdu = fits.open(stem + ".pbcor.fits")[0]
    image = np.squeeze(hdu.data).astype(float)
    pb = np.squeeze(fits.open(stem + ".pb.fits.gz")[0].data)

    flat = image * pb
    _, _, rms = sigma_clipped_stats(flat[np.isfinite(flat)], sigma=3.0,
                                    maxiters=5)
    image = np.where(pb > DISPLAY_PB_FLOOR, image, np.nan)

    header = WCS(hdu.header).celestial.to_header()
    for key in ("BMAJ", "BMIN", "BPA"):
        header[key] = hdu.header[key]
    frequency = hdu.header["CRVAL3"] * u.Hz

    print(f"  {label}: {frequency.to(u.GHz):.2f}, beam "
          f"{hdu.header['BMAJ'] * 3600:.3f}\" x {hdu.header['BMIN'] * 3600:.3f}\" "
          f"PA {hdu.header['BPA']:.1f}, RMS {rms * 1e6:.1f} uJy/beam")
    return {"image": image, "flat": flat, "wcs": WCS(header),
            "header": header, "frequency": frequency, "rms": rms}


def find_nuclei(band_map):
    """Locate the two nuclei as the two brightest separated peaks.

    The search runs on the flat-noise map rather than the primary-beam
    corrected one: the correction inflates the noise by 1/pb towards the
    field edge, so on the corrected map an edge fluctuation can outrank a
    real source.
    """
    image, wcs = band_map["flat"], band_map["wcs"]
    first = np.unravel_index(np.nanargmax(image), image.shape)
    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    rows, cols = np.mgrid[:image.shape[0], :image.shape[1]]
    masked = np.where(np.hypot(rows - first[0], cols - first[1]) * scale
                      < PEAK_EXCLUSION.value, np.nan, image)
    second = np.unravel_index(np.nanargmax(masked), masked.shape)
    return [wcs.pixel_to_world(p[1], p[0]) for p in (first, second)]


def beam_covariance(header, pixel_scale_deg):
    """Covariance matrix of a beam, in pixels, from its FITS keywords.

    Position angle is measured from north toward east. North is +y and east
    is -x on these north-up images, so that rotation runs counterclockwise
    from +y, putting the major axis at 90 + BPA from the +x axis.
    """
    to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    major = header["BMAJ"] / pixel_scale_deg * to_sigma
    minor = header["BMIN"] / pixel_scale_deg * to_sigma
    theta = np.radians(90.0 + header["BPA"])

    rotation = np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta), np.cos(theta)]])
    return rotation @ np.diag([major ** 2, minor ** 2]) @ rotation.T


def covariance_to_beam(covariance):
    """Major/minor sigma and position angle of a beam from its covariance."""
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    major, minor = np.sqrt(values)
    theta = np.arctan2(vectors[1, 0], vectors[0, 0])
    return major, minor, np.degrees(theta) - 90.0


def kernel_from_header(header, pixel_scale_deg):
    """A unit-area Gaussian kernel matching one map's beam."""
    covariance = beam_covariance(header, pixel_scale_deg)
    major, minor, position_angle = covariance_to_beam(covariance)
    return Gaussian2DKernel(x_stddev=major, y_stddev=minor,
                            theta=np.radians(90.0 + position_angle))


def match_and_ratio(low, high, frame_factor, snr):
    """Put both bands on one grid and effective beam, then take the index.

    Returns the index map, the convolved low-band map for contouring, the
    cutout carrying the shared WCS, and the common beam.
    """
    coords = find_nuclei(low)
    separation = coords[0].separation(coords[1]).to(u.arcsec)
    centre = coords[0].directional_offset_by(
        coords[0].position_angle(coords[1]), separation / 2)
    projected = (separation.to(u.rad).value
                 * ANGULAR_DIAMETER_DISTANCE).to(u.kpc)
    print(f"  nuclei separated by {separation.value:.2f}\" "
          f"({projected.value:.2f} kpc projected)")

    # The low-resolution band defines the grid; the other is interpolated on
    # to it. Both are Jy/beam, which interpolation leaves alone.
    cutout = Cutout2D(low["image"], centre, separation * 2 * frame_factor,
                      wcs=low["wcs"])
    high_on_grid, _ = reproject_interp(
        (high["image"], high["wcs"]), cutout.wcs, shape_out=cutout.shape)

    pixel_scale = abs(cutout.wcs.proj_plane_pixel_scales()[0].to(u.deg).value)
    to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    maps, beams = {}, {}
    for label, band, data in (("low", low, cutout.data),
                              ("high", high, high_on_grid)):
        header = band["header"]
        beams[label] = beam_covariance(header, pixel_scale)
        # Jy/beam -> Jy/pixel, so that a unit-area kernel conserves flux.
        beam_pixels = (2.0 * np.pi * np.sqrt(np.linalg.det(beams[label])))
        maps[label] = (data / beam_pixels, beam_pixels)

    # Each map convolved with the other's beam: both end on B_low (*) B_high.
    convolved = {}
    for label, other in (("low", "high"), ("high", "low")):
        kernel = kernel_from_header(
            (low if other == "low" else high)["header"], pixel_scale)
        convolved[label] = convolve_fft(
            maps[label][0], kernel, normalize_kernel=True,
            nan_treatment="fill", fill_value=0.0, allow_huge=True)

    common = beams["low"] + beams["high"]
    major, minor, position_angle = covariance_to_beam(common)
    print(f"  common beam {major / to_sigma * pixel_scale * 3600:.3f}\" x "
          f"{minor / to_sigma * pixel_scale * 3600:.3f}\" "
          f"PA {position_angle:.1f}")

    # Convolution correlates the noise, so the input RMS no longer applies;
    # measure it again on each convolved map, where sigma clipping rejects
    # the two nuclei and leaves the sky.
    index = np.full(convolved["low"].shape, np.nan)
    detected = np.ones(convolved["low"].shape, dtype=bool)
    for label in ("low", "high"):
        _, _, rms = sigma_clipped_stats(convolved[label], sigma=3.0,
                                        maxiters=5)
        detected &= convolved[label] > snr * rms
        print(f"  {label}-band convolved RMS {rms * 1e6:.2f} uJy/pixel, "
              f"{(convolved[label] > snr * rms).sum()} pixels above {snr:g} sigma")

    ratio = np.log(convolved["high"][detected] / convolved["low"][detected])
    index[detected] = ratio / np.log(high["frequency"] / low["frequency"])
    print(f"  index measured in {detected.sum()} pixels")
    report_aperture_dependence(convolved, cutout, coords,
                               np.log(high["frequency"] / low["frequency"]))
    return index, convolved["low"], cutout, (major, minor, position_angle)


def report_aperture_dependence(convolved, cutout, coords, log_frequency_ratio):
    """Print alpha integrated in a few apertures around each nucleus.

    A spectral index that shifts with the aperture it is measured in is not
    measuring a spectrum; it is measuring the difference in what the two
    observations filtered out. Printing it every run keeps that visible
    instead of leaving it to be discovered later.
    """
    scale = abs(cutout.wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    rows, cols = np.mgrid[:cutout.data.shape[0], :cutout.data.shape[1]]
    for name, coord in zip(("south", "north"), coords):
        x, y = cutout.wcs.world_to_pixel(coord)
        radius = np.hypot(rows - float(y), cols - float(x)) * scale
        values = []
        for aperture in (0.5, 1.0, 1.5):
            inside = radius < aperture
            low_flux = np.nansum(convolved["low"][inside])
            high_flux = np.nansum(convolved["high"][inside])
            if low_flux <= 0 or high_flux <= 0:
                values.append(f'r={aperture}": n/a')
                continue
            alpha = np.log(high_flux / low_flux) / log_frequency_ratio
            values.append(f'r={aperture}": {alpha:+.2f}')
        print(f"  {name} nucleus alpha by aperture -- {', '.join(values)}")


def render(index, low_convolved, cutout, common_beam, low, high,
           pltparams, plot_root):
    """Draw the index map with low-band contours over it."""
    valid = np.isfinite(index)
    print(f"  alpha: median {np.nanmedian(index):.2f}, "
          f"range {np.nanmin(index):.2f} to {np.nanmax(index):.2f}")

    fig = plt.figure(figsize=(9, 8))
    ax = plt.subplot(projection=cutout.wcs)
    ax.set_facecolor(pltparams.background)
    image = ax.imshow(index, origin="lower", cmap=pltparams.colormap(COLORMAP),
                      vmin=ALPHA_RANGE[0], vmax=ALPHA_RANGE[1])

    ax.coords[0].set_axislabel("RA (J2000)")
    ax.coords[1].set_axislabel("Dec (J2000)")
    ax.coords[0].set_major_formatter("hh:mm:ss.s")
    ax.coords[1].set_major_formatter("dd:mm:ss")
    ax.set_title(f"{low['frequency'].to(u.GHz).value:.0f} / "
                 f"{high['frequency'].to(u.GHz).value:.0f} GHz spectral index",
                 loc="right")

    # Contours of the lower band locate the nuclei under the index colours.
    _, _, rms = sigma_clipped_stats(low_convolved, sigma=3.0, maxiters=5)
    ax.contour(low_convolved, levels=np.array([5, 20, 80]) * rms,
               colors=ON_IMAGE_COLOR, linewidths=0.6, alpha=0.5)

    pixel_scale = abs(cutout.wcs.proj_plane_pixel_scales()[0].to(u.deg).value)
    to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    inset = 0.08 * cutout.data.shape[0]
    ax.add_patch(Ellipse((inset, inset),
                         width=common_beam[1] / to_sigma,
                         height=common_beam[0] / to_sigma,
                         angle=common_beam[2] + 90,
                         facecolor=ON_IMAGE_COLOR, edgecolor=ON_IMAGE_COLOR,
                         alpha=0.8))

    scalebar_angle = (SCALEBAR_LENGTH / ANGULAR_DIAMETER_DISTANCE).to(
        u.deg, equivalencies=u.dimensionless_angles())
    add_scalebar(ax, scalebar_angle, label=f"{SCALEBAR_LENGTH:latex_inline}",
                 color=ON_IMAGE_COLOR,
                 fontproperties=fm.FontProperties(
                     size=pltparams.annotation_size, family=pltparams.font))

    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    colorbar.ax.tick_params(labelsize=pltparams.tick_size)
    pltparams.label_colorbar(colorbar, r"$\alpha$  ($S_\nu \propto \nu^\alpha$)")

    os.makedirs(os.path.join(plot_root, "images"), exist_ok=True)
    out = os.path.join(plot_root, "images",
                       f"alma_spectral_index{pltparams.file_suffix}.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def build_parser():
    parser = argparse.ArgumentParser(
        description="Map the ALMA band 6 / band 7 spectral index of IR 23128.")
    parser.add_argument("--palette", choices=list(PlotParams.PALATTES),
                        default=PlotParams.DEFAULT_PALATTE,
                        help="Figure palette; tags the filename "
                             f"(default: {PlotParams.DEFAULT_PALATTE})")
    parser.add_argument("--scaling", choices=list(PlotParams.SCALINGS),
                        default="paper",
                        help="Font/DPI scaling (default: paper)")
    parser.add_argument("--band7", choices=sorted(BAND7_STEMS),
                        default=DEFAULT_BAND7,
                        help="Which band 7 image to difference against; "
                             "'aggregate' is deeper but carries CO(3-2) "
                             f"(default: {DEFAULT_BAND7})")
    parser.add_argument("--snr", type=float, default=DEFAULT_SNR,
                        help="Signal-to-noise required in both bands "
                             f"(default: {DEFAULT_SNR})")
    parser.add_argument("--frame-factor", type=float,
                        default=DEFAULT_FRAME_FACTOR,
                        help="Frame width as a multiple of the nuclear "
                             f"separation (default: {DEFAULT_FRAME_FACTOR})")
    parser.add_argument("--plot-root", default=None,
                        help="Root for rendered figures; this goes to "
                             "<plot-root>/images/ "
                             "(default: plots/ under the irspec output root)")
    return parser


def main():
    args = build_parser().parse_args()
    if args.plot_root is None:
        args.plot_root = str(paths.plots_dir())
    pltparams = PlotParams(palatte=args.palette, scaling=args.scaling)

    low = load_band(BAND6_STEM, "band 6")
    high = load_band(BAND7_STEMS[args.band7], f"band 7 ({args.band7})")
    if args.band7 == "aggregate":
        print("  WARNING: the aggregate band 7 continuum carries CO(3-2); "
              "alpha runs about 1.4 high. See the module docstring.")
    index, low_convolved, cutout, common_beam = match_and_ratio(
        low, high, args.frame_factor, args.snr)
    out = render(index, low_convolved, cutout, common_beam, low, high,
                 pltparams, args.plot_root)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
