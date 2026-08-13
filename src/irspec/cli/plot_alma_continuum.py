"""
Renders the ALMA 1.2mm continuum map of the IR 23128 merger.

The map is the aggregate continuum of project 2017.1.00759.S: four spectral
windows combined at 247.6 GHz, restored with a 0.298" x 0.239" beam. Both
nuclei of the merger fall inside the one 40" field, so a single cutout centred
between them frames the whole system, and the framing is taken from the data
rather than hardcoded -- the two brightest separated peaks are located, the
cutout is centred on their midpoint and sized from their separation.

The figure lands in <plot-root>/images/ as

    alma_1p2mm_continuum_<palette>.png

Usage examples (from anywhere; the script anchors its own paths):
    python plot_alma_continuum.py
    python plot_alma_continuum.py --palette dark
    python plot_alma_continuum.py --softening 0.02 --frame-factor 1.6

Dynamic range
-------------
The map runs from a 20 uJy/beam noise floor to a 1.8 mJy/beam peak, so the
full span is about 90:1 and a linear scale set by the peak leaves everything
below ~10 sigma black. The stretch is asinh rather than log: it stays linear
through the noise, so the background keeps its texture, and compresses only
once the signal is well clear of it.

`--softening` is the real control, and it trades one thing against the other.
Smaller values compress the bright end harder, lifting fainter emission but
lifting the noise along with it. At the 0.08 default the 1 sigma noise lands
near 8% of the colour range -- present as texture, not competing with the
sources. At 0.02 it lands near 27% and the map reads as pure speckle.
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
from astropy.visualization import ImageNormalize, AsinhStretch
from astropy.visualization.wcsaxes import add_scalebar

from irspec.plotparams import PlotParams
from irspec import paths

warnings.filterwarnings("ignore", category=FITSFixedWarning)

REDSHIFT = 0.044601                  # IR 23128
LUMINOSITY_DISTANCE = 194.99 * u.Mpc # IR 23128, as quoted elsewhere here

# An angular scale wants the angular-diameter distance, not the luminosity
# distance the rest of the repository passes to add_scalebar. Deriving one
# from the other rather than from a cosmology keeps this figure tied to
# whatever cosmology produced the 194.99 Mpc: D_A = D_L / (1+z)^2 is the
# Etherington reciprocity relation and holds in any FLRW metric.
#
# D_A is the smaller of the two (178.7 Mpc here), so a given physical length
# subtends a *larger* angle than the luminosity distance implies -- the bars
# drawn from D_L elsewhere in the repository are about 9% too short for their
# labels, and a size read off one of them is overestimated by that much.
ANGULAR_DIAMETER_DISTANCE = LUMINOSITY_DISTANCE / (1 + REDSHIFT) ** 2

# Data and output locations resolve lazily through irspec.paths.
def _alma_dir():
    return str(paths.image_data_dir() / "alma")

# The continuum maps this can draw. b7 is the pipeline's four-window aggregate,
# which is the deepest band 7 image but carries CO(3-2) at 331.0 GHz; b7-linefree
# is the one window clear of every bright line in that tuning, a quarter of the
# bandwidth and correspondingly noisier. For looking at the morphology either
# is fine, but the contaminated one is not a pure dust map -- see
# plot_spectral_index.py, which measures how far off it runs.
#
# The suffixes differ between deliveries. Older cycles ship a single-term
# image as <stem>.pbcor.fits with <stem>.pb.fits.gz beside it; the band 4 data
# were cleaned with multi-term mfs, so the continuum is the zeroth Taylor term,
# <stem>.tt0.pbcor.fits, and the primary beam is <stem>.pb.tt0.fits. Carrying
# both per band beats guessing from the filename.
BAND_MAPS = {
    "b4": {"stem": ("member.uid___A001_X3788_X7bda.F23128-5919_sci"
                    ".spw17_19_21_23.cont.I"),
           "image": ".tt0.pbcor.fits", "pb": ".pb.tt0.fits"},
    "b6": {"stem": ("member.uid___A001_X1289_X55c.ESO148-IG002_sci"
                    ".spw25_27_29_31.cont.I"),
           "image": ".pbcor.fits", "pb": ".pb.fits.gz"},
    "b7": {"stem": ("member.uid___A001_X2ca_X36.ESO_148-IG002_sci"
                    ".spw23_25_27_29.cont.I"),
           "image": ".pbcor.fits", "pb": ".pb.fits.gz"},
# NOTE: "line-free" overstates it. 13CO(3-2) (rest 330.5880 GHz) lands at
# 316.47 GHz, inside this window's 315.93-317.80 GHz coverage. 13CO is far
# weaker than 12CO -- this window shows the *least* extended emission of the
# four, so the contamination is evidently small -- but the window is
# line-poor rather than line-free.
    "b7-linefree": {"stem": ("member.uid___A001_X2ca_X36.ari_l"
                             ".ESO_148-IG002_sci.spw3_316867MHz.12m.mfs.I"),
                    "image": ".pbcor.fits", "pb": ".pb.fits.gz"},
}
DEFAULT_BAND = "b6"
DEFAULT_SUBDIR = "images"

# Round lengths the scalebar may take. The frame varies from 6" when only one
# nucleus is detected to 12" when both are, and a bar fixed at 5 kpc spans
# almost the entire short frame, so the longest of these that still fits in
# half the frame is chosen instead.
SCALEBAR_CHOICES = (0.2, 0.5, 1, 2, 5, 10) * u.kpc

# How much sky to leave around the two nuclei, as a multiple of their
# separation. 1.1 seats both well inside the frame while letting the system
# fill it; larger values pad the figure out with blank noise.
DEFAULT_FRAME_FACTOR = 1.1

# Asinh softening, as a fraction of the displayed range. See the module
# docstring for what this trades.
DEFAULT_SOFTENING = 0.08

# Where the colour scale starts, in units of the map RMS. Slightly negative so
# the noise reads as noise -- a floor at zero clips the negative half of it and
# the background turns into a flat, false-looking wall of colormap zero.
FLOOR_SIGMA = -1.0

# Contour levels, in units of the RMS. These pick the faint structure back out
# of the compressed bright end, where the stretch has flattened the contrast.
CONTOUR_SIGMAS = (5, 10, 20, 40, 80)

# Primary-beam level below which pixels are blanked for display. Lower than
# the 0.5 used when measuring, purely so the cutout corners stay filled: the
# half-power radius falls just inside the framed region and would otherwise
# cut a visible arc across it.
DISPLAY_PB_FLOOR = 0.3

# Radius around the brightest peak excluded when looking for the second one,
# so the search does not return two pixels of the same nucleus.
PEAK_EXCLUSION = 2.0 * u.arcsec

# ...and the radius beyond which a peak is not taken to be the companion at
# all. The two nuclei are 5.6" apart, so anything past this is a different
# object. Without it the band 4 map, where the north nucleus is largely
# resolved out by a 0.87" maximum recoverable scale, pairs the south nucleus
# with an unrelated source 15" away and frames a field containing neither
# nucleus at its centre.
PEAK_MAX_SEPARATION = 10.0 * u.arcsec

# Significance a second peak must clear before it is treated as the companion
# nucleus rather than a noise excursion. The real companion is 25 sigma in
# band 6 and 13 in band 7; in band 4 the north nucleus is largely resolved out
# and the best candidate is 6 sigma, so that map is framed on the south
# nucleus alone instead of on a pair that is not there.
COMPANION_SNR = 10.0

# Frame width used when only one nucleus is detected, where there is no
# separation to scale from.
SINGLE_FRAME_SIZE = 6.0 * u.arcsec

# Colour for everything drawn *on* the image -- contours, the beam, the
# scalebar. Deliberately not pltparams.foreground: that contrasts against the
# figure background, which is white under the light palette, whereas these
# annotations sit on the near-black low end of the colormap in either palette.
ON_IMAGE_COLOR = "white"

COLORMAP = "inferno"


def load_continuum(band):
    """Load one band's continuum map, its WCS, frequency and noise level.

    The pipeline image carries four axes (RA, Dec, FREQ, STOKES), so the
    degenerate axes are squeezed away and the WCS reduced to its celestial
    part. Only the primary-beam-corrected image ships in this cycle's product
    set, and its noise climbs toward the field edge, so the flat-noise map
    whose RMS is quoted has to be rebuilt as pbcor * pb.
    """
    spec = BAND_MAPS[band]
    stem = os.path.join(_alma_dir(), spec["stem"])
    hdu = fits.open(stem + spec["image"])[0]
    image = np.squeeze(hdu.data)
    pb = np.squeeze(fits.open(stem + spec["pb"])[0].data)

    flat_noise = image * pb
    _, _, rms = sigma_clipped_stats(flat_noise[np.isfinite(flat_noise)],
                                    sigma=3.0, maxiters=5)

    # Blank the far field. Out there the correction amplifies noise faster
    # than signal, and those pixels would otherwise set the display maximum
    # and flatten the whole stretch. The RMS above is measured first, on the
    # flat-noise map, so this only affects what is drawn.
    image = np.where(pb > DISPLAY_PB_FLOOR, image, np.nan)

    header = WCS(hdu.header).celestial.to_header()
    for key in ("BMAJ", "BMIN", "BPA"):
        header[key] = hdu.header[key]
    frequency = hdu.header["CRVAL3"] * u.Hz
    return image, flat_noise, WCS(header), header, rms, frequency


def find_nuclei(flat_noise, image, wcs, rms):
    """Locate the two merger nuclei as the two brightest separated peaks.

    The search runs on the flat-noise map, not the primary-beam-corrected one.
    Dividing by the primary beam inflates the noise towards the field edge by
    1/pb, so on the corrected map an edge fluctuation can outshine a real
    source: on the band 7 image the second peak found that way is a ~4 sigma
    blob at pb ~ 0.3, fifteen arcseconds off, rather than the north nucleus.
    Brightness is still quoted from the corrected map, which is the calibrated
    one, and at the nuclei pb is essentially unity so the two agree.
    """
    first = np.unravel_index(np.nanargmax(flat_noise), flat_noise.shape)

    scale = abs(wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value)
    rows, cols = np.mgrid[:flat_noise.shape[0], :flat_noise.shape[1]]
    radius = np.hypot(rows - first[0], cols - first[1]) * scale
    masked = np.where((radius < PEAK_EXCLUSION.value)
                      | (radius > PEAK_MAX_SEPARATION.value),
                      np.nan, flat_noise)
    second = np.unravel_index(np.nanargmax(masked), masked.shape)

    peaks = [first]
    if flat_noise[second] >= COMPANION_SNR * rms:
        peaks.append(second)
    else:
        print(f"  no companion above {COMPANION_SNR:g} sigma within "
              f"{PEAK_MAX_SEPARATION}; framing on the one detection "
              f"(best candidate {flat_noise[second] / rms:.0f} sigma)")

    coords = [wcs.pixel_to_world(peak[1], peak[0]) for peak in peaks]
    for coord, peak in zip(coords, [image[p] for p in peaks]):
        print(f"  nucleus at {coord.to_string('hmsdms', precision=2)}: "
              f"{peak * 1e3:.3f} mJy/beam ({peak / rms:.0f} sigma)")
    return coords


def frame_on_system(image, wcs, coords, frame_factor):
    """Cut out a square centred between the nuclei and sized to hold both."""
    if len(coords) == 1:
        print(f"  framing {SINGLE_FRAME_SIZE} on "
              f"{coords[0].to_string('hmsdms', precision=2)}")
        return Cutout2D(image, coords[0], SINGLE_FRAME_SIZE, wcs=wcs)

    separation = coords[0].separation(coords[1]).to(u.arcsec)
    centre = coords[0].directional_offset_by(
        coords[0].position_angle(coords[1]), separation / 2)

    projected = (separation.to(u.rad).value
                 * ANGULAR_DIAMETER_DISTANCE).to(u.kpc)

    size = separation * 2 * frame_factor
    print(f"  nuclei separated by {separation.value:.2f}\" "
          f"({projected.value:.2f} kpc projected); framing "
          f"{size.value:.1f}\" on {centre.to_string('hmsdms', precision=2)}")
    return Cutout2D(image, centre, size, wcs=wcs)


def render(cutout, header, rms, softening, pltparams, plot_root, subdir,
           band, frequency):
    """Draw the cutout with an asinh stretch, beam, contours and scalebar."""
    data = cutout.data * 1e3          # Jy/beam -> mJy/beam
    rms_mjy = rms * 1e3
    peak = np.nanmax(data)

    norm = ImageNormalize(data, vmin=FLOOR_SIGMA * rms_mjy, vmax=peak,
                          stretch=AsinhStretch(softening))
    print(f"  displaying {FLOOR_SIGMA:+.0f} sigma to peak, "
          f"dynamic range {peak / rms_mjy:.0f}:1")

    fig = plt.figure(figsize=(9, 8))
    ax = plt.subplot(projection=cutout.wcs)
    image = ax.imshow(data, origin="lower", norm=norm,
                      cmap=pltparams.colormap(COLORMAP))

    ax.coords[0].set_axislabel("RA (J2000)")
    ax.coords[1].set_axislabel("Dec (J2000)")
    ax.coords[0].set_major_formatter("hh:mm:ss.s")
    ax.coords[1].set_major_formatter("dd:mm:ss")
    wavelength = frequency.to(u.mm, equivalencies=u.spectral())
    ax.set_title(f"ALMA {wavelength.value:.2f} mm continuum "
                 f"({frequency.to(u.GHz).value:.0f} GHz)", loc="right")

    ax.contour(data, levels=np.array(CONTOUR_SIGMAS) * rms_mjy,
               colors=ON_IMAGE_COLOR, linewidths=0.5, alpha=0.4)

    draw_beam(ax, cutout, header)

    frame = (cutout.data.shape[1]
             * abs(cutout.wcs.proj_plane_pixel_scales()[0].to(u.arcsec)))
    length = choose_scalebar(frame)
    scalebar_angle = (length / ANGULAR_DIAMETER_DISTANCE).to(
        u.deg, equivalencies=u.dimensionless_angles())
    add_scalebar(ax, scalebar_angle, label=f"{length:latex_inline}",
                 color=ON_IMAGE_COLOR,
                 fontproperties=fm.FontProperties(
                     size=pltparams.annotation_size, family=pltparams.font))

    # Ticks at sigma multiples rather than round flux values: on an asinh
    # scale the eye cannot read significance off the colours, and significance
    # is what the contour levels mean.
    ticks = np.array((0,) + CONTOUR_SIGMAS) * rms_mjy
    ticks = ticks[ticks <= peak]
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02,
                            ticks=ticks)
    colorbar.ax.set_yticklabels([f"{tick:.2f}" for tick in ticks],
                                fontsize=pltparams.tick_size)
    pltparams.label_colorbar(colorbar, r"$S_\nu$ (mJy beam$^{-1}$)")

    os.makedirs(os.path.join(plot_root, subdir), exist_ok=True)
    out = os.path.join(plot_root, subdir,
                       f"alma_{band}_continuum{pltparams.file_suffix}.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def choose_scalebar(frame_width):
    """Longest round scalebar that still fits in half the frame."""
    half = (frame_width.to(u.rad).value / 2
            * ANGULAR_DIAMETER_DISTANCE).to(u.kpc)
    fits_in = [c for c in SCALEBAR_CHOICES if c <= half]
    return max(fits_in) if fits_in else min(SCALEBAR_CHOICES)


def draw_beam(ax, cutout, header):
    """Draw the restoring beam in the lower-left corner, as radio maps do."""
    scale = abs(cutout.wcs.proj_plane_pixel_scales()[0].to(u.deg).value)

    # Position angle is measured from north toward east; matplotlib measures
    # counterclockwise from the x-axis, hence the offset.
    inset = 0.08 * cutout.data.shape[0]
    ax.add_patch(Ellipse((inset, inset),
                         width=header["BMIN"] / scale,
                         height=header["BMAJ"] / scale,
                         angle=header["BPA"] + 90,
                         facecolor=ON_IMAGE_COLOR,
                         edgecolor=ON_IMAGE_COLOR, alpha=0.8))


def build_parser():
    parser = argparse.ArgumentParser(
        description="Render the ALMA 1.2mm continuum map of IR 23128.")
    parser.add_argument("--band", choices=sorted(BAND_MAPS),
                        default=DEFAULT_BAND,
                        help="Which continuum map to draw; 'b7' is the "
                             "aggregate band 7, which carries CO(3-2) "
                             f"(default: {DEFAULT_BAND})")
    parser.add_argument("--subdir", default=DEFAULT_SUBDIR,
                        help="Subdirectory of <plot-root> to write into "
                             f"(default: {DEFAULT_SUBDIR})")
    parser.add_argument("--palette", choices=list(PlotParams.PALATTES),
                        default=PlotParams.DEFAULT_PALATTE,
                        help="Figure palette; tags the filename "
                             f"(default: {PlotParams.DEFAULT_PALATTE})")
    parser.add_argument("--scaling", choices=list(PlotParams.SCALINGS),
                        default="paper",
                        help="Font/DPI scaling (default: paper)")
    parser.add_argument("--frame-factor", type=float,
                        default=DEFAULT_FRAME_FACTOR,
                        help="Frame width as a multiple of the nuclear "
                             f"separation (default: {DEFAULT_FRAME_FACTOR})")
    parser.add_argument("--softening", type=float, default=DEFAULT_SOFTENING,
                        help="Asinh softening; smaller lifts fainter "
                             "emission and more noise with it "
                             f"(default: {DEFAULT_SOFTENING})")
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

    image, flat_noise, wcs, header, rms, frequency = load_continuum(args.band)
    print(f"  {args.band}: {frequency.to(u.GHz):.2f}, RMS {rms * 1e6:.1f} "
          f"uJy/beam, beam {header['BMAJ'] * 3600:.3f}\" x "
          f"{header['BMIN'] * 3600:.3f}\"")

    coords = find_nuclei(flat_noise, image, wcs, rms)
    cutout = frame_on_system(image, wcs, coords, args.frame_factor)

    out = render(cutout, header, rms, args.softening, pltparams,
                 args.plot_root, args.subdir, args.band, frequency)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
