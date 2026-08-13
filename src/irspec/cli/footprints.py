"""
Overlays the IFU footprints on the NIRCam and MIRI imaging.

For every broad-band image under irspec/image_data/, this draws the sky
footprint of each MRS channel and of the NIRSpec cube on top of it, so the
spectroscopic coverage can be read against the morphology the imaging
shows. A channel is represented by its short sub-band; see SUBBANDS for
why one stands in for all three.

Two versions of every image are rendered:

    flux     the image itself, log-stretched, with the footprints over it
    contour  log-spaced intensity contours of the same image, no greyscale

Both land in <plot-root>/footprints/ as

    <instrument>_<filter>_<pointing>_<style>_<palette>.png

e.g. plots/footprints/nircam_f200w_south_flux_light.png

Usage examples (from anywhere; the script anchors its own paths):
    python footprints.py
    python footprints.py --pointing north
    python footprints.py --pointing both --palette dark
    python footprints.py --style contour --scaling presentation

Which footprint is drawn
------------------------
The footprint of a cube is the set of spaxels that carry a spectrum --
`np.any(np.isfinite(...))` down the wavelength axis -- outlined along the
spaxel edges. It is *not* the S_REGION polygon in the cube header: these
cubes are resampled onto a north-up grid, and their S_REGION is the
axis-aligned bounding box of that grid, which is a good deal larger than
the rotated square the dither pattern actually covers (only about half the
array of a typical MRS cube is measured at all).

The outline is built by the same code the per-spaxel maps use for their own
footprint (`Spaxelcube._footprint_edges`), so the boundary drawn here is the
same boundary drawn there, by construction rather than by agreement.

Segmentation maps (*_segm.fits) in the MIRI directory are skipped: they hold
integer source labels, not flux, and nothing about a log stretch or a set of
intensity contours means anything on them.
"""
import argparse
import glob
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.visualization.wcsaxes import add_scalebar

from irspec.datacube import Datacube
from irspec.plotparams import PlotParams
from irspec import paths
from irspec.spaxel_fit import Spaxelcube, contour_field

warnings.filterwarnings("ignore", category=FITSFixedWarning)

DEFAULT_REDSHIFT = 0.044601          # IR 23128
LUMINOSITY_DISTANCE = 194.99 * u.Mpc # IR 23128

# An angular scale wants the angular-diameter distance, not the luminosity
# distance. Deriving one from the other rather than from a cosmology keeps
# this tied to whatever cosmology produced the 194.99 Mpc: D_A = D_L/(1+z)^2
# is the Etherington reciprocity relation and holds in any FLRW metric.
#
# D_A is the smaller of the two (178.7 Mpc here), so a given physical length
# subtends a *larger* angle than the luminosity distance implies: bars drawn
# from D_L come out about 9% too short for their labels, and a size read off
# one is overestimated by that much.
DISTANCE = LUMINOSITY_DISTANCE / (1 + DEFAULT_REDSHIFT) ** 2

# Data and output locations resolve lazily through irspec.paths
# (IRSPEC_DATA_ROOT / IRSPEC_OUTPUT_ROOT).
def image_dirs():
    root = paths.image_data_dir()
    return {"NIRCam": str(root / "nircam"), "MIRI": str(root / "miri")}


def data_dir_for(pointing):
    return str(paths.cube_dir(pointing))


NIRSPEC_CUBE = "g395h-f290lp_s3d.fits"

CHANNELS = (1, 2, 3, 4)
# One sub-band per channel. The three sub-bands of a channel are the same
# pointing read out over different wavelengths, so their footprints sit on
# top of each other to within a spaxel or two -- drawing all three triples
# the lines on the panel and the entries in the legend without showing
# anything the one does not. Widen to ("short", "medium", "long") to draw
# them all; the dash patterns below are what keeps that case legible.
SUBBANDS = ("short",)

# One hue per MRS channel, plus one for NIRSpec, and one dash pattern per
# sub-band.
CHANNEL_COLORS = {1: "#0bb4ff", 2: "#00bfa0", 3: "#ffa300", 4: "#e60049"}
NIRSPEC_COLOR = "#9b19f5"
SUBBAND_STYLES = {"short": "solid", "medium": (0, (6, 2)),
                  "long": (0, (1.5, 1.5))}

# Greyscale for the log-stretched image, per palette: bright source on a
# light page, bright source on a dark one.
IMAGE_CMAPS = {"light": "bone_r", "dark": "bone"}
# Sequential map the contour levels are coloured through. Not keyed to the
# palette: the contour version draws no image, so the levels need a map that
# reads against bare background either way, which viridis does.
CONTOUR_CMAP = "viridis"

# Blank sky left around the outermost footprint, as a fraction of the span
# the footprints themselves cover.
VIEW_MARGIN = 0.35
# Number of log-spaced contour levels, and the percentile of the in-view
# pixels the lowest one sits at. Starting well up the distribution keeps the
# faint end -- which is noise at this depth -- from filling the panel with
# contours that trace the read pattern rather than the galaxy.
CONTOUR_LEVELS = 8
CONTOUR_FLOOR_PERCENTILE = 92.0
# Percentiles setting the ends of the log stretch on the flux version.
FLUX_PERCENTILES = (40.0, 99.9)


def cube_specs(pointing):
    """The cubes to outline for one pointing, in drawing order.

    Returns a list of dicts with the cube's `label`, `path`, `color` and
    `linestyle`. Cubes that are not on disk are dropped with a note: the
    north pointing has the twelve MRS bands but no NIRSpec cube, and that
    is a fact about the observations, not an error.
    """
    data_dir = data_dir_for(pointing)
    wanted = [(f"CH{channel}-{subband}",
               f"ch{channel}-{subband}_s3d.fits",
               CHANNEL_COLORS[channel],
               SUBBAND_STYLES[subband])
              for channel in CHANNELS for subband in SUBBANDS]
    wanted.append(("NIRSpec G395H", NIRSPEC_CUBE, NIRSPEC_COLOR, "solid"))

    specs = []
    for label, filename, color, linestyle in wanted:
        path = os.path.join(data_dir, filename)
        if not os.path.exists(path):
            print(f"  no {label} cube for the {pointing} pointing ({path})")
            continue
        specs.append({"label": label, "path": path,
                      "color": color, "linestyle": linestyle})
    return specs


def footprint_edges_sky(cube_path, redshift):
    """The footprint outline of one cube, as sky-coordinate line segments.

    Returns a SkyCoord of shape (n_segments, 2): the two endpoints of each
    edge. Working in sky coordinates rather than pixels is what lets one
    footprint be drawn onto any image, whatever grid that image is on.
    """
    cube = Datacube(cube_path, redshift=redshift)
    # Every spaxel holding at least one finite channel. See the module
    # docstring for why this, and not the header's S_REGION.
    measured = np.any(np.isfinite(cube.science_data), axis=0)
    # The same edge-walk the per-spaxel maps outline themselves with. It is
    # reached through the class rather than copied so the two cannot drift.
    segments = np.asarray(Spaxelcube._footprint_edges(measured), dtype=float)
    if segments.size == 0:
        return None
    # (n, 2, 2) -> (2n, 2) so the whole outline transforms in one call.
    corners = segments.reshape(-1, 2)
    sky = cube.wcs.pixel_to_world(corners[:, 0], corners[:, 1])
    return sky.reshape(segments.shape[:-1])


def to_image_segments(sky_edges, image_wcs):
    """Projects sky-coordinate footprint edges into an image's pixel grid."""
    x, y = image_wcs.world_to_pixel(sky_edges.ravel())
    return np.stack([x, y], axis=1).reshape(-1, 2, 2)


# INSTRUME as the pipeline writes it, cased as it is written in prose.
INSTRUMENT_NAMES = {"NIRCAM": "NIRCam", "MIRI": "MIRI", "NIRSPEC": "NIRSpec"}


def load_image(path):
    """Reads the science array, WCS, unit, and identifying labels of an image."""
    with fits.open(path) as hdul:
        raw = hdul[0].header.get("INSTRUME", "")
        instrument = INSTRUMENT_NAMES.get(raw.upper(), raw.title())
        band = hdul[0].header.get("FILTER", "")
        sci = hdul["SCI"]
        # Left at its stored dtype: these mosaics run to a couple of hundred
        # megabytes as float64, and nothing here needs the extra precision.
        data = sci.data
        wcs = WCS(sci.header).celestial
        unit = sci.header.get("BUNIT", "")
    return data, wcs, instrument, band, unit


def image_files():
    """Every broad-band image to plot, as (path, instrument, filter).

    Only the headers are read here; the arrays are loaded one at a time in
    the render loop. Segmentation maps are excluded -- they carry source
    labels rather than flux, so neither rendering means anything on them.
    """
    found = []
    for directory in image_dirs().values():
        for path in sorted(glob.glob(os.path.join(directory, "*.fits"))):
            if path.endswith("_segm.fits"):
                continue
            found.append(path)
    return found


def view_limits(segment_sets, data_shape):
    """Frames the view on the footprints, with a margin of sky around them.

    The images are whole-mosaic products tens of arcseconds across and the
    footprints span a few, so shown whole the overlay would be a speck.
    The box is squared off so the two renderings of an image, and the
    different images, are all directly comparable.
    """
    points = np.concatenate([segments.reshape(-1, 2)
                             for segments in segment_sets])
    (x_min, y_min), (x_max, y_max) = points.min(axis=0), points.max(axis=0)

    half = 0.5 * (1 + 2 * VIEW_MARGIN) * max(x_max - x_min, y_max - y_min)
    x_mid, y_mid = 0.5 * (x_min + x_max), 0.5 * (y_min + y_max)
    rows, cols = data_shape
    # Clipped to the array: asking for sky the image does not cover would
    # pad the panel with blank margin and shrink the part that has data.
    return ((max(x_mid - half, -0.5), min(x_mid + half, cols - 0.5)),
            (max(y_mid - half, -0.5), min(y_mid + half, rows - 0.5)))


def view_slice(shape, xlim, ylim):
    """The array slice covering the framed box, as (y0, y1, x0, x1)."""
    rows, cols = shape
    x_lo, x_hi = (int(np.floor(v + 0.5)) for v in xlim)
    y_lo, y_hi = (int(np.floor(v + 0.5)) for v in ylim)
    return (max(y_lo, 0), min(y_hi + 1, rows),
            max(x_lo, 0), min(x_hi + 1, cols))


def in_view(data, xlim, ylim):
    """The finite, positive pixels inside the framed box.

    Stretches and contour levels are set from these rather than from the
    whole mosaic: the panel shows one galaxy out of a wide field, and a
    stretch fitted to the field is set by sky that is never displayed.
    """
    y0, y1, x0, x1 = view_slice(data.shape, xlim, ylim)
    window = data[y0:y1, x0:x1]
    values = window[np.isfinite(window) & (window > 0)]
    if values.size == 0:
        # Nothing positive in frame -- an all-negative or all-blank cutout.
        # Fall back to the whole array so the stretch is still defined, and
        # let the caller say so rather than dying on an empty percentile.
        values = data[np.isfinite(data) & (data > 0)]
    return values


def draw_flux(fig, ax, data, values, pltparams, unit):
    """Draws the image itself, log-stretched, and its colorbar."""
    low, high = np.percentile(values, FLUX_PERCENTILES)
    image = ax.imshow(data, origin="lower",
                      norm=LogNorm(vmin=low, vmax=high),
                      cmap=pltparams.colormap(IMAGE_CMAPS[pltparams.palatte]))
    pltparams.label_colorbar(fig.colorbar(image, ax=ax, pad=0.02), unit)


def draw_contours(fig, ax, data, xlim, ylim, pltparams, unit):
    """Draws log-spaced intensity contours, and a colorbar keyed to them.

    Contours only the framed box, passing explicit pixel coordinates so the
    cropped array still lands in the right place. Cropping first is what
    makes `contour_field` cheap: it runs a hole-fill over the array it is
    given, and these mosaics are tens of millions of pixels of which the
    panel shows a few hundred thousand.
    """
    y0, y1, x0, x1 = view_slice(data.shape, xlim, ylim)
    # `contour_field` masks the non-positive pixels a log norm cannot hold,
    # then closes the level sets over the saturated cores and caps the
    # levels below them -- see its docstring for why both are needed.
    field, levels = contour_field(data[y0:y1, x0:x1], CONTOUR_LEVELS,
                                  CONTOUR_FLOOR_PERCENTILE)
    if field is None:
        return
    norm = LogNorm(vmin=levels[0], vmax=levels[-1])
    cmap = CONTOUR_CMAP
    yy, xx = np.mgrid[y0:y1, x0:x1]
    ax.contour(xx, yy, field, levels=levels, norm=norm,
               cmap=cmap, linewidths=1.2)

    # The bar is built from a ScalarMappable rather than from the ContourSet.
    # A colorbar taken straight off contour *lines* paints only the lines,
    # leaving a near-empty white strip, and labels every level to full
    # precision. This draws the ramp the levels are coloured from, ticked at
    # the levels themselves.
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                       pad=0.02, ticks=levels)
    bar.ax.set_yticklabels([_level_label(level) for level in levels])
    bar.ax.minorticks_off()
    pltparams.label_colorbar(bar, unit)


def _level_label(level):
    """Formats a contour level for the colorbar: no more digits than the
    value warrants, and no exponent until the number is genuinely long."""
    if level >= 1000:
        return f"{level:.0f}"
    if level >= 10:
        return f"{level:.0f}"
    if level >= 1:
        return f"{level:.1f}"
    return f"{level:.2f}"


# Where the compass sits in the panel, and how long its arms are, as
# fractions of the framed box.
_COMPASS_ANCHOR = (0.13, 0.84)
_COMPASS_LENGTH = 0.08
_COMPASS_LABEL_OFFSET = 1.35


def add_compass(ax, wcs, xlim, ylim, pltparams):
    """Draws N and E arrows, oriented by the image's own WCS."""
    x_span, y_span = xlim[1] - xlim[0], ylim[1] - ylim[0]
    x_anchor = xlim[0] + _COMPASS_ANCHOR[0] * x_span
    y_anchor = ylim[0] + _COMPASS_ANCHOR[1] * y_span
    anchor = wcs.pixel_to_world(x_anchor, y_anchor)
    length = _COMPASS_LENGTH * min(x_span, y_span)

    color = pltparams.foreground
    halo = [pe.withStroke(linewidth=3, foreground=pltparams.background)]
    for position_angle, label in ((0 * u.deg, "N"), (90 * u.deg, "E")):
        tip = wcs.world_to_pixel(
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
                                    path_effects=halo), zorder=6)
        ax.text(x_anchor + dx * _COMPASS_LABEL_OFFSET,
                y_anchor + dy * _COMPASS_LABEL_OFFSET, label, color=color,
                fontsize=pltparams.annotation_size, ha="center", va="center",
                path_effects=halo, zorder=6)


def render(data, wcs, instrument, band, unit, pointing, specs, segment_sets,
           style, pltparams, plot_root):
    """Renders one image in one style and writes it out. Returns the path."""
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(projection=wcs)

    xlim, ylim = view_limits(segment_sets, data.shape)
    # Set before the image: `set_xlim` turns autoscaling off, so `imshow`
    # inherits the framing instead of snapping back to the full mosaic.
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    label = f"Surface Brightness [{unit}]" if unit else "Surface Brightness"
    if style == "flux":
        draw_flux(fig, ax, data, in_view(data, xlim, ylim), pltparams, label)
    else:
        draw_contours(fig, ax, data, xlim, ylim, pltparams, label)

    for spec, segments in zip(specs, segment_sets):
        # autolim=False: the framing is already decided above, and an
        # outline should not be able to argue with it.
        ax.add_collection(
            LineCollection(segments, colors=spec["color"],
                           linestyles=spec["linestyle"], linewidths=1.6,
                           label=spec["label"], zorder=5),
            autolim=False)

    scalebar_angle = (1 * u.kpc / DISTANCE).to(
        u.deg, equivalencies=u.dimensionless_angles())
    fontprops = fm.FontProperties(size=pltparams.annotation_size,
                                  family=pltparams.font)
    add_scalebar(ax, scalebar_angle, label="1 kpc",
                 color=pltparams.foreground, fontproperties=fontprops)
    add_compass(ax, wcs, xlim, ylim, pltparams)

    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.set_title(f"{instrument} {band}", loc="left")
    ax.set_title(f"{pointing} nucleus", loc="right")
    # Below the panel rather than inside it, where it would sit on the very
    # footprints it is labelling. The size comes from `legend.fontsize`,
    # which the scaling has already set on the rcParams.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=5,
              frameon=False, handlelength=2.6, columnspacing=1.4)

    os.makedirs(os.path.join(plot_root, "footprints"), exist_ok=True)
    name = (f"{instrument}_{band}_{pointing}_{style}"
            f"{pltparams.file_suffix}.png").lower()
    out = os.path.join(plot_root, "footprints", name)
    fig.savefig(out)
    plt.close(fig)
    return out


def build_parser():
    """Builds the argument parser."""
    parser = argparse.ArgumentParser(
        description="Overlay the MRS and NIRSpec IFU footprints on the "
                    "NIRCam and MIRI imaging.")
    parser.add_argument("--pointing", choices=["south", "north", "both"],
                        default="south",
                        help="Which nucleus's cubes to outline. The north "
                             "pointing has the twelve MRS bands but no "
                             "NIRSpec cube (default: south)")
    parser.add_argument("--style", choices=["flux", "contour", "both"],
                        default="both",
                        help="Render the log-stretched image, the intensity "
                             "contours, or both (default: both)")
    parser.add_argument("--redshift", type=float, default=DEFAULT_REDSHIFT,
                        help=f"Source redshift (default: {DEFAULT_REDSHIFT})")
    parser.add_argument("--palette", choices=list(PlotParams.PALATTES),
                        default=PlotParams.DEFAULT_PALATTE,
                        help="Figure palette; tags the filename "
                             f"(default: {PlotParams.DEFAULT_PALATTE})")
    parser.add_argument("--scaling", choices=list(PlotParams.SCALINGS),
                        default="paper",
                        help="Font/DPI scaling (default: paper)")
    parser.add_argument("--plot-root", default=None,
                        help="Root for rendered figures; these go to "
                             "<plot-root>/footprints/ "
                             "(default: plots/ under the irspec output root)")
    return parser


def main():
    args = build_parser().parse_args()
    if args.plot_root is None:
        args.plot_root = str(paths.plots_dir())
    pltparams = PlotParams(palatte=args.palette, scaling=args.scaling)

    pointings = (["south", "north"] if args.pointing == "both"
                 else [args.pointing])
    styles = (["flux", "contour"] if args.style == "both" else [args.style])

    images = image_files()
    if not images:
        sys.exit("No images found under " + ", ".join(image_dirs().values()))

    # Each cube is read once and reduced to its sky outline here, before any
    # image is opened: the cubes are the expensive part to read, and an
    # outline on the sky does not depend on what it will be drawn over.
    outlines = {}
    for pointing in pointings:
        print(f"{pointing} pointing:")
        kept = []
        for spec in cube_specs(pointing):
            edges = footprint_edges_sky(spec["path"], args.redshift)
            if edges is None:
                print(f"  {spec['label']} has no measured spaxels; skipping")
                continue
            kept.append((spec, edges))
        if not kept:
            print("  no cubes to outline; skipping")
            continue
        outlines[pointing] = kept
        print(f"  outlined {len(kept)} cubes")

    written = []
    for path in images:
        data, wcs, instrument, band, unit = load_image(path)
        for pointing, kept in outlines.items():
            segment_sets = [to_image_segments(edges, wcs)
                            for _, edges in kept]
            specs = [spec for spec, _ in kept]
            for style in styles:
                out = render(data, wcs, instrument, band, unit, pointing,
                             specs, segment_sets, style, pltparams,
                             args.plot_root)
                written.append(out)
                print(f"  {out}")

    print(f"\nWrote {len(written)} figure(s).")


if __name__ == "__main__":
    main()
