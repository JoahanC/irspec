"""
Fits spaxel maps for one or more emission lines and renders the standard
set of diagnostic maps (component count, total flux, relative velocity,
and velocity dispersion) for each.

Data products and figures go to two separate trees.

Fit tables, FITS parameter cubes, and q3dfit's own products are written to
<output-root>/<line>/ for the south source and <output-root>/north/<line>/
for the north source, with the output root defaulting to the
repository-level spaxs/ directory.

Rendered maps are written to <plot-root>/line_maps/<line>/<run-name>/,
with the plot root defaulting to the repository-level plots/ directory.
The run name is chosen per invocation (--run-name), so several fitting
runs of the same line -- different fitters, parameters, or sources -- keep
their figures side by side instead of overwriting each other. It defaults
to "<source>_<fitter>" (e.g. "south_twogauss").

Figures render in the light palette by default, or the dark one with
--palette dark. The palette tags the filename (..._light.png /
..._dark.png), so rendering a line both ways into the same run directory
keeps both rather than overwriting.

Usage examples (from anywhere; the script anchors its own paths):
    python generate_spaxel_map.py "[NeIII]"
    python generate_spaxel_map.py "[NeII]" "[NeIII]" "[SIV]" --workers 8
    python generate_spaxel_map.py "[KVII]" --instrument NIRSpec
    python generate_spaxel_map.py "[NeIII]" --skip-fit   # re-render only
    python generate_spaxel_map.py "[SIV]" --fitter q3dfit --workers 8
    python generate_spaxel_map.py "[NeIII]" --run-name snr5_alpha01
    python generate_spaxel_map.py "[NeIII]" --skip-fit --palette dark
    python generate_spaxel_map.py "[NeIII]" --skip-fit --contour-image f200w
    python generate_spaxel_map.py --config runs.ini      # batch, see below

The default two-Gaussian fitter writes the fit table and rendered maps
described above. The alternative q3dfit backend (--fitter q3dfit) instead
writes q3dfit's own per-spaxel/collected products into the same per-line
data directory, and renders the analogous maps (component count, total
flux, and per-component velocity offset and dispersion) into the same
per-run plot directory.


Batch runs (--config)
---------------------
Command-line options apply to every line named on the command line, which
is limiting when a batch needs different settings per line. `--config`
takes an INI file (the same format as the CAFE parfiles) listing the lines
to run, each with its own settings, and works through them in file order:

    # runs.ini -- comments start with #
    [DEFAULT]                 # optional; inherited by every entry below
    source  = south
    workers = 8
    palette = light

    [[NeIII]]                 # section name is the emission line
    scaling = paper

    [[SIV]]                   # inherits every DEFAULT above
    fitter  = q3dfit
    workers = 4

    [h2s1_dark]               # a label, because it sets `line` explicitly
    line     = [H_2_S_1]
    palette  = dark
    skip-fit = true

Keys are the long option names without the leading dashes (`run-name` and
`run_name` are both accepted), and values are validated exactly as they
are on the command line, so an unknown fitter or a bad palette is rejected
with the same message. Boolean flags (`skip-fit`, `no-render`) take
true/false/yes/no/1/0.

A section name is the emission line itself, unless the section sets `line`
-- then the section name is just a label, which lets the same line appear
more than once (different fitters, palettes, or parameters) and becomes
that entry's default `run-name` so the two runs' figures do not collide.

Options given on the command line alongside `--config` act as the base for
every entry, so precedence runs: command line < [DEFAULT] < section.

Every entry is attempted even if an earlier one fails; failures are
reported in a summary at the end and make the script exit non-zero.
"""
import argparse
import configparser
import glob
import os
import sys
import traceback
import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning
from irspec.datacube import Datacube
from irspec.spaxel_fit import Spaxelcube
from irspec.plotparams import PlotParams
from irspec.emission_io import read_line_params
from irspec import paths

warnings.filterwarnings("ignore", category=AstropyDeprecationWarning, append=True)

DEFAULT_REDSHIFT = 0.044601  # IR 23128
# Data and output roots resolve lazily through irspec.paths
# (IRSPEC_DATA_ROOT / IRSPEC_OUTPUT_ROOT).
NIRSPEC_CUBE = "g395h-f290lp_s3d.fits"


def cube_path(line_params, source, instrument):
    """Resolves the datacube file for a line, or raises FileNotFoundError."""
    data_dir = str(paths.cube_dir(source))
    if instrument == "NIRSpec":
        pattern = os.path.join(data_dir, NIRSPEC_CUBE)
    else:
        channel, subchannel = line_params[0], line_params[1]
        pattern = os.path.join(data_dir, f"ch{channel}-{subchannel}_s3d.fits")
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No datacube matching {pattern}")
    return matches[0]


def output_dir(line, source, output_root):
    """Builds (and creates if needed) the per-line data-product directory.
    Returned with a trailing separator, as Spaxelcube joins by string
    concatenation."""
    if output_root is None:
        output_root = str(paths.spaxs_dir())
    if source == "north":
        path = os.path.join(output_root, "north", line)
    else:
        path = os.path.join(output_root, line)
    os.makedirs(path, exist_ok=True)
    return os.path.join(path, "")


def plot_dir(line, run_name, plot_root):
    """Builds (and creates if needed) the per-line, per-run plot directory
    under <plot-root>/line_maps/.

    The run name segment is what separates one fitting run's figures from
    another's; the source is not a path component here, so a run name that
    distinguishes it (the default does) is what keeps north and south maps
    from overwriting each other.
    """
    if plot_root is None:
        plot_root = str(paths.plots_dir())
    path = os.path.join(plot_root, "line_maps", line, run_name)
    os.makedirs(path, exist_ok=True)
    return path


# Batch configuration
#
# A config entry is turned back into a command line and re-parsed by the
# same argparse parser that handles a direct invocation. That is what keeps
# the two entry points from drifting: every option's type, `choices`, and
# error message is defined once, and a config file can never accept a
# setting the command line would reject.


def _parser_option_specs(parser):
    """Returns (valid dests, dests that are bare on/off flags) for `parser`.

    `nargs == 0` is what distinguishes a store_true/store_const flag, which
    a config file expresses as a boolean, from an option that carries a
    value.
    """
    valid, flags = set(), set()
    for action in parser._actions:
        if not action.option_strings or action.dest == "help":
            continue
        valid.add(action.dest)
        if action.nargs == 0:
            flags.add(action.dest)
    return valid, flags


def _as_bool(value, key, path):
    """Parses a config boolean, accepting the usual INI spellings."""
    try:
        return configparser.ConfigParser.BOOLEAN_STATES[str(value).strip().lower()]
    except KeyError:
        raise ValueError(
            f"{path}: {key} expects a boolean (true/false, yes/no, 1/0), "
            f"got {value!r}") from None


def _settings_to_argv(line, settings, flag_dests, valid_dests, path, label):
    """Renders one entry's settings as the command line equivalent."""
    argv = [line]
    for key, value in settings.items():
        if key not in valid_dests:
            raise ValueError(
                f"{path}: [{label}] has unknown setting {key!r}. "
                f"Valid settings: {sorted(valid_dests)}")
        option = "--" + key.replace("_", "-")
        if key in flag_dests:
            if _as_bool(value, key, path):
                argv.append(option)
        else:
            argv.append(option)
            argv.append(str(value))
    return argv


def read_run_config(path, parser, base_settings):
    """Reads a batch config file into an ordered list of parsed arg
    namespaces, one per entry, in the order the file lists them.

    Args:
        path (str): The INI file to read.
        parser (argparse.ArgumentParser): The parser each entry is validated
            through.
        base_settings (dict): Settings from the command line that seed every
            entry, i.e. the lowest-precedence layer.

    Returns:
        list: ``(label, line, argparse.Namespace)`` triples.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No config file at {path}")

    # interpolation=None because line names and paths may contain '%', which
    # the default interpolation would treat as a reference and choke on.
    config = configparser.ConfigParser(interpolation=None,
                                       inline_comment_prefixes=("#", ";"))
    # Keys are option names, so normalize the dash/underscore spellings onto
    # argparse's dest. Section names keep their case ('[NeIII]' is a line).
    config.optionxform = lambda option: option.strip().lower().replace("-", "_")
    config.read(path)

    valid_dests, flag_dests = _parser_option_specs(parser)
    # A config cannot name another config, and the lines to run come from
    # the section headers rather than from a setting.
    valid_dests -= {"config"}

    entries = []
    for label in config.sections():
        settings = dict(config[label])
        # A section that names its line explicitly is a labelled run, so the
        # same line can appear repeatedly; the label then names the run, and
        # hence its plot directory, keeping the entries' figures apart.
        line = settings.pop("line", None)
        if line is not None:
            settings.setdefault("run_name", label)
        else:
            line = label
        merged = {**base_settings, **settings}
        argv = _settings_to_argv(line, merged, flag_dests, valid_dests,
                                 path, label)
        try:
            entries.append((label, line, parser.parse_args(argv)))
        except SystemExit:
            # argparse has already printed what was wrong; name the section
            # it came from, which the reconstructed command line does not say.
            print(f"{path}: invalid settings in section [{label}]",
                  file=sys.stderr)
            raise

    if not entries:
        raise ValueError(f"{path}: no line sections found")
    return entries


def render_maps(spaxel_fit, fig_dir):
    """Renders the standard set of spaxel maps from a completed fit."""
    spaxel_fit.render_multicomponent_plot(fig_dir, savefig=True)
    spaxel_fit.render_totflux_plot(fig_dir, savefig=True)
    spaxel_fit.render_rel_vel_plot("G1CEN", fig_dir, savefig=True)
    spaxel_fit.render_rel_vel_plot("G2CEN", fig_dir, savefig=True)
    spaxel_fit.render_sigma_plot("G1SIGMA", fig_dir, savefig=True)
    spaxel_fit.render_sigma_plot("G2SIGMA", fig_dir, savefig=True)
    spaxel_fit.render_moment_panels(fig_dir, savefig=True)
    spaxel_fit.render_chisqr_plot(fig_dir, savefig=True)


def process_line(line, args, line_dict):
    """Fits (or reloads) one emission line's spaxel map and renders it."""
    print(f"=== {line} ({args.instrument}, {args.source}, {args.fitter}, "
          f"{args.palette}) ===")
    path = cube_path(line_dict[line], args.source, args.instrument)
    print(f"datacube: {path}")
    datacube = Datacube(path, redshift=args.redshift)
    out_dir = output_dir(line, args.source, args.output_root)
    fig_dir = plot_dir(line, args.run_name, args.plot_root)
    spaxel_fit = Spaxelcube(datacube, line, out_dir, plot_output=fig_dir,
                            palette=args.palette, scaling=args.scaling,
                            contour_image=args.contour_image)

    if args.fitter == "q3dfit":
        # q3dfit writes its own per-spaxel and collected products (not
        # spaxel_fit.fitparams), so --workers maps to q3dfit's MPI core
        # count and the maps come from the render_q3d_* methods.
        spaxel_fit.q3dfit_fit(ncores=args.workers or 1)
        if not args.no_render:
            spaxel_fit.render_q3d_maps(fig_dir, savefig=True)
            print(f"maps written to {fig_dir}")
        print(f"q3dfit products written to {out_dir}")
        return

    if args.skip_fit:
        spaxel_fit.load_fit(out_dir + "twogaussian_raw.dat")
    else:
        spaxel_fit.two_gaussian_fit(n_workers=args.workers)
    if not args.no_render:
        render_maps(spaxel_fit, fig_dir)
        print(f"maps written to {fig_dir}")
    print(f"outputs written to {out_dir}")


def build_parser():
    """Builds the argument parser."""
    parser = argparse.ArgumentParser(
        description="Fit spaxel maps for one or more emission lines.")
    parser.add_argument("lines", nargs="*",
                        help="Emission line name(s) as listed in fitparams.dat, "
                             "e.g. '[NeIII]' '[H_2_S_1]'. Omit when using "
                             "--config")
    parser.add_argument("--source", choices=["south", "north"], default="south",
                        help="Which nucleus's datacubes to fit (default: south)")
    parser.add_argument("--instrument", choices=["MRS", "NIRSpec"], default="MRS",
                        help="Instrument whose cube covers the line (default: MRS)")
    parser.add_argument("--redshift", type=float, default=DEFAULT_REDSHIFT,
                        help=f"Source redshift (default: {DEFAULT_REDSHIFT})")
    parser.add_argument("--fitter", choices=["twogauss", "q3dfit"],
                        default="twogauss",
                        help="Fitting backend: the built-in two-Gaussian "
                             "fitter (default) or q3dfit")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallelism for fitting. twogauss: worker "
                             "processes (default all CPUs; 1 disables "
                             "multiprocessing). q3dfit: MPI core count "
                             "(needs mpiexec; default 1)")
    parser.add_argument("--output-root", default=None,
                        help="Root directory for per-line data products "
                             "(default: spaxs/ under the irspec output root)")
    parser.add_argument("--plot-root", default=None,
                        help="Root directory for rendered figures; maps go to "
                             "<plot-root>/line_maps/<line>/<run-name>/ "
                             "(default: plots/ under the irspec output root)")
    parser.add_argument("--run-name", default=None,
                        help="Name of this fitting run, used as the directory "
                             "the line's maps are written to. Pick a name that "
                             "identifies the run's configuration so repeat runs "
                             "don't overwrite each other "
                             "(default: <source>_<fitter>)")
    parser.add_argument("--palette", choices=list(PlotParams.PALATTES),
                        default=PlotParams.DEFAULT_PALATTE,
                        help="Figure palette. Saved maps are tagged with it "
                             "(..._light.png / ..._dark.png), so rendering both "
                             "into one run directory keeps both "
                             f"(default: {PlotParams.DEFAULT_PALATTE})")
    parser.add_argument("--scaling", choices=list(PlotParams.SCALINGS),
                        default="presentation",
                        help="Font/DPI scaling for the rendered figures "
                             "(default: presentation)")
    parser.add_argument("--contour-image", default=None, metavar="IMAGE",
                        help="Overlay the maps with an image's intensity "
                             "contours instead of the cube's footprint "
                             "outline. Takes a filter name ('f200w', 'f356w', "
                             "'f560w', 'f770w', 'f1500w') or a path to an "
                             "imaging FITS file. The tag is added to the "
                             "saved filenames, so the overlaid maps sit "
                             "beside the footprint ones rather than "
                             "replacing them (default: footprint outline)")
    parser.add_argument("--skip-fit", action="store_true",
                        help="Skip fitting and re-render maps from the saved "
                             "fit table")
    parser.add_argument("--no-render", action="store_true",
                        help="Fit only; skip rendering the maps")
    parser.add_argument("--config", metavar="FILE", default=None,
                        help="INI file listing lines to run in order, each "
                             "with its own settings (see the module docstring). "
                             "Command-line options given alongside it become "
                             "the base every entry inherits")
    return parser


def finalize_args(args, parser, line_dict):
    """Applies the defaults and validation shared by both entry points.

    Runs for every config entry as well as for a direct invocation, so a
    batch cannot slip past a check the command line enforces.
    """
    if args.fitter == "q3dfit" and args.skip_fit:
        parser.error("--skip-fit applies only to the twogauss fitter")

    if args.run_name is None:
        args.run_name = f"{args.source}_{args.fitter}"
    # The run name is a single directory under the line, so reject anything
    # that would escape it (a stray "../" or a nested path).
    if args.run_name != os.path.basename(args.run_name.rstrip("/\\")) \
            or args.run_name in (".", ".."):
        parser.error(f"--run-name must be a single directory name, "
                     f"got {args.run_name!r}")

    unknown = [line for line in args.lines if line not in line_dict]
    if unknown:
        parser.error(f"Unknown line(s) {unknown}. "
                     f"Available: {sorted(line_dict)}")
    return args


def run_entries(entries, line_dict):
    """Works through (label, line, args) entries in order, running every one
    even if an earlier entry fails.

    A batch is long-running and its entries are independent, so one line
    without a datacube should not discard the work queued behind it. The
    traceback is still printed as it happens, and the failures are collected
    for the closing summary.

    Returns:
        list: the labels that failed.
    """
    failed = []
    for label, line, args in entries:
        try:
            process_line(line, args, line_dict)
        except Exception:
            traceback.print_exc()
            print(f"!!! {label} failed; continuing")
            failed.append(label)
    return failed


def main():
    parser = build_parser()
    args = parser.parse_args()
    line_dict = read_line_params()

    if args.config:
        if args.lines:
            parser.error("pass lines either as arguments or in --config, "
                         "not both")
        # Recover just the options the user actually typed, by diffing the
        # parsed values against the parser's own defaults. (Asking argparse
        # via argument_default=SUPPRESS does not work: an action's explicit
        # `default=` takes precedence, so those options come back anyway --
        # and `--run-name`'s None would reach a config entry as the literal
        # string "None".) These seed every entry, under the file's settings.
        defaults = vars(parser.parse_args([]))
        base = {key: value for key, value in vars(args).items()
                if key not in ("lines", "config")
                and value != defaults.get(key)}
        try:
            entries = read_run_config(args.config, parser, base)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        entries = [(label, line, finalize_args(entry_args, parser, line_dict))
                   for label, line, entry_args in entries]
        print(f"{len(entries)} entries from {args.config}")
    else:
        if not args.lines:
            parser.error("no lines given; pass line names or --config FILE")
        finalize_args(args, parser, line_dict)
        entries = [(line, line, args) for line in args.lines]

    failed = run_entries(entries, line_dict)
    if failed:
        print(f"\n{len(failed)}/{len(entries)} entries failed: "
              f"{', '.join(failed)}")
        sys.exit(1)
    print(f"\nall {len(entries)} entries completed")


if __name__ == "__main__":
    main()
