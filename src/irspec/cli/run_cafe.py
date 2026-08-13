"""Run a CAFE fit through the :class:`~irspec.cubespec.CubeSpec` interface.

Fits a CRETA extraction with CAFE and writes the result into a matching,
cleanly-named directory. Two extraction types are supported:

``single``
    Fits the one-aperture spectrum with :meth:`CubeSpec.perform_fit`
    (``specmod``), producing ``{NAME}_cafefit.asdf`` and friends.
``grid``
    Fits every spaxel of a grid extraction individually with
    :meth:`CubeSpec.perform_grid_fit` (``cubemod``), producing a per-spaxel
    parameter cube ``{NAME}_parcube.fits`` and ``{NAME}_fitpars.ini``.

Outputs are organised by type (see ``runner_utils``)::

    cafe_output/
      single/
        runs.txt                <- manifest of every single fit
        {NAME}_single/          <- this fit's files, rebased on NAME
      grid/
        runs.txt                <- manifest of every grid fit
        {NAME}_grid/            <- this fit's parameter cube, rebased on NAME

Plots mirror that layout under ``plots/cafe`` so a run's figures sit at the
same relative path as its numbers::

    plots/cafe/
      single/
        {NAME}_single/
          {NAME}_fit.png           <- the fitted spectrum + residuals
          {NAME}_fit_nirspec.png   <- the same fit and residuals, one figure
          {NAME}_fit_ch1.png          per channel, over that channel's
          ...                         wavelengths only
          {NAME}_fit_ch4.png

A single fit is rendered automatically once the fit finishes (``--no-plot``
suppresses everything, ``--no-channel-plot`` just the per-channel figures).
Because the figures are rebuilt from the run's saved
``{NAME}_cafefit.asdf``, an existing fit can be re-rendered as often as needed
without refitting::

    pixi run python run_cafe.py --param-file ir23128s_single1S_param.txt \
        --name single1S_r1.0as_SB --mode SB --plot-only

Use the same ``--name`` (or the same ``--param-file``, which yields the same
default name) as the CRETA run so the fit reads from
``creta_extractions/{single,grid}/{NAME}_{single,grid}/`` and writes to the
matching CAFE directory. The extraction type is taken from the parameter-file
name ("single"/"grid"), or forced with ``--extraction``.

CAFE derives its file names from the input spectrum/cube file name. For single
extractions ``CubeSpec.perform_fit`` reads a fixed
``<target>_SingleExt_r<r_ap>as.fits``; since ``run_creta.py`` renamed that to
``{NAME}.fits``, this script briefly bridges the expected name (via a symlink).
For grid extractions ``CubeSpec.perform_grid_fit`` locates the ``{NAME}_cube.fits``
directly (no bridge needed). In both cases CAFE's output directory and files are
renamed to the ``{NAME}`` scheme afterwards.

A grid fit runs for hours, so it is checkpointed as it goes: every
``--checkpoint-every`` spaxels (10 by default) the parameter cube so far is
written to a sidecar beside the eventual output. If the fit is interrupted --
Ctrl-C, ``kill``, a scheduler reclaiming the node -- it flushes that sidecar on
the way out and prints the command to continue::

    pixi run python run_cafe.py --param-file ir23128s_grid1_param.txt \
        --name southern_grid --mode SB --resume

The resumed fit skips every spaxel already fitted (and every one already found
unfittable), while still using the fitted ones to seed their neighbours, so it
takes the same path through the grid an uninterrupted run would have. Only the
spaxel in progress when the fit stopped is redone. The checkpoint is deleted
once the fit completes, and one belonging to a different cube, mode or redshift
is refused rather than mixed into the run.

Run this from ``irspec/src`` so the default relative paths resolve, e.g.::

    pixi run python run_cafe.py --param-file ir23128s_single1S_param.txt \
        --name southern_nucleus --mode SB

    pixi run python run_cafe.py --param-file ir23128s_grid1_param.txt \
        --name southern_grid --mode SB
"""

import argparse
import glob
import os
import shutil
import sys

from irspec.cli import runner_utils as ru
from irspec import paths
from irspec.cubespec import CubeSpec


def build_cubespec(args, name, creta_output, cafe_output, plot_output=None):
    """Construct a :class:`CubeSpec`. Both output directories must exist before
    construction (CubeSpec verifies them); the plot directory is created lazily
    by CubeSpec when something is actually rendered."""
    os.makedirs(creta_output, exist_ok=True)
    os.makedirs(cafe_output, exist_ok=True)

    datacube = None
    if args.datacube is not None:
        from irspec.datacube import Datacube
        datacube = Datacube(args.datacube, redshift=args.redshift,
                            verbose=args.verbose)

    return CubeSpec(
        input_path=args.input_path,
        param_path=args.param_path,
        param_file=args.param_file,
        creta_output_path=creta_output,
        cafe_output_path=cafe_output,
        redshift=args.redshift,
        datacube=datacube,
        mode=args.mode,
        name=name,
        plot_output_path=plot_output,
        verbose=args.verbose,
    )


def _rebase_prefix(directory, old_prefix, new_prefix):
    """Rename every file in ``directory`` whose name starts with ``old_prefix``
    so that prefix becomes ``new_prefix`` (e.g. ``grid1_cube_parcube.fits`` ->
    ``grid1_parcube.fits``). Returns the sorted list of resulting file names.
    """
    result = []
    for path in sorted(glob.glob(os.path.join(directory, "*"))):
        if os.path.isdir(path):
            continue
        old = os.path.basename(path)
        new = new_prefix + old[len(old_prefix):] if old.startswith(old_prefix) else old
        if new != old:
            os.replace(path, os.path.join(directory, new))
        result.append(new)
    return sorted(result)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a CAFE single or grid fit on a CRETA extraction via the CubeSpec interface.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--param-file", required=True,
                        help="Parameter file name (e.g. ir23128s_single1S_param.txt).")
    parser.add_argument("--name", default=None,
                        help="Run name; must match the CRETA run. Drives the "
                             "{NAME}_{single|grid} output directory and file "
                             "names. Defaults to the id parsed from the param file.")
    parser.add_argument("--creta-name", default=None,
                        help="Name of the CRETA run to read, when it differs "
                             "from --name (the fit's name often carries a mode "
                             "tag the extraction does not). Defaults to --name.")
    parser.add_argument("--extraction", default="auto",
                        choices=["auto", "single", "grid"],
                        help="Extraction type; 'auto' reads it from the param file name.")
    parser.add_argument("--input-path", default=None,
                        help="Directory containing the input datacubes "
                             "(default: cubes/ under the irspec data root).")
    parser.add_argument("--param-path", default=None,
                        help="Directory containing the parameter file "
                             "(default: param_files/ under the irspec data root).")
    parser.add_argument("--creta-root", default=None,
                        help="Root holding the CRETA single/ or grid/ group with the extraction "
                             "(default: creta_extractions/ under the irspec output root).")
    parser.add_argument("--cafe-root", default=None,
                        help="Root under which the single/ or grid/ group and fit are created "
                             "(default: cafe_output/ under the irspec output root).")
    parser.add_argument("--plots-root", default=None,
                        help="Root under which this run's plot directory is created, "
                             "mirroring the cafe_output layout "
                             "(default: plots/cafe/ under the irspec output root).")
    parser.add_argument("--redshift", type=float, default=ru.DEFAULT_REDSHIFT,
                        help="Systemic redshift of the target.")
    parser.add_argument("--mode", default="AGN",
                        choices=["AGN", "SB", "SB_and_AGN"],
                        help="CAFE spectral parameter mode "
                             "(selects inpars_jwst_nirspec-miri_<mode>.ini).")
    parser.add_argument("--extract", default="Flux_st",
                        help="Grid only: which CRETA flux column to fit "
                             "('Flux_st' stitched, or 'Flux_ap').")
    parser.add_argument("--force-all-lines", action="store_true",
                        help="Grid only: fit every catalog line at every spaxel.")
    parser.add_argument("--resume", action="store_true",
                        help="Grid only: continue an interrupted fit from its "
                             "checkpoint instead of refitting the whole grid. "
                             "Spaxels already done are kept and still seed "
                             "their neighbours; only the spaxel that was in "
                             "progress when the fit stopped is refitted.")
    parser.add_argument("--checkpoint-every", type=int, default=10,
                        metavar="N",
                        help="Grid only: checkpoint the parameter cube every N "
                             "spaxels, bounding what an interrupt can cost. "
                             "0 disables checkpointing.")
    parser.add_argument("--datacube", default=None,
                        help="Optional path to a datacube to attach (unused by the fit).")
    parser.add_argument("--plot-only", action="store_true",
                        help="Skip fitting and (re)render the plots of an existing "
                             "run from its saved CAFE products.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Fit without rendering any plots.")
    parser.add_argument("--no-channel-plot", action="store_true",
                        help="Skip the per-channel (NIRSpec, CH1-4) breakdown "
                             "of the fit.")
    parser.add_argument("--plot-dpi", type=int, default=250,
                        help="Resolution of the rendered figures.")
    parser.add_argument("--residual", default="percent",
                        choices=["percent", "sigma"],
                        help="Residual panel units. 'percent' of the flux, or "
                             "'sigma' to divide by the measurement error -- "
                             "worth using only once the extraction carries "
                             "trustworthy uncertainties.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


def plot_single(args, spec, name):
    """Render the plots of a single-aperture run into the run's plot directory
    (``spec.plot_output_path``, set in :func:`build_cubespec`): the full-range
    fit, plus a per-channel breakdown of the same fit.

    Reads the fit back from ``{NAME}_cafefit.asdf``, so this works both right
    after a fit and on its own (``--plot-only``) against an earlier one.

    Returns:
        list: The paths written.
    """
    title = f'{name} ({spec.target}, r = {spec.asec}")'
    paths = [spec.plot_fit(save=True, dpi=args.plot_dpi, title=title,
                           residual=args.residual)]
    if not args.no_channel_plot:
        paths.extend(spec.plot_channels(save=True, dpi=args.plot_dpi,
                                        title=title, residual=args.residual))
    for path in paths:
        print(f"[run_cafe] plot:     {path}")
    return paths


def run_single(args, spec, name, cafe_group, creta_name):
    """Fit a single-aperture extraction and organise its CAFE outputs."""
    # perform_fit reads a fixed "<target>_SingleExt_r<r_ap>as.fits". run_creta
    # renamed that extraction to "<creta_name>.fits", so bridge the expected
    # name. The extraction is named for the CRETA run, not this fit: one
    # extraction feeds several fits that differ only by --mode.
    expected_fn = f"{spec.target}_SingleExt_r{spec.asec}as.fits"
    expected_path = os.path.join(spec.creta_output_path, expected_fn)
    clean_source = os.path.join(spec.creta_output_path, f"{creta_name}.fits")
    bridged = False
    if not os.path.exists(expected_path):
        if os.path.exists(clean_source):
            os.symlink(os.path.basename(clean_source), expected_path)
            bridged = True
        else:
            raise FileNotFoundError(
                f"No CRETA extraction found in {spec.creta_output_path}:\n"
                f"  expected {os.path.basename(clean_source)} (or {expected_fn}).\n"
                f"Run run_creta.py for this run first."
            )

    # CAFE writes into a subdir named after the spectrum stem; clean any stale
    # copy and the destination run dir before (re)fitting.
    result_stem = "".join(expected_fn.split(".")[:-1])  # <target>_SingleExt_r<..>as
    intermediate = os.path.join(cafe_group, result_stem)
    cafe_run_dir = os.path.join(cafe_group, ru.run_dirname(name, "single"))
    shutil.rmtree(intermediate, ignore_errors=True)
    shutil.rmtree(cafe_run_dir, ignore_errors=True)

    print(f"[run_cafe] name={name} target={spec.target} mode={args.mode} "
          f"type=single z={spec.redshift}")
    print(f"[run_cafe] spectrum: {clean_source if bridged else expected_path}")
    print(f"[run_cafe] output:   {cafe_run_dir}")

    try:
        spec.perform_fit()
    finally:
        if bridged:
            os.remove(expected_path)  # remove the symlink, not its target

    # Move CAFE's output into the clean {NAME}_single dir and rebase file names.
    os.replace(intermediate, cafe_run_dir)
    files = ru.rename_outputs(cafe_run_dir, name, target=spec.target)

    plots = [] if args.no_plot else plot_single(args, spec, name)

    props = {
        "name": name,
        "type": "single",
        "param_file": args.param_file,
        "target": spec.target,
        "redshift": spec.redshift,
        "mode": args.mode,
        "inpars": f"inpars_jwst_nirspec-miri_{args.mode}.ini",
        "source": os.path.relpath(clean_source, args.creta_root),
        "output_dir": os.path.relpath(cafe_run_dir, args.cafe_root),
        "files": ", ".join(files),
    }
    if plots:
        props["plots"] = ", ".join(os.path.relpath(p, args.plots_root)
                                   for p in plots)
    manifest = ru.write_manifest(args.cafe_root, name, "single", props)
    print(f"[run_cafe] wrote {len(files)} files; manifest -> {manifest}")


def run_grid(args, spec, name, cafe_group):
    """Fit every spaxel of a grid extraction and organise its CAFE outputs."""
    # The CRETA grid cube is "<name>_cube.fits" in the run dir; perform_grid_fit
    # locates it via CubeSpec.run_name, so no filename bridge is needed.
    cube_path = spec._grid_cube_path()
    # CAFE names its outputs after the cube stem (dots -> 'p'); for a dot-free
    # "<name>_cube.fits" that is simply "<name>_cube".
    result_stem = "p".join(os.path.basename(cube_path).split(".")[:-1])
    intermediate = os.path.join(cafe_group, result_stem)
    cafe_run_dir = os.path.join(cafe_group, ru.run_dirname(name, "grid"))
    # The intermediate directory is where CAFE keeps its checkpoint, so a
    # resume must not start by deleting it -- that is precisely the state being
    # resumed from. A fresh run still clears it, so a stale checkpoint from an
    # abandoned attempt can never bleed into one.
    if not args.resume:
        shutil.rmtree(intermediate, ignore_errors=True)
    shutil.rmtree(cafe_run_dir, ignore_errors=True)

    print(f"[run_cafe] name={name} target={spec.target} mode={args.mode} "
          f"type=grid z={spec.redshift}")
    print(f"[run_cafe] cube:   {cube_path}")
    print(f"[run_cafe] output: {cafe_run_dir}")
    if args.resume:
        print(f"[run_cafe] resume: {intermediate}")

    spec.perform_grid_fit(extract=args.extract,
                          force_all_lines=args.force_all_lines,
                          resume=args.resume,
                          checkpoint_every=args.checkpoint_every)

    # Move CAFE's output into the clean {NAME}_grid dir and drop the "_cube"
    # infix CAFE inherited from the cube file name.
    os.replace(intermediate, cafe_run_dir)
    files = _rebase_prefix(cafe_run_dir, result_stem, name)

    props = {
        "name": name,
        "type": "grid",
        "param_file": args.param_file,
        "target": spec.target,
        "redshift": spec.redshift,
        "mode": args.mode,
        "inpars": f"inpars_jwst_nirspec-miri_{args.mode}.ini",
        "extract": args.extract,
        "grid": (f"{spec.param_dict.get('nx_steps')}x"
                 f"{spec.param_dict.get('ny_steps')} "
                 f"spax={spec.param_dict.get('spax_size')}as"),
        "cube": os.path.relpath(cube_path, args.creta_root),
        "output_dir": os.path.relpath(cafe_run_dir, args.cafe_root),
        "files": ", ".join(files),
    }
    manifest = ru.write_manifest(args.cafe_root, name, "grid", props)
    print(f"[run_cafe] wrote {len(files)} files; manifest -> {manifest}")


def main():
    args = parse_args()
    if args.input_path is None:
        args.input_path = str(paths.cube_dir()) + "/"
    if args.param_path is None:
        args.param_path = str(paths.param_dir()) + "/"
    if args.creta_root is None:
        args.creta_root = str(paths.creta_dir()) + "/"
    if args.cafe_root is None:
        args.cafe_root = str(paths.cafe_dir()) + "/"
    if args.plots_root is None:
        args.plots_root = str(paths.plots_dir() / "cafe") + "/"
    name = args.name or ru.extraction_id(args.param_file)
    run_type = ru.detect_type(args.param_file, args.extraction)
    if run_type is None:
        raise ValueError(
            "Could not determine extraction type: the parameter file name "
            "contains neither 'single' nor 'grid'. Pass --extraction explicitly."
        )

    creta_name = args.creta_name or name
    creta_run_dir = ru.run_output_dir(args.creta_root, creta_name, run_type)
    cafe_group = ru.group_dir(args.cafe_root, run_type)
    plot_run_dir = ru.run_output_dir(args.plots_root, name, run_type)
    if not os.path.isdir(creta_run_dir):
        if args.plot_only:
            # Replotting reads the CAFE products only, so a missing extraction
            # dir (e.g. a mode-tagged fit name) must not block the render.
            print(f"[run_cafe] note: no CRETA run dir at {creta_run_dir}; "
                  f"replotting reads the CAFE products only.")
            creta_run_dir = ru.group_dir(args.creta_root, run_type)
        else:
            raise FileNotFoundError(
                f"CRETA run directory not found:\n  {creta_run_dir}\n"
                f"Run run_creta.py with the same --name/--param-file first "
                f"(or point at it with --creta-name)."
            )

    spec = build_cubespec(args, name, creta_run_dir, cafe_group, plot_run_dir)

    # Guard against a run-type/param mismatch (e.g. --extraction grid on a
    # "single" param file), which CubeSpec would otherwise reject deeper in.
    detected = spec.param_dict.get("type")
    if detected is not None and detected != run_type:
        raise ValueError(
            f"Requested fit '{run_type}' but the param file looks like "
            f"'{detected}'. Fix --extraction or the param file."
        )

    if args.plot_only:
        if run_type != "single":
            raise NotImplementedError(
                "--plot-only currently renders single-aperture fits only; "
                "grid runs are plotted from their parameter cube elsewhere."
            )
        cafe_run_dir = ru.run_output_dir(args.cafe_root, name, run_type)
        if not os.path.isdir(cafe_run_dir):
            raise FileNotFoundError(
                f"CAFE run directory not found:\n  {cafe_run_dir}\n"
                f"Fit this run first (drop --plot-only)."
            )
        print(f"[run_cafe] name={name} type=single (plot only)")
        print(f"[run_cafe] fit:      {cafe_run_dir}")
        plot_single(args, spec, name)
    elif run_type == "single":
        run_single(args, spec, name, cafe_group, creta_name)
    else:
        run_grid(args, spec, name, cafe_group)
        if not args.no_plot:
            print("[run_cafe] note: grid runs are not plotted by this script yet.")

    print("[run_cafe] Done.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        # A grid fit runs for hours and stopping it early is a normal thing to
        # do, not a crash. CAFE has already flushed its checkpoint by the time
        # this lands, so all that is left is to say how to pick it back up.
        argv = [a for a in sys.argv[1:] if a != "--resume"] + ["--resume"]
        print("\n[run_cafe] interrupted. Resume this fit with:\n"
              f"  python {os.path.basename(sys.argv[0])} {' '.join(argv)}",
              file=sys.stderr)
        sys.exit(130)
