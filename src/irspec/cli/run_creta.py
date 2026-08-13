"""Run a CRETA extraction through the :class:`~irspec.cubespec.CubeSpec` interface.

Drives CRETA's single- or grid-aperture extraction on a set of MIRI/MRS
datacubes and (optionally) reformats the resulting spectrum CSV for Thomas Lai's
web viewer.

Outputs are organised for easy navigation (see ``runner_utils``):

    creta_extractions/
      single/
        runs.txt                <- manifest of every single run
        {NAME}_single/          <- this run's files, rebased on NAME
      grid/
        runs.txt
        {NAME}_grid/

The extraction type is taken from the parameter-file name ("single"/"grid"),
or forced with ``--extraction``. ``NAME`` defaults to the identifier derived
from the parameter file (e.g. ``single1S``); override it once with ``--name``
and pass the same value to ``run_cafe.py`` to keep a run's CRETA and CAFE
outputs together.

Run this from ``irspec/src`` so the default relative paths resolve, e.g.::

    pixi run python run_creta.py --param-file ir23128s_single1S_param.txt \
        --name southern_nucleus --mode SB --rewrite-csv

    pixi run python run_creta.py --param-file ir23128s_grid4_param.txt \
        --name pah_grid

    # also emit a resolution- + flux-matched stitched cube (hybrid by default):
    pixi run python run_creta.py --param-file ir23128s_single1S_param.txt \
        --name southern_nucleus --stitch --stitch-method hybrid
"""

import argparse
import glob
import os
import shutil

from irspec.cli import runner_utils as ru
from irspec import paths
from irspec.cubespec import CubeSpec


def build_stitched_cube(spec, name, run_dir, args):
    """Resolution-match and flux-match the MIRI/MRS sub-band cubes into a single
    stitched cube with :class:`~irspec.cube_stitch.CubeStitcher`, and write it
    into the run directory as ``{name}_stitched.fits``.

    NIRSpec ``g395h`` is excluded on purpose: the MIRI PSF model that drives the
    resolution match does not describe it.

    Args:
        spec (CubeSpec): The constructed run (used for ``input_path``).
        name (str): Run name; drives the output file name.
        run_dir (str): Directory to write the stitched cube into.
        args: Parsed CLI args (``stitch_method``, ``stitch_per_spaxel``,
            ``stitch_reference``, ``verbose``).

    Returns:
        tuple: ``(filename, props)`` -- the written file's base name and a dict
        of stitch metadata for the run manifest.
    """
    from irspec.datacube import Datacube
    from irspec.cube_stitch import CubeStitcher

    files = sorted(glob.glob(os.path.join(spec.input_path, "ch[1-4]-*_s3d.fits")))
    if not files:
        raise FileNotFoundError(
            "No MIRI sub-band cubes (ch[1-4]-*_s3d.fits) found in "
            f"{spec.input_path}; cannot stitch."
        )
    cubes = [Datacube(f, redshift=0, verbose=args.verbose) for f in files]

    # Resolve --stitch-err-rescale: 'none'->None, 'auto'->'auto', else a float.
    er = args.stitch_err_rescale
    if er.lower() in ("none", "off"):
        err_rescale = None
    elif er.lower() == "auto":
        err_rescale = "auto"
    else:
        err_rescale = float(er)

    stitched = CubeStitcher(
        cubes,
        flux_method=args.stitch_method,
        per_spaxel=args.stitch_per_spaxel,
        reference=args.stitch_reference,
        err_rescale=err_rescale,
    ).build()

    out_name = f"{name}_stitched.fits"
    stitched.write(os.path.join(run_dir, out_name))

    bands = list(dict.fromkeys(stitched.bandnames.tolist()))
    anchor = next((e["band"] for e in (stitched.stitch_log or [])
                   if e["overlap"] is None), None)
    # Compact per-band ERR-rescale summary (only the bands actually inflated).
    er_applied = ", ".join(
        f"{e['band']}:x{e['factor']:.2f}"
        for e in (stitched.err_rescale_log or []) if e.get("factor", 1.0) != 1.0)
    props = {
        "stitch_method": args.stitch_method,
        "stitch_per_spaxel": args.stitch_per_spaxel,
        "stitch_reference": args.stitch_reference,
        "stitch_anchor": anchor,
        "stitch_bands": len(bands),
        "stitch_channels": int(stitched.flux.shape[0]),
        "stitch_wave_um": f"{stitched.waves.min():.2f}-{stitched.waves.max():.2f}",
        "stitch_err_rescale": args.stitch_err_rescale,
        "stitch_err_factors": er_applied or "none",
    }
    return out_name, props


def build_cubespec(args, name, creta_output, cafe_output):
    """Construct a :class:`CubeSpec`. Both output directories must exist before
    construction (CubeSpec verifies them), so create them first."""
    os.makedirs(creta_output, exist_ok=True)
    os.makedirs(cafe_output, exist_ok=True)

    # The extraction code path does not use the datacube object, so only build
    # one when a cube is explicitly provided.
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
        verbose=args.verbose,
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a CRETA single/grid extraction via the CubeSpec interface.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--param-file", required=True,
                        help="Parameter file name (e.g. ir23128s_single1S_param.txt).")
    parser.add_argument("--name", default=None,
                        help="Run name; drives the {NAME}_{single|grid} output "
                             "directory and file names. Defaults to the id parsed "
                             "from the param file (e.g. single1S).")
    parser.add_argument("--input-path", default=None,
                        help="Directory containing the input datacubes "
                             "(default: cubes/ under the irspec data root).")
    parser.add_argument("--param-path", default=None,
                        help="Directory containing the parameter file "
                             "(default: param_files/ under the irspec data root).")
    parser.add_argument("--creta-root", default=None,
                        help="Root under which single/ and grid/ groups are created "
                             "(default: creta_extractions/ under the irspec output root).")
    parser.add_argument("--cafe-root", default=None,
                        help="CAFE root (only its group dir is created so the "
                             "CubeSpec constructor is satisfied; default: "
                             "cafe_output/ under the irspec output root).")
    parser.add_argument("--redshift", type=float, default=ru.DEFAULT_REDSHIFT,
                        help="Systemic redshift of the target.")
    parser.add_argument("--mode", default="AGN",
                        choices=["AGN", "SB", "SB_and_AGN"],
                        help="CAFE spectral parameter mode (recorded in the manifest).")
    parser.add_argument("--extraction", default="auto",
                        choices=["auto", "single", "grid"],
                        help="Extraction type; 'auto' reads it from the param file name.")
    parser.add_argument("--rewrite-csv", action="store_true",
                        help="After a single extraction, also write the web-viewer CSV.")
    parser.add_argument("--stitch", action="store_true",
                        help="Also build a resolution- and flux-matched stitched "
                             "cube from the MIRI sub-bands with CubeStitcher "
                             "(NIRSpec g395h excluded).")
    parser.add_argument("--stitch-method", default="hybrid",
                        choices=["chained", "joint", "hybrid"],
                        help="CubeStitcher flux-matching method used by --stitch.")
    parser.add_argument("--stitch-per-spaxel", action="store_true",
                        help="Apply the stitch flux-matching correction per spaxel.")
    parser.add_argument("--stitch-reference", default="snr",
                        help="Anchor band for stitching: 'snr'/'blue'/'red', a "
                             "band name (e.g. ch2_LONG), or an integer index.")
    parser.add_argument("--stitch-err-rescale", default="auto",
                        help="Inflate input ERR before stitching to correct the "
                             "known pipeline ERR under-estimation: 'auto' "
                             "(measure per band, apply only to pre-Build-11.0 "
                             "cubes; no-op on reprocessed data), 'none', or a "
                             "fixed float factor (e.g. 2.5).")
    parser.add_argument("--datacube", default=None,
                        help="Optional path to a datacube to attach (unused by extraction).")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    return parser.parse_args()


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
    name = args.name or ru.extraction_id(args.param_file)

    run_type = ru.detect_type(args.param_file, args.extraction)
    if run_type is None:
        raise ValueError(
            "Could not determine extraction type: the parameter file name "
            "contains neither 'single' nor 'grid'. Pass --extraction explicitly."
        )

    creta_run_dir = ru.run_output_dir(args.creta_root, name, run_type)
    # Start from a clean run directory so stale files from a previous run of the
    # same name are not left behind and then renamed alongside the new outputs.
    shutil.rmtree(creta_run_dir, ignore_errors=True)

    spec = build_cubespec(
        args, name, creta_run_dir, ru.group_dir(args.cafe_root, run_type)
    )

    # Guard against a name/param mismatch (e.g. --extraction grid on a "single"
    # param file), which CubeSpec would otherwise reject deeper in.
    detected = spec.param_dict.get("type")
    if detected is not None and detected != run_type:
        raise ValueError(
            f"Requested extraction '{run_type}' but the param file looks like "
            f"'{detected}'. Fix --extraction or the param file."
        )

    print(f"[run_creta] name={name} target={spec.target} mode={args.mode} "
          f"type={run_type}")
    print(f"[run_creta] input:  {spec.input_path}")
    print(f"[run_creta] output: {creta_run_dir}")

    if run_type == "single":
        spec.perform_single_extraction()
        if args.rewrite_csv:
            spec.rewrite_spec_csv()  # reads the raw CSV before we rename below
    else:
        spec.perform_grid_extraction()

    # Strip aperture/grid-size tokens from the file names and rebase on NAME.
    files = ru.rename_outputs(creta_run_dir, name, target=spec.target)

    # Optionally add a CubeStitcher stitched cube (resolution + flux matched).
    stitch_props = None
    if args.stitch:
        print(f"[run_creta] stitching MIRI sub-bands "
              f"(method={args.stitch_method}, per_spaxel={args.stitch_per_spaxel}, "
              f"reference={args.stitch_reference})...")
        stitch_file, stitch_props = build_stitched_cube(
            spec, name, creta_run_dir, args)
        files.append(stitch_file)
        print(f"[run_creta] wrote stitched cube -> {stitch_file} "
              f"({stitch_props['stitch_channels']} channels, "
              f"{stitch_props['stitch_wave_um']} um, "
              f"anchor {stitch_props['stitch_anchor']})")
        print(f"[run_creta] ERR rescale ({args.stitch_err_rescale}): "
              f"{stitch_props['stitch_err_factors']}")

    # Record the run in the group manifest.
    props = {
        "name": name,
        "type": run_type,
        "param_file": args.param_file,
        "target": spec.target,
        "redshift": spec.redshift,
        "mode": args.mode,
        "output_dir": os.path.relpath(creta_run_dir, args.creta_root),
    }
    if run_type == "single":
        props["aperture_arcsec"] = spec.asec
    else:
        props["grid"] = (f"{spec.param_dict.get('nx_steps')}x"
                         f"{spec.param_dict.get('ny_steps')} "
                         f"step={spec.param_dict.get('step_size')}as")
    if stitch_props is not None:
        props.update(stitch_props)
    props["files"] = ", ".join(files)
    manifest = ru.write_manifest(args.creta_root, name, run_type, props)

    print(f"[run_creta] wrote {len(files)} files; manifest -> {manifest}")
    print("[run_creta] Done.")


if __name__ == "__main__":
    main()
