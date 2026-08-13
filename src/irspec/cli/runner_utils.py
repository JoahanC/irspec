"""Shared helpers for the ``run_creta.py`` / ``run_cafe.py`` runner scripts.

Centralises the output-organisation policy so it is defined in exactly one
place:

* CRETA and CAFE outputs are grouped by extraction type under ``single/`` and
  ``grid/`` subdirectories of their respective roots.
* Each run gets its own directory named ``{NAME}_{single|grid}`` inside the
  matching group directory.
* Output file names have the CRETA/CAFE machine tokens (``_SingleExt_r1.0as``,
  ``_GridExt_35x35_s0.39as`` and the like — which bake in the aperture/grid
  size) stripped and are rebased on ``NAME``.
* A ``runs.txt`` manifest in each group directory records every run and its
  defining properties.

``NAME`` defaults to the identifier derived from the parameter-file name, so
passing the same ``--param-file`` to both scripts yields a consistent name
without repeating yourself; ``--name`` overrides it.
"""

import glob
import os
import re
from datetime import datetime


# Default systemic redshift of IR23128-S, the object these param files target.
DEFAULT_REDSHIFT = 0.044601

# Manifest file written into each group directory (single/ and grid/).
MANIFEST_NAME = "runs.txt"

# Matches the CRETA/CAFE token that embeds the aperture radius or grid size,
# e.g. "IR23128-S_SingleExt_r1.0as" or "IR23128-S_GridExt_35x35_s0.39as".
_SIZE_TOKEN = re.compile(
    r"^(?P<base>.+?)_(?:SingleExt_r[\d.]+as|GridExt_\d+x\d+_s[\d.]+as)"
)


def extraction_id(param_file):
    """Derive a short run identifier (e.g. ``single1S``) from the param file
    name, mirroring ``CubeSpec.extraction_id``. Falls back to the file stem."""
    base = os.path.basename(param_file)
    parts = base.split("_")
    return parts[1] if len(parts) > 1 else os.path.splitext(base)[0]


def detect_type(param_file, override="auto"):
    """Return 'single' or 'grid' for a run, from ``override`` or, when it is
    'auto', from the parameter-file name. Returns ``None`` if undetermined."""
    if override != "auto":
        return override
    base = os.path.basename(param_file).lower()
    if "single" in base:
        return "single"
    if "grid" in base:
        return "grid"
    return None


def run_dirname(name, run_type):
    """The per-run directory name, e.g. ``mynucleus_single``."""
    return f"{name}_{run_type}"


def group_dir(root, run_type):
    """The ``single/`` or ``grid/`` group directory under ``root``."""
    return os.path.join(root, run_type)


def run_output_dir(root, name, run_type):
    """Full per-run output directory: ``{root}/{run_type}/{name}_{run_type}``."""
    return os.path.join(group_dir(root, run_type), run_dirname(name, run_type))


def clean_filename(filename, name, target=None):
    """Return ``filename`` with any CRETA/CAFE size token stripped and the base
    rebased on ``name``. Files that carry no size token but start with the
    target prefix are rebased too; anything else is returned unchanged."""
    match = _SIZE_TOKEN.match(filename)
    if match:
        return name + filename[match.end():]
    if target and filename.startswith(target + "_"):
        return name + "_" + filename[len(target) + 1:]
    return filename


def rename_outputs(directory, name, target=None):
    """Rename every file in ``directory`` to the cleaned scheme in place.

    Returns the sorted list of resulting file names.
    """
    result = []
    for path in sorted(glob.glob(os.path.join(directory, "*"))):
        if os.path.isdir(path):
            continue
        old = os.path.basename(path)
        new = clean_filename(old, name, target=target)
        if new != old:
            os.replace(path, os.path.join(directory, new))
        result.append(new)
    return sorted(result)


def _strip_block(text, header):
    """Remove an existing ``[header]`` block (through its trailing separator)
    so a re-run replaces rather than duplicates its manifest entry."""
    if header not in text:
        return text
    out, skipping = [], False
    for line in text.splitlines(keepends=True):
        if line.strip() == header:
            skipping = True
            continue
        if skipping:
            # A separator line (all dashes) ends the block.
            if line.strip() and set(line.strip()) == {"-"}:
                skipping = False
            continue
        out.append(line)
    return "".join(out)


def write_manifest(root, name, run_type, props):
    """Append (or replace) this run's entry in ``{root}/{run_type}/runs.txt``.

    ``props`` is an ordered mapping of property name -> value. Returns the
    manifest path.
    """
    gdir = group_dir(root, run_type)
    os.makedirs(gdir, exist_ok=True)
    manifest = os.path.join(gdir, MANIFEST_NAME)

    header = f"[{run_dirname(name, run_type)}]"
    lines = [header, f"{'timestamp':<16}= {datetime.now().isoformat(timespec='seconds')}"]
    for key, value in props.items():
        lines.append(f"{key:<16}= {value}")
    block = "\n".join(lines) + "\n" + "-" * 64 + "\n"

    existing = ""
    if os.path.exists(manifest):
        with open(manifest) as fh:
            existing = fh.read()
    existing = _strip_block(existing, header)
    if existing and not existing.endswith("\n"):
        existing += "\n"

    with open(manifest, "w") as fh:
        fh.write(existing + block)
    return manifest
