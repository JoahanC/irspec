"""Filesystem roots for irspec data and outputs.

Science data (cubes, imaging, parameter files) lives outside the package.
Every location is resolved lazily, in this order (first hit wins):

  1. An explicit path passed by the caller.
  2. The IRSPEC_DATA_ROOT / IRSPEC_OUTPUT_ROOT environment variables.
  3. ~/.config/irspec/config.toml with keys ``data_root`` / ``output_root``.
  4. Data: a descriptive error. Outputs: the current working directory.

The pixi workspace sets both environment variables on activation, so inside
`pixi run` nothing needs configuring.
"""
from pathlib import Path
import os
import tomllib

_CONFIG_FILE = Path.home() / ".config" / "irspec" / "config.toml"


def _config_value(key):
    if _CONFIG_FILE.is_file():
        with open(_CONFIG_FILE, "rb") as f:
            return tomllib.load(f).get(key)
    return None


def data_root(explicit=None) -> Path:
    for cand in (explicit, os.environ.get("IRSPEC_DATA_ROOT"),
                 _config_value("data_root")):
        if cand:
            return Path(cand).expanduser()
    raise FileNotFoundError(
        "irspec data root not configured: set IRSPEC_DATA_ROOT, add "
        'data_root = "/path/to/data" to ~/.config/irspec/config.toml, '
        "or pass an explicit path.")


def output_root(explicit=None) -> Path:
    for cand in (explicit, os.environ.get("IRSPEC_OUTPUT_ROOT"),
                 _config_value("output_root")):
        if cand:
            return Path(cand).expanduser()
    return Path.cwd()


def image_data_dir(explicit=None) -> Path:
    return data_root(explicit) / "image_data"


def cube_dir(region="south", explicit=None) -> Path:
    return data_root(explicit) / ("cubes" if region == "south" else "cubes_north")


def param_dir(explicit=None) -> Path:
    return data_root(explicit) / "param_files"


def spaxs_dir(explicit=None) -> Path:
    return output_root(explicit) / "spaxs"


def plots_dir(explicit=None) -> Path:
    return output_root(explicit) / "plots"


def creta_dir(explicit=None) -> Path:
    return output_root(explicit) / "creta_extractions"


def cafe_dir(explicit=None) -> Path:
    return output_root(explicit) / "cafe_output"
