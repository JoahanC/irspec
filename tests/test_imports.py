"""Smoke tests: the package imports from any cwd and resources resolve."""
import importlib

import pytest

MODULES = ["datacube", "plotparams", "fitfuncs", "spec_helpers", "cube_stitch",
           "emission_io", "spaxel_fit", "spaxel_fit_render", "sort_dynamic",
           "map_ratios", "false_color", "paths", "cubespec"]


@pytest.mark.parametrize("mod", MODULES)
def test_module_imports(mod):
    # cubespec must import even without CAFE/CRETA installed (lazy imports).
    importlib.import_module(f"irspec.{mod}")


def test_line_params_resource():
    from irspec.emission_io import (read_line_params, read_line_params2,
                                    read_line_params3)
    assert read_line_params() and read_line_params2() and read_line_params3()
