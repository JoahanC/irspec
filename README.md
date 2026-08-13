# irspec

## Install

```
pip install -e .
```

CAFE/CRETA (spectral fitting) are not on PyPI; install separately with
`pip install -e /path/to/CAFE-master`. The q3dfit backend is optional:
`pip install irspec[q3dfit]`.

## Data locations

Science data lives outside the package. Point irspec at it with the
`IRSPEC_DATA_ROOT` environment variable (a directory containing
`image_data/`, `cubes/`, `cubes_north/`, `param_files/`) and set
`IRSPEC_OUTPUT_ROOT` for where products (`spaxs/`, `plots/`,
`cafe_output/`, `creta_extractions/`) are written. Alternatively put
`data_root` / `output_root` in `~/.config/irspec/config.toml`. The pixi
workspace sets both variables automatically.

## Command-line tools

`irspec-run-creta`, `irspec-run-cafe`, `irspec-spaxel-map`,
`irspec-footprints`, `irspec-plot-line-maps`, `irspec-plot-pah-maps`,
`irspec-plot-alma`, `irspec-spectral-index`, `irspec-false-color` —
all accept `--help`.
