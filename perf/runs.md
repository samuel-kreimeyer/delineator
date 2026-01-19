# Run Log

## Test cases
- tests/data/sample.xml (use COGO points)
- tests/data/080724_Drainage.xml (use COGO points)

## Baseline
- Date:
- Command:
  - `python -m delineator tests/data/sample.xml -o perf/runs/sample --use-cogo-points`
  - `python -m delineator tests/data/080724_Drainage.xml -o perf/runs/drainage --use-cogo-points`
- Timing:
- Output compare:
  - `diff -u tests/output/report.json perf/runs/sample/report.json`
  - For rasters, compare array values (example):
    - `python - <<'PY'`
      `import rasterio, numpy as np`
      `a=rasterio.open('tests/output/dem.tif').read(1)`
      `b=rasterio.open('perf/runs/sample/dem.tif').read(1)`
      `print(np.array_equal(a,b))`
      `PY`
- Notes:
  - `tests/data/sample.xml` failed in cumulative mode: "No valid watershed polygons generated" (see `perf/flamegraphs/sample.svg` for partial capture).
  - `tests/data/080724_Drainage.xml` produced `perf/flamegraphs/drainage.svg`.

## After changes
- Date: 2026-01-18
- Command:
  - `.venv/bin/python -m delineator tests/data/080724_Drainage.xml -o perf/runs/drainage --use-cogo-points --keep-intermediates`
  - `/usr/bin/time -p .venv/bin/python -m delineator tests/data/080724_Drainage.xml -o perf/runs/drainage_iter --use-cogo-points`
- Timing: ~76s wall time (single run, prior to iterparse change)
- Timing (iterparse): real 86.90s, user 86.54s, sys 3.71s (includes full run; not directly comparable to earlier rough timing)
- Timing (tree parse baseline, same command): real 96.56s, user 93.76s, sys 5.65s
- Chunked rasterization validation:
  - Command: `perf/compare_runs.py --input tests/data/080724_Drainage.xml --output-dir perf/runs/drainage_chunk --reference-dir tests/output --use-cogo-points --cell-size 1.0 --snap-distance 10.0 --keep-intermediates`
  - Result: all outputs match reference with 2-decimal float tolerance; shapefile DBF includes extra columns `lfp_length`, `lfp_slope` not present in reference.
- Chunked rasterization A/B timing (same command, `--cell-size 1.0 --snap-distance 10.0`):
  - Non-chunked (legacy): real 156.66s, user 157.75s, sys 11.44s
  - Chunked (current): real 154.69s, user 156.53s, sys 10.98s
- Output compare:
  - `perf/compare_runs.py` vs `tests/output` reported raster profile differences and shapefile/report diffs; likely due to dependency/tooling differences (GDAL/Whitebox/pyogrio versions).
- Notes:
  - Flamegraph: `perf/flamegraphs/drainage.svg`
  - Matching `tests/output` required `--cell-size 1.0` and `--snap-distance 10.0`; with those flags, report/summary matched except for a floating-point rounding difference in `elev_mean_ft`.
