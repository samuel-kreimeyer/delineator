# Performance Findings

## Implemented
- src/delineator/output/boundaries.py: vectorized watershed statistics computation to avoid per-watershed full-array scans.
  - Uses numpy grouping to compute area and elevation stats in one pass.
  - Expected to improve exclusive-mode performance when many watersheds exist.
- src/delineator/parsers/landxml.py: switched to iterparse streaming to avoid building a full XML tree.
  - Lowers memory usage and reduces parsing overhead for large LandXML files.

## Flamegraph highlights (tests/data/080724_Drainage.xml)
- perf/flamegraphs/drainage.svg shows the biggest CPU share in TIN rasterization setup:
  - matplotlib triangulation trifinder initialization (`matplotlib/tri/_trifinder.py`) ~39%.
  - LandXML parsing of faces (`delineator/parsers/landxml.py:_parse_faces`) ~5-6%.
  - TIN average edge length computation (`delineator/geometry/tin.py:compute_average_edge_length`) ~7%.
  - TIN interpolation (`matplotlib/tri/_triinterpolate.py:_interpolate_multikeys`) ~13-14%.

## Opportunities (not yet implemented)
- src/delineator/geometry/rasterize.py: large DEM interpolation allocates full grids; consider windowed/chunked interpolation to reduce peak memory and allow streaming writes.
- src/delineator/geometry/rasterize.py: consider bypassing `matplotlib.tri` for large TINs (e.g., CGAL/scipy spatial Delaunay + custom interpolator or Whitebox/PDAL gridding) to avoid heavy trifinder setup costs seen in flamegraph.
- src/delineator/parsers/landxml.py: for huge LandXML files, `iterparse` could reduce peak XML memory and speed up point/face parsing.
- src/delineator/hydrology/watershed.py and src/delineator/output/boundaries.py: cumulative mode runs Whitebox once per point; parallelizing per-point runs (if thread-safe) or batching could reduce wall-clock time.
- src/delineator/geometry/tin.py: cache point arrays/triangle index mapping after parsing to avoid repeated sorting on large TINs.
