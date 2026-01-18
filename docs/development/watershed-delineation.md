Watershed delineation: missing points
====================================

Problem
-------
In cumulative mode, points are snapped with Whitebox `snap_pour_points`. Any point
that does not snap within the snap distance is skipped, so no watershed is
computed for it. The boundaries step then iterates only the snapped pour points,
so those points never appear downstream.

Recommendations
---------------
1) Fallback to original point geometry when snapping fails
   - In `src/delineator/hydrology/watershed.py`, if a snapped point is not created,
     run `watershed` using the original (unsnapped) point instead of skipping.
   - This ensures every study point produces some tributary area.

2) Make snapping more forgiving
   - Increase the default `snap_distance` (or make it user-configurable per run)
     so points are less likely to be dropped.

3) Add explicit reporting for dropped points
   - If snapping fails and you choose to keep it mandatory, log and emit a report
     (or warning list) of the point IDs/names that were skipped so the user can
     adjust inputs or snap distance.
