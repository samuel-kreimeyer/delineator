# Understanding Results

This guide explains how to interpret delineator's output files and reports.

## Output Files Overview

After running delineator, your output directory will contain:

```
output/
├── dem.tif              # Raw DEM (if --keep-intermediates)
├── conditioned_dem.tif  # Depression-handled DEM (if --keep-intermediates)
├── flow_dir.tif         # Flow direction (if --keep-intermediates)
├── flow_acc.tif         # Flow accumulation (if --keep-intermediates)
├── watershed.tif        # Watershed raster (if --keep-intermediates)
├── pour_points.shp      # Snapped pour points
├── watersheds.shp       # Watershed boundaries (or .geojson/.gpkg)
└── report.json          # Statistics report (or .csv/.txt)
```

## Watershed Boundaries

The `watersheds.shp` (or equivalent) file contains polygon features representing each watershed.

### Attributes

| Field | Description |
|-------|-------------|
| `watershed_id` | Unique identifier matching the pour point ID |
| `point_name` | Name of the pour point |
| `area_sqft` | Watershed area in square feet |
| `area_acres` | Watershed area in acres |
| `area_sqmi` | Watershed area in square miles |
| `elev_min` | Minimum elevation in the watershed (ft) |
| `elev_max` | Maximum elevation in the watershed (ft) |
| `elev_mean` | Mean elevation in the watershed (ft) |

### Cumulative vs Exclusive Watersheds

**Cumulative mode (default)**: Watersheds represent the total contributing area to each pour point. If pour points are nested (one upstream of another), watersheds will overlap.

```
Pour Point A (upstream)     Pour Point B (downstream)
     ↓                           ↓
  [====]                    [==========]
  Watershed A               Watershed B (includes A)
```

**Exclusive mode** (`--exclusive`): Each cell belongs to only one watershed. Upstream areas are "carved out" of downstream watersheds.

```
Pour Point A (upstream)     Pour Point B (downstream)
     ↓                           ↓
  [====]                    [    ======]
  Watershed A               Watershed B (excludes A)
```

## Pour Points

The `pour_points.shp` file contains the **snapped** pour point locations.

### Why Points Move

Pour points are snapped to nearby cells with high flow accumulation to ensure they land on actual drainage channels rather than flat areas or ridges.

### Attributes

| Field | Description |
|-------|-------------|
| `id` | Point identifier |
| `name` | Point name from input |

### Checking Snap Quality

If a watershed has unexpectedly small area, the pour point may not have snapped to a real channel. Solutions:

1. Increase `--snap-distance`
2. Decrease `--cell-size` for better channel resolution
3. Check if the original point location is correct

## Statistics Report

### JSON Format

```json
{
  "metadata": {
    "generated_at": "2024-01-15T10:30:00",
    "total_watersheds": 5
  },
  "summary": {
    "total_area_sqft": 4363527.56,
    "total_area_acres": 100.17,
    "total_area_sqmi": 0.156,
    "min_elevation_ft": 333.52,
    "max_elevation_ft": 368.46
  },
  "watersheds": [
    {
      "watershed_id": 1,
      "point_name": "Outlet_1",
      "area_sqft": 11037.67,
      "area_acres": 0.25,
      "area_sqmi": 0.0004,
      "elev_min_ft": 355.27,
      "elev_max_ft": 362.50,
      "elev_mean_ft": 359.39
    }
  ]
}
```

### CSV Format

```csv
watershed_id,point_name,area_sqft,area_acres,area_sqmi,elev_min,elev_max,elev_mean
1,Outlet_1,11037.67,0.25,0.0004,355.27,362.50,359.39
```

### Text Format

```
============================================================
WATERSHED DELINEATION REPORT
============================================================

Generated: 2024-01-15T10:30:00
Total watersheds: 5

------------------------------------------------------------
SUMMARY STATISTICS
------------------------------------------------------------

Total area: 0.1565 mi²
           100.17 acres
           4363528 ft²

Elevation range: 333.52 - 368.46 ft

------------------------------------------------------------
INDIVIDUAL WATERSHEDS
------------------------------------------------------------

Watershed 1: Outlet_1
  Area: 0.25 acres (0.0004 mi²)
  Elevation: min=355.27 ft, max=362.50 ft, mean=359.39 ft
```

## Intermediate Rasters

When using `--keep-intermediates`, these rasters are preserved:

### dem.tif
Raw DEM interpolated from the TIN surface. Useful for:
- Visualizing the terrain
- Checking interpolation quality
- Identifying data gaps (NoData areas)

### conditioned_dem.tif
DEM after depression handling. Compare with `dem.tif` to see where breaching/filling occurred.

### flow_dir.tif
D8 flow direction grid. Values indicate the direction water flows from each cell:
- 1 = East
- 2 = Northeast
- 4 = North
- 8 = Northwest
- 16 = West
- 32 = Southwest
- 64 = South
- 128 = Southeast

### flow_acc.tif
Flow accumulation grid. Each cell's value represents the number of upstream cells that drain to it. High values indicate streams/channels.

### watershed.tif
Raster showing watershed membership. Each cell's value is the ID of the watershed it belongs to.

## Visualizing Results

### In QGIS

1. Add `dem.tif` as base layer with hillshade styling
2. Add `watersheds.shp` with semi-transparent fill
3. Add `pour_points.shp` as markers
4. Optionally add `flow_acc.tif` with logarithmic scaling to see drainage network

### In ArcGIS

1. Add DEM and apply hillshade
2. Add watershed polygons with unique values symbology
3. Add pour points with labels

## Common Issues

### Watershed Too Small

**Symptoms**: Watershed area is near zero or just a few cells.

**Causes**:
1. Pour point not snapping to drainage channel
2. Grid resolution too coarse to capture narrow channels
3. Pour point located on a ridge or flat area

**Solutions**:
1. Increase `--snap-distance`
2. Decrease `--cell-size`
3. Verify pour point coordinates

### Unexpected Watershed Shapes

**Symptoms**: Watershed boundary doesn't follow expected terrain.

**Causes**:
1. Depression handling altered flow paths
2. TIN surface has errors or gaps
3. Coordinate system mismatch

**Solutions**:
1. Compare `dem.tif` with `conditioned_dem.tif`
2. Check source TIN for issues
3. Verify EPSG code matches your data

### Missing Watersheds

**Symptoms**: Some pour points have no watershed output.

**Causes**:
1. Pour point outside TIN extent
2. Pour point in NoData area
3. Snapping failed

**Solutions**:
1. Check pour point coordinates against TIN bounds
2. Increase snap distance
3. Review verbose output for warnings
