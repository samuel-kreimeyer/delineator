# Configuration Guide

This guide covers all command-line options available in delineator.

## Command Syntax

```bash
delineator <input_file> -o <output_dir> [options]
```

## Required Arguments

| Argument | Description |
|----------|-------------|
| `input_file` | Path to LandXML file containing TIN surface |
| `-o, --output-dir` | Output directory for results |

## Study Points (one required)

You must specify study points using one of these options:

| Option | Description |
|--------|-------------|
| `--use-cogo-points` | Extract study points from COGO points in LandXML |
| `--study-points FILE` | Load study points from external file (CSV, GeoJSON, Shapefile, GeoPackage) |

## Processing Options

### Cell Size

```bash
--cell-size <float>
```

Raster cell size in map units (typically feet or meters, depending on your coordinate system).

- **Default**: Auto-computed as half the average TIN triangle edge length
- **Recommendation**: Use smaller values (e.g., 1 ft) for better resolution of narrow features like ditches

Example:
```bash
delineator input.xml -o output/ --use-cogo-points --cell-size 1
```

### Snap Distance

```bash
--snap-distance <float>
```

Maximum distance (in map units) to snap pour points to the nearest high-flow-accumulation cell.

- **Default**: 2x the cell size
- **Purpose**: Ensures pour points land on actual drainage channels, not flat areas
- **Recommendation**: Increase if pour points have very small watersheds (may indicate they're not snapping to channels)

Example:
```bash
delineator input.xml -o output/ --use-cogo-points --snap-distance 10
```

### Depression Handling Method

```bash
--depression-method {fill,breach,breach-least-cost}
```

Method for handling depressions (sinks) in the DEM before flow routing.

| Method | Description | Speed | Quality |
|--------|-------------|-------|---------|
| `fill` | Fill depressions to overflow level | Fast | Lower |
| `breach` | Carve channels through barriers (default) | Medium | Good |
| `breach-least-cost` | Optimal breach paths | Slow | Best |

- **Default**: `breach`
- **Recommendation**: Use `breach` for most cases. Use `breach-least-cost` for final production runs where accuracy is critical.

### Watershed Mode

```bash
--exclusive
```

By default, delineator computes **cumulative watersheds** where each watershed includes all upstream contributing area. Watersheds can overlap if pour points are nested along the same flow path.

Use `--exclusive` to compute non-overlapping watersheds where each cell is assigned to only one watershed (the first pour point it drains to).

| Mode | Overlapping | Use Case |
|------|-------------|----------|
| Cumulative (default) | Yes | Standard drainage analysis |
| Exclusive | No | When you need non-overlapping boundaries |

## Output Options

### Output Format

```bash
--output-format {shapefile,geojson,geopackage}
```

Format for watershed boundary output.

- **Default**: `shapefile`
- **Options**:
  - `shapefile` - ESRI Shapefile (.shp)
  - `geojson` - GeoJSON (.geojson)
  - `geopackage` - OGC GeoPackage (.gpkg)

### Report Format

```bash
--report-format {json,csv,text,all}
```

Format for watershed statistics report.

- **Default**: `json`
- **Options**:
  - `json` - Machine-readable JSON
  - `csv` - Spreadsheet-compatible CSV
  - `text` - Human-readable plain text
  - `all` - Generate all three formats

### Keep Intermediate Files

```bash
--keep-intermediates
```

Retain intermediate raster files for debugging or visualization:
- `conditioned_dem.tif` - Depression-handled DEM
- `flow_dir.tif` - D8 flow direction grid
- `flow_acc.tif` - Flow accumulation grid
- `watershed.tif` - Watershed raster

## Other Options

### EPSG Override

```bash
--epsg <code>
```

Override the coordinate system EPSG code (normally read from LandXML).

Example:
```bash
delineator input.xml -o output/ --use-cogo-points --epsg 3441
```

### Verbose Output

```bash
-v, --verbose
```

Enable detailed logging output for debugging.

## Complete Example

```bash
delineator survey.xml -o results/ \
    --use-cogo-points \
    --cell-size 1 \
    --snap-distance 10 \
    --depression-method breach \
    --output-format geopackage \
    --report-format all \
    --keep-intermediates \
    --verbose
```

## Units

Delineator uses **imperial units** for all output:

| Measurement | Units |
|-------------|-------|
| Area | Square feet (ft²), Acres, Square miles (mi²) |
| Elevation | Feet (ft) |
| Cell size | Map units (typically feet) |
| Snap distance | Map units (typically feet) |

The coordinate system units are determined by the EPSG code in your LandXML file.
