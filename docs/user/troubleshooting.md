# Troubleshooting

This guide covers common issues and their solutions.

## Installation Issues

### WhiteBox Tools Download Fails

**Error**: `Downloading WhiteboxTools pre-compiled binary...` hangs or fails.

**Solution**: WhiteBox Tools downloads its binary on first use. If this fails:

1. Check your internet connection
2. Try running with `--verbose` to see detailed errors
3. Manually download from [WhiteboxGeo](https://www.whiteboxgeo.com/) and place in the whitebox package directory

### Missing Dependencies

**Error**: `ModuleNotFoundError: No module named 'rasterio'`

**Solution**: Reinstall with all dependencies:

```bash
pip install -e ".[dev]"
```

Or install missing package:

```bash
pip install rasterio
```

## Input File Issues

### LandXML Parse Error

**Error**: `LandXMLParseError: No Surface elements found`

**Causes**:
- File is not valid LandXML
- Surface is in unexpected location
- Namespace issues

**Solutions**:
1. Verify file is valid XML: `xmllint --noout input.xml`
2. Check that file contains `<Surface>` element
3. Try opening in Civil 3D or other software to verify

### No COGO Points Found

**Error**: `No study points found` when using `--use-cogo-points`

**Cause**: LandXML file doesn't contain `<CgPoints>` section.

**Solutions**:
1. Check if file has COGO points: search for `<CgPoint` in the file
2. Use external study points file instead: `--study-points points.csv`

### Invalid Coordinate System

**Warning**: `No EPSG code found, using 0 (undefined CRS)`

**Cause**: LandXML doesn't specify coordinate system.

**Solutions**:
1. Add `--epsg <code>` to specify manually
2. Edit LandXML to add: `<CoordinateSystem epsgCode="3441"/>`

## Processing Issues

### Very Small Watersheds

**Symptom**: Watershed has area of only a few square feet.

**Cause**: Pour point not snapping to drainage channel.

**Diagnosis**:
```bash
# Keep intermediate files and check flow accumulation
delineator input.xml -o output/ --use-cogo-points --keep-intermediates -v
```

Then check flow accumulation at pour point location.

**Solutions**:
1. Increase snap distance:
   ```bash
   --snap-distance 20
   ```
2. Increase grid resolution:
   ```bash
   --cell-size 1
   ```
3. Both:
   ```bash
   --cell-size 1 --snap-distance 10
   ```

### WhiteBox Tools Fails Silently

**Symptom**: Process completes but output files are missing.

**Cause**: WhiteBox Tools requires absolute paths internally.

**Solution**: This should be fixed in current version. If you encounter this:
1. Use `--verbose` to see detailed output
2. Ensure output directory exists and is writable
3. Check disk space

### Out of Memory

**Symptom**: Process killed or `MemoryError`

**Cause**: Very large TIN or very fine cell size creating huge rasters.

**Solutions**:
1. Increase cell size:
   ```bash
   --cell-size 2
   ```
2. Process a smaller area
3. Use a machine with more RAM

### Slow Processing

**Symptom**: Process takes very long time.

**Causes**:
- Very fine cell size
- Large TIN surface
- `breach-least-cost` depression method

**Solutions**:
1. Use coarser cell size for initial testing
2. Use `--depression-method breach` instead of `breach-least-cost`
3. Process subset of data first

## Output Issues

### Shapefile Field Name Truncation

**Warning**: `Column names longer than 10 characters will be truncated`

**Cause**: ESRI Shapefile format limits field names to 10 characters.

**Solutions**:
1. Use GeoPackage format (no limit):
   ```bash
   --output-format geopackage
   ```
2. Use GeoJSON format:
   ```bash
   --output-format geojson
   ```
3. Ignore warning (data is still valid, just truncated names)

### Watershed Polygons Have Holes

**Symptom**: Watershed polygon has unexpected holes.

**Cause**: NoData areas in DEM or disconnected drainage areas.

**Solutions**:
1. Check DEM for gaps: visualize `dem.tif` in GIS software
2. Verify TIN coverage
3. May be correct if terrain has internal drainage

### Report Shows Wrong Units

**Note**: All output is in imperial units (feet, acres, square miles).

If your input is in meters, the numbers will be incorrect. Ensure your coordinate system uses feet (check EPSG code).

## Getting Help

### Enable Verbose Mode

Always run with `-v` when troubleshooting:

```bash
delineator input.xml -o output/ --use-cogo-points -v
```

### Keep Intermediate Files

Preserve intermediate files for debugging:

```bash
delineator input.xml -o output/ --use-cogo-points --keep-intermediates
```

### Check Intermediate Results

1. **dem.tif**: Does it look correct? Any gaps?
2. **flow_acc.tif**: Do streams appear where expected?
3. **pour_points.shp**: Did points snap to channels?

### Report Issues

If you encounter a bug:

1. Run with `--verbose` and capture output
2. Note your Python version: `python --version`
3. Note your OS and version
4. Report at: https://github.com/example/delineator/issues
