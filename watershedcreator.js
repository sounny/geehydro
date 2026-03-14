// Earth Engine watershed creator with optional pour point snapping and basin metrics

// DEM options (first band used for routing calculations)
var demList = {
  'HydroSHEDS 03VFDEM': ee.Image('WWF/HydroSHEDS/03VFDEM'),
  'JAXA AW3D30 (mosaic)': ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1').mosaic(),
  'NASA ASTER GED': ee.Image('NASA/ASTER_GED/AG100_003'),
  'Copernicus GLO-30 (mosaic)': ee.ImageCollection('COPERNICUS/DEM/GLO30').mosaic(),
  'USGS SRTMGL1': ee.Image('USGS/SRTMGL1_003')
};

// Compute D8 flow direction from a single-band DEM.
// Each pixel is assigned one of the 8 powers-of-two D8 direction codes.
function computeD8FlowDirection(dem) {
  var band = ee.String(dem.bandNames().get(0));
  var neighborhood = dem.neighborhoodToBands(ee.Kernel.square(1));
  var center = neighborhood.select(band.cat('_0_0'));

  // Offsets are evaluated in row/column order around the center pixel.
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var directions = [32, 64, 128, 16, 1, 8, 4, 2];
  var distances = [Math.SQRT2, 1, Math.SQRT2, 1, 1, Math.SQRT2, 1, Math.SQRT2];

  var maxSlope = ee.Image(-9999);
  var flowDir = ee.Image(0).byte();

  for (var i = 0; i < 8; i++) {
    // Build neighborhood band names (example: elevation_-1_1) dynamically.
    var row = offsets[i][0];
    var col = offsets[i][1];
    var bandName = band.cat('_').cat(ee.Number(row).format('%d')).cat('_').cat(ee.Number(col).format('%d'));
    var neighbor = neighborhood.select(bandName);

    // Slope = elevation drop / cell-to-cell distance; only downslope values can route flow.
    var slope = center.subtract(neighbor).divide(distances[i]);
    var dirCode = ee.Image.constant(directions[i]).byte();
    var shouldReplace = slope.gt(maxSlope).and(slope.gt(0));

    // Track steepest downslope neighbor and corresponding D8 direction code.
    maxSlope = maxSlope.where(shouldReplace, slope);
    flowDir = flowDir.where(shouldReplace, dirCode);
  }

  return flowDir.rename('flowDirection');
}

// Delineate the upstream drainage area by repeatedly adding cells that flow into
// already accepted watershed cells.
function delineateWatershed(flowDir, outletMask, iterations) {
  var maxIterations = iterations || 120;
  var watershed = outletMask.rename('watershed').selfMask();
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var inflowDirections = [2, 4, 8, 1, 16, 128, 64, 32];

  for (var step = 0; step < maxIterations; step++) {
    // Neighborhood bands expose each surrounding pixel relative to current pixel.
    var wsNeigh = watershed.unmask(0).neighborhoodToBands(ee.Kernel.square(1));
    var dirNeigh = flowDir.neighborhoodToBands(ee.Kernel.square(1));
    var newCells = ee.Image(0).byte();

    for (var j = 0; j < 8; j++) {
      // If a neighbor routes toward the center and that neighbor is already in
      // the watershed, then the center should be included too.
      var row = offsets[j][0];
      var col = offsets[j][1];
      var wsBand = 'watershed_' + row + '_' + col;
      var dirBand = 'flowDirection_' + row + '_' + col;
      var routesIn = dirNeigh.select(dirBand).eq(inflowDirections[j]);
      var neighborIsWs = wsNeigh.select(wsBand).eq(1);
      newCells = newCells.or(routesIn.and(neighborIsWs));
    }

    watershed = watershed.unmask(0).or(newCells).rename('watershed').selfMask();
  }

  return watershed;
}

// Compute a simple D8 flow accumulation image by iterative upstream propagation.
function computeFlowAccumulation(flowDir, iterations) {
  var maxIterations = iterations || 120;
  var accumulation = ee.Image(1).toFloat();
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var inflowDirections = [2, 4, 8, 1, 16, 128, 64, 32];

  for (var step = 0; step < maxIterations; step++) {
    var accNeigh = accumulation.neighborhoodToBands(ee.Kernel.square(1));
    var dirNeigh = flowDir.neighborhoodToBands(ee.Kernel.square(1));
    var inflowSum = ee.Image(0).toFloat();

    for (var j = 0; j < 8; j++) {
      // Add the neighbor accumulation only when that neighbor drains into center.
      var row = offsets[j][0];
      var col = offsets[j][1];
      var accBand = 'constant_' + row + '_' + col;
      var dirBand = 'flowDirection_' + row + '_' + col;
      var routesIn = dirNeigh.select(dirBand).eq(inflowDirections[j]);
      var neighborAcc = accNeigh.select(accBand).multiply(routesIn);
      inflowSum = inflowSum.add(neighborAcc);
    }

    // A cell contributes its own area (+1) in addition to all routed inflow.
    accumulation = inflowSum.add(1).toFloat();
  }

  return accumulation.rename('flowAccumulation');
}

// Snap pour point to highest accumulation cell within radius to improve outlet placement.
function snapPourPoint(pointGeom, accumulation, radiusMeters) {
  var radius = ee.Number(radiusMeters || 1500);
  var searchArea = pointGeom.buffer(radius);

  // Find the largest accumulation value near the clicked outlet.
  var maxAcc = ee.Number(accumulation.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: searchArea,
    scale: 90,
    maxPixels: 1e8
  }).get('flowAccumulation'));

  // Convert the max-value pixel(s) to a point feature and return first match.
  var snappedFeature = accumulation.eq(maxAcc).selfMask().reduceToVectors({
    geometry: searchArea,
    scale: 90,
    geometryType: 'centroid',
    reducer: ee.Reducer.first(),
    maxPixels: 1e8
  }).first();

  return ee.Algorithms.If(snappedFeature, ee.Feature(snappedFeature).geometry(), pointGeom);
}

var demNames = Object.keys(demList);
var demSelect = ui.Select({items: demNames, value: demNames[0]});
var iterationsSlider = ui.Slider({min: 40, max: 220, value: 120, step: 10});
var snapToggle = ui.Checkbox({label: 'Snap pour point to highest accumulation', value: true});
var snapRadiusBox = ui.Textbox({value: '1500', placeholder: 'Snap radius meters'});
var runButton = ui.Button('Create Watershed');

var panel = ui.Panel([
  ui.Label('Watershed Creator', {fontWeight: 'bold'}),
  ui.Label('1) Draw one point outlet with the point drawing tool.'),
  ui.Label('2) Select DEM and click Create Watershed.'),
  ui.Label('DEM:'),
  demSelect,
  ui.Label('Iterations (higher handles larger basins):'),
  iterationsSlider,
  snapToggle,
  ui.Label('Snap radius (meters):'),
  snapRadiusBox,
  runButton
]);
panel.style().set({position: 'top-right', width: '360px'});
ui.root.add(panel);

Map.setOptions('SATELLITE');
Map.drawingTools().setDrawModes(['point']);
Map.drawingTools().draw();

runButton.onClick(function() {
  var layers = Map.drawingTools().layers();
  if (layers.length() === 0) {
    ui.alert('Draw one outlet point before creating the watershed.');
    return;
  }

  var outlet = ee.Feature(layers.get(0)).geometry();
  var dem = demList[demSelect.getValue()];
  var iterations = iterationsSlider.getValue();
  var snapRadius = parseInt(snapRadiusBox.getValue(), 10);
  snapRadius = isNaN(snapRadius) ? 1500 : snapRadius;

  var flowDir = computeD8FlowDirection(dem);
  var accumulation = computeFlowAccumulation(flowDir, iterations);

  // Optionally adjust outlet to a nearby flow path so delineation follows channels.
  var snappedOutlet = ee.Geometry(ee.Algorithms.If(
    snapToggle.getValue(),
    snapPourPoint(outlet, accumulation, snapRadius),
    outlet
  ));

  var outletMask = ee.Image().byte().paint(snappedOutlet, 1);
  var watershed = delineateWatershed(flowDir, outletMask, iterations);

  // Estimate watershed area in square kilometers for quick QA.
  var areaKm2 = ee.Image.pixelArea().divide(1e6).updateMask(watershed).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: snappedOutlet.buffer(500000),
    scale: 90,
    maxPixels: 1e9
  }).get('area');

  Map.clear();
  Map.addLayer(dem, {min: 0, max: 3000}, 'DEM');
  Map.addLayer(accumulation.log10(), {min: 0, max: 6, palette: ['081d58', '225ea8', '1d91c0', '41b6c4', 'c7e9b4']}, 'Log10 Flow Accumulation');
  Map.addLayer(watershed.selfMask(), {palette: ['2c7fb8']}, 'Watershed');
  Map.addLayer(outlet, {color: 'ff0000'}, 'Original Outlet');
  Map.addLayer(snappedOutlet, {color: '00ff00'}, 'Snapped Outlet');
  Map.centerObject(snappedOutlet, 9);

  print('Watershed area (km^2):', areaKm2);
});
