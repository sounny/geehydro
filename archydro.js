// Earth Engine Arc Hydro-style workflow for hydrologic modeling.
//
// This script extends the existing D8 tools into a more complete Arc Hydro
// chain: DEM fill, flow direction, flow accumulation, stream definition,
// stream segmentation, Strahler order, drainage lines, and catchment polygons.

// DEM options (first band used for analysis).
var demList = {
  'HydroSHEDS 03VFDEM': ee.Image('WWF/HydroSHEDS/03VFDEM'),
  'JAXA AW3D30 (mosaic)': ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1').mosaic(),
  'NASA ASTER GED': ee.Image('NASA/ASTER_GED/AG100_003'),
  'Copernicus GLO-30 (mosaic)': ee.ImageCollection('COPERNICUS/DEM/GLO30').mosaic(),
  'USGS SRTMGL1': ee.Image('USGS/SRTMGL1_003'),
  'USGS GTOPO30': ee.Image('USGS/GTOPO30'),
  'NOAA ETOPO1': ee.Image('NOAA/NGDC/ETOPO1')
};

// Compute D8 flow direction from a single-band DEM.
// Direction coding follows Arc Hydro style powers of 2.
function d8FlowDirection(dem) {
  var band = ee.String(dem.bandNames().get(0));
  var neighborhood = dem.neighborhoodToBands(ee.Kernel.square(1));
  var center = neighborhood.select(band.cat('_0_0'));

  // Analyze all 8 neighboring cells around the center pixel.
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
    var row = offsets[i][0];
    var col = offsets[i][1];
    var bandName = band
      .cat('_')
      .cat(ee.Number(row).format('%d'))
      .cat('_')
      .cat(ee.Number(col).format('%d'));

    // Positive slope means the neighbor is lower and can receive flow.
    var neighbor = neighborhood.select(bandName);
    var slope = center.subtract(neighbor).divide(distances[i]);
    var dirImg = ee.Image.constant(directions[i]).byte();
    var replace = slope.gt(maxSlope).and(slope.gt(0));

    maxSlope = maxSlope.where(replace, slope);
    flowDir = flowDir.where(replace, dirImg);
  }

  return flowDir.rename('flowDirection');
}

// Fill local pits by repeatedly applying a focal maximum.
// This is a simplified depression-fill to enforce downslope routing.
function fillDem(dem, iterations) {
  iterations = iterations || 12;
  var filled = dem;
  for (var i = 0; i < iterations; i++) {
    var focalMax = filled.focal_max({radius: 1, units: 'pixels'});
    filled = filled.where(filled.lt(focalMax), focalMax);
  }
  return filled.rename('filledDem');
}

// Compute flow accumulation by iterative upstream contribution transfer.
function flowAccumulation(flowDir, iterations) {
  iterations = iterations || 120;

  // Start with a unit contribution per cell, then route it downstream.
  var acc = ee.Image.constant(1).rename('flowAccum');
  var kernel = ee.Kernel.square(1);
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var inflowDir = [2, 4, 8, 1, 16, 128, 64, 32];

  for (var step = 0; step < iterations; step++) {
    var accNeigh = acc.neighborhoodToBands(kernel);
    var dirNeigh = flowDir.neighborhoodToBands(kernel);
    var contrib = ee.Image(0);

    // Add each neighbor's accumulation only if it drains into center pixel.
    for (var j = 0; j < 8; j++) {
      var row = offsets[j][0];
      var col = offsets[j][1];
      var dirBand = 'flowDirection_' + row + '_' + col;
      var accBand = 'flowAccum_' + row + '_' + col;
      var mask = dirNeigh.select(dirBand).eq(inflowDir[j]);
      contrib = contrib.add(accNeigh.select(accBand).updateMask(mask));
    }

    acc = acc.add(contrib);
  }

  return acc;
}

// Build stream links by assigning an ID to each connected stream segment.
function streamLinks(streams) {
  return streams
    .selfMask()
    .connectedComponents(ee.Kernel.plus(1), 1024)
    .rename('streamLink');
}

// Compute Strahler stream order by counting neighbor stream junctions.
function strahlerOrder(streams, flowDir, iterations) {
  iterations = iterations || 60;
  var order = streams.selfMask().multiply(0).add(1).rename('streamOrder');
  var kernel = ee.Kernel.square(1);

  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var inflowDir = [2, 4, 8, 1, 16, 128, 64, 32];

  for (var i = 0; i < iterations; i++) {
    var streamNeigh = streams.neighborhoodToBands(kernel);
    var orderNeigh = order.unmask(0).neighborhoodToBands(kernel);
    var dirNeigh = flowDir.neighborhoodToBands(kernel);

    var maxInOrder = ee.Image(0);
    var sameCount = ee.Image(0);

    // Gather orders from upstream stream neighbors that flow into center.
    for (var j = 0; j < 8; j++) {
      var row = offsets[j][0];
      var col = offsets[j][1];
      var streamBand = 'constant_' + row + '_' + col;
      var orderBand = 'streamOrder_' + row + '_' + col;
      var dirBand = 'flowDirection_' + row + '_' + col;
      var valid = streamNeigh.select(streamBand).eq(1).and(dirNeigh.select(dirBand).eq(inflowDir[j]));
      var upstreamOrder = orderNeigh.select(orderBand).updateMask(valid);

      maxInOrder = maxInOrder.max(upstreamOrder.unmask(0));
      sameCount = sameCount.add(upstreamOrder.eq(maxInOrder).updateMask(valid).unmask(0));
    }

    // Strahler rule: increment only when 2+ upstream tributaries share max order.
    var updated = maxInOrder.where(sameCount.gte(2), maxInOrder.add(1));
    order = order.where(streams.eq(1).and(updated.gt(order)), updated);
  }

  return order.rename('streamOrder');
}

// Propagate stream IDs upstream to create catchment grids.
function catchmentGrid(flowDir, streams, iterations) {
  iterations = iterations || 140;
  var links = streamLinks(streams);
  var basins = links.unmask(0).rename('catchment');
  var kernel = ee.Kernel.square(1);

  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var inflowDir = [2, 4, 8, 1, 16, 128, 64, 32];

  for (var step = 0; step < iterations; step++) {
    var basinNeigh = basins.neighborhoodToBands(kernel);
    var dirNeigh = flowDir.neighborhoodToBands(kernel);
    var newPix = ee.Image(0);

    // Transfer stream/catchment IDs from downstream neighbor to upstream pixels.
    for (var j = 0; j < 8; j++) {
      var row = offsets[j][0];
      var col = offsets[j][1];
      var dirBand = 'flowDirection_' + row + '_' + col;
      var basinBand = 'catchment_' + row + '_' + col;
      var mask = dirNeigh.select(dirBand).eq(inflowDir[j]);
      var contrib = basinNeigh.select(basinBand).updateMask(mask);
      newPix = newPix.where(contrib.neq(0), contrib);
    }

    basins = basins.where(basins.eq(0).and(newPix.neq(0)), newPix);
  }

  return basins.rename('catchment');
}

// Compute downstream flow length by walking flow directions iteratively.
function downstreamFlowLength(flowDir, iterations) {
  iterations = iterations || 150;
  var length = ee.Image(0).rename('flowLength');
  var kernel = ee.Kernel.square(1);

  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [0, -1],            [0, 1],
    [1, -1],  [1, 0],   [1, 1]
  ];
  var outflowDir = [32, 64, 128, 16, 1, 8, 4, 2];
  var stepDist = [Math.SQRT2, 1, Math.SQRT2, 1, 1, Math.SQRT2, 1, Math.SQRT2];

  for (var step = 0; step < iterations; step++) {
    var lenNeigh = length.neighborhoodToBands(kernel);
    var candidate = ee.Image(0);

    // For each direction, pull the downstream pixel length and add step distance.
    for (var j = 0; j < 8; j++) {
      var row = offsets[j][0];
      var col = offsets[j][1];
      var dirCode = outflowDir[j];
      var bandName = 'flowLength_' + row + '_' + col;
      var routed = lenNeigh.select(bandName).add(stepDist[j]);
      candidate = candidate.max(routed.updateMask(flowDir.eq(dirCode)).unmask(0));
    }

    length = length.max(candidate);
  }

  return length.rename('flowLength');
}

// User interface controls.
var demNames = Object.keys(demList);
var demSelect = ui.Select({items: demNames, value: demNames[0]});
var streamThreshold = ui.Slider({min: 100, max: 15000, value: 1500, step: 100, style: {stretch: 'horizontal'}});
var fillIterations = ui.Slider({min: 4, max: 30, value: 12, step: 1, style: {stretch: 'horizontal'}});
var routeIterations = ui.Slider({min: 40, max: 260, value: 140, step: 10, style: {stretch: 'horizontal'}});
var runBtn = ui.Button('Run Arc Hydro Workflow');

var panel = ui.Panel([
  ui.Label('Arc Hydro-style D8 Modeling', {fontWeight: 'bold'}),
  ui.Label('Select DEM:'), demSelect,
  ui.Label('Depression fill iterations:'), fillIterations,
  ui.Label('Routing iterations (accumulation/catchment/length):'), routeIterations,
  ui.Label('Stream definition threshold (cells):'), streamThreshold,
  ui.Label('Draw an AOI polygon, then run the workflow.'),
  runBtn
]);
panel.style().set({position: 'top-right', width: '340px'});
ui.root.add(panel);

Map.drawingTools().setDrawModes(['polygon']);
Map.drawingTools().draw();

runBtn.onClick(function() {
  var layers = Map.drawingTools().layers();
  if (layers.length() === 0) {
    ui.alert('Please draw an AOI polygon.');
    return;
  }

  var aoi = ee.Feature(layers.get(0).getEeObject()).geometry();
  var rawDem = demList[demSelect.getValue()].clip(aoi);
  var fillIters = fillIterations.getValue();
  var routeIters = routeIterations.getValue();

  // Execute Arc Hydro chain from terrain conditioning to vector products.
  var conditionedDem = fillDem(rawDem, fillIters).clip(aoi);
  var fdir = d8FlowDirection(conditionedDem).clip(aoi);
  var facc = flowAccumulation(fdir, routeIters).clip(aoi);
  var streams = facc.gte(streamThreshold.getValue()).rename('streams');
  var links = streamLinks(streams).clip(aoi);
  var order = strahlerOrder(streams, fdir, 50).clip(aoi);
  var catchment = catchmentGrid(fdir, streams, routeIters).clip(aoi);
  var flowLength = downstreamFlowLength(fdir, routeIters).clip(aoi);

  // Convert raster products to Arc Hydro-like vector outputs.
  var drainageLines = streams.selfMask().reduceToVectors({
    geometry: aoi,
    geometryType: 'polygon',
    scale: 90,
    labelProperty: 'stream',
    maxPixels: 1e8
  });

  var catchmentPolygons = catchment.selfMask().reduceToVectors({
    geometry: aoi,
    geometryType: 'polygon',
    scale: 90,
    labelProperty: 'catchment',
    maxPixels: 1e8
  });

  Map.clear();
  Map.centerObject(aoi, 10);
  Map.addLayer(conditionedDem, {min: 0, max: 3000, palette: ['0d0887', '7e03a8', 'cc4778', 'f89540', 'f0f921']}, 'Conditioned DEM');
  Map.addLayer(fdir, {min: 1, max: 128, palette: ['f7fbff', '6baed6', '08519c']}, 'Flow Direction');
  Map.addLayer(facc.log10(), {min: 0, max: 6, palette: ['ffffff', '9ecae1', '3182bd', '08519c']}, 'Flow Accumulation (log10)');
  Map.addLayer(streams.selfMask(), {palette: ['0000ff']}, 'Stream Raster');
  Map.addLayer(order.selfMask(), {min: 1, max: 6, palette: ['ffffcc', 'a1dab4', '41b6c4', '2c7fb8', '253494']}, 'Strahler Order');
  Map.addLayer(links.randomVisualizer(), {}, 'Stream Links');
  Map.addLayer(catchment.randomVisualizer(), {}, 'Catchment Grid');
  Map.addLayer(flowLength, {min: 0, max: 600, palette: ['edf8fb', '66c2a4', '238b45', '00441b']}, 'Downstream Flow Length');
  Map.addLayer(drainageLines.style({color: '00ffff', fillColor: '00000000', width: 1}), {}, 'Drainage Features');
  Map.addLayer(catchmentPolygons.style({color: 'ff8800', fillColor: '00000000', width: 1}), {}, 'Catchment Polygons');
  Map.addLayer(aoi, {color: 'ff0000'}, 'AOI');

  print('Arc Hydro products', {
    flowDirection: fdir,
    flowAccumulation: facc,
    streams: streams,
    streamLinks: links,
    streamOrder: order,
    catchmentGrid: catchment,
    downstreamFlowLength: flowLength,
    drainageFeatures: drainageLines,
    catchmentPolygons: catchmentPolygons
  });
});
