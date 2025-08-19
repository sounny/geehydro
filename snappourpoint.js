// Earth Engine snap pour point to highest accumulation within search radius

// DEM options (first band used)
var demList = {
  'HydroSHEDS 03VFDEM': ee.Image('WWF/HydroSHEDS/03VFDEM'),
  'JAXA AW3D30 (mosaic)': ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1').mosaic(),
  'NASA ASTER GED': ee.Image('NASA/ASTER_GED/AG100_003'),
  'Copernicus GLO-30 (mosaic)': ee.ImageCollection('COPERNICUS/DEM/GLO30').mosaic(),
  'USGS SRTMGL1': ee.Image('USGS/SRTMGL1_003'),
  'USGS GTOPO30': ee.Image('USGS/GTOPO30'),
  'NOAA ETOPO1': ee.Image('NOAA/NGDC/ETOPO1')
};

// D8 flow direction from a single-band DEM
function D8Algorithm(dem) {
  var band = ee.String(dem.bandNames().get(0));
  var neighborhood = dem.neighborhoodToBands(ee.Kernel.square(1));
  var center = neighborhood.select(band.cat('_0_0'));
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [ 0, -1],           [ 0, 1],
    [ 1, -1], [ 1, 0],  [ 1, 1]
  ];
  var directions = [32, 64, 128, 16, 1, 8, 4, 2];
  var distances = [Math.SQRT2, 1, Math.SQRT2, 1, 1, Math.SQRT2, 1, Math.SQRT2];
  var maxSlope = ee.Image(-9999);
  var flowDir = ee.Image(0).byte();
  for (var i = 0; i < 8; i++) {
    var row = offsets[i][0];
    var col = offsets[i][1];
    var bandName = band.cat('_').cat(ee.Number(row).format('%d')).cat('_').cat(ee.Number(col).format('%d'));
    var neighbor = neighborhood.select(bandName);
    var slope = center.subtract(neighbor).divide(distances[i]);
    var dirImg = ee.Image.constant(directions[i]).byte();
    var replace = slope.gt(maxSlope).and(slope.gt(0));
    maxSlope = maxSlope.where(replace, slope);
    flowDir = flowDir.where(replace, dirImg);
  }
  return flowDir.rename('flowDirection');
}

// Flow accumulation using an iterative D8 approach
function flowAccumulation(flowDir, aoi, iterations) {
  iterations = iterations || 100;
  var acc = ee.Image.constant(1).rename('flowAccum');
  var kernel = ee.Kernel.square(1);
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [ 0, -1],           [ 0, 1],
    [ 1, -1], [ 1, 0],  [ 1, 1]
  ];
  var inflowDir = [2, 4, 8, 1, 16, 128, 64, 32];
  for (var step = 0; step < iterations; step++) {
    var accNeigh = acc.neighborhoodToBands(kernel);
    var dirNeigh = flowDir.neighborhoodToBands(kernel);
    var contrib = ee.Image(0);
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
  return acc.clip(aoi);
}

// Snap pour point to highest flow accumulation within radius (m)
function snapPourPoint(acc, point, radius, scale) {
  var buffer = point.buffer(radius);
  var maxVal = ee.Number(acc.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: buffer,
    scale: scale,
    bestEffort: true
  }).get('flowAccum'));
  var maxMask = acc.eq(maxVal).selfMask();
  var lonLat = ee.Image.pixelLonLat().updateMask(maxMask);
  var coords = lonLat.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: buffer,
    scale: scale,
    bestEffort: true
  });
  return ee.Geometry.Point([coords.get('longitude'), coords.get('latitude')]);
}

// UI elements
var demNames = Object.keys(demList);
var demSelect = ui.Select({items: demNames, value: demNames[0]});
var radiusSlider = ui.Slider({min: 50, max: 1000, value: 200, step: 50, style: {stretch: 'horizontal'}});
var runBtn = ui.Button('Snap Pour Point');
var panel = ui.Panel([
  ui.Label('Select DEM:'), demSelect,
  ui.Label('Search radius (m):'), radiusSlider,
  ui.Label('Draw a pour point then click the button.'),
  runBtn
]);
panel.style().set({position: 'top-right', width: '300px'});
ui.root.add(panel);

Map.drawingTools().setDrawModes(['point']);
Map.drawingTools().draw();

runBtn.onClick(function() {
  var layers = Map.drawingTools().layers();
  if (layers.length() === 0) {
    ui.alert('Please draw a pour point.');
    return;
  }
  var pt = ee.Feature(layers.get(0)).geometry();
  var dem = demList[demSelect.getValue()];
  var area = pt.buffer(radiusSlider.getValue() * 5);
  var fdir = D8Algorithm(dem.clip(area));
  var acc = flowAccumulation(fdir, area, 100);
  var scale = dem.projection().nominalScale();
  var snapped = snapPourPoint(acc, pt, radiusSlider.getValue(), scale);
  Map.clear();
  Map.addLayer(acc.log10(), {min: 0, max: 5, palette: ['ffffff', '0000ff']}, 'Flow Accumulation (log10)');
  Map.addLayer(snapped, {color: 'red'}, 'Snapped Point');
  Map.addLayer(pt, {color: 'yellow'}, 'Original Point');
});

