// Earth Engine flow length (downstream) calculation

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

// Downstream flow length using an iterative D8 approach
function flowLength(flowDir, aoi, iterations) {
  iterations = iterations || 100;
  var flen = ee.Image.constant(0).rename('flowLength');
  var kernel = ee.Kernel.square(1);
  var offsets = [
    [-1, -1], [-1, 0], [-1, 1],
    [ 0, -1],           [ 0, 1],
    [ 1, -1], [ 1, 0],  [ 1, 1]
  ];
  var directions = [32, 64, 128, 16, 1, 8, 4, 2];
  var distances = [Math.SQRT2, 1, Math.SQRT2, 1, 1, Math.SQRT2, 1, Math.SQRT2];
  for (var step = 0; step < iterations; step++) {
    var flNeigh = flen.neighborhoodToBands(kernel);
    var newLen = ee.Image(0);
    for (var j = 0; j < 8; j++) {
      var row = offsets[j][0];
      var col = offsets[j][1];
      var dirBand = 'flowLength_' + row + '_' + col;
      var mask = flowDir.eq(directions[j]);
      var contrib = flNeigh.select(dirBand).add(distances[j]);
      newLen = newLen.where(mask, contrib);
    }
    flen = flen.max(newLen);
  }
  return flen.clip(aoi);
}

// UI elements
var demNames = Object.keys(demList);
var demSelect = ui.Select({items: demNames, value: demNames[0]});
var runBtn = ui.Button('Compute Flow Length');
var panel = ui.Panel([
  ui.Label('Select DEM:'), demSelect,
  ui.Label('Draw an AOI polygon then click the button.'),
  runBtn
]);
panel.style().set({position: 'top-right', width: '300px'});
ui.root.add(panel);

Map.drawingTools().setDrawModes(['polygon']);
Map.drawingTools().draw();

runBtn.onClick(function() {
  var layers = Map.drawingTools().layers();
  if (layers.length() === 0) {
    ui.alert('Please draw an AOI polygon.');
    return;
  }
  var aoi = ee.Feature(layers.get(0)).geometry();
  var dem = demList[demSelect.getValue()];
  var fdir = D8Algorithm(dem);
  var flen = flowLength(fdir, aoi, 100);
  Map.clear();
  Map.addLayer(flen, {min: 0, max: 50000, palette: ['ffffff', 'ff0000']}, 'Flow Length (m)');
  Map.addLayer(aoi, {color: 'blue'}, 'AOI');
});

