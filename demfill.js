// Earth Engine DEM depression filling demonstration

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

// UI elements
var demNames = Object.keys(demList);
var demSelect = ui.Select({items: demNames, value: demNames[0]});
var runBtn = ui.Button('Fill DEM Depressions');
var panel = ui.Panel([
  ui.Label('Select DEM:'), demSelect,
  ui.Label('Click the button to compute a pit-filled DEM.'),
  runBtn
]);
panel.style().set({position: 'top-right'});
ui.root.add(panel);

runBtn.onClick(function() {
  var dem = demList[demSelect.getValue()];
  // Fill depressions to enforce drainage
  var filled = ee.Algorithms.Hydrology.fill(dem);
  var depth = filled.subtract(dem);
  Map.clear();
  Map.addLayer(dem, {min: 0, max: 3000, palette: ['ffffff', '000000']}, 'Original DEM');
  Map.addLayer(filled, {min: 0, max: 3000, palette: ['ffffff', '0000ff']}, 'Filled DEM');
  Map.addLayer(depth.selfMask(), {min: 0, max: 50, palette: ['ffff00', 'ff0000']}, 'Fill Depth');
});
