// Run in the signed-in Earth Engine Code Editor. Launch the resulting Drive task.
// The output preserves the source MODIS pixels and carries raw DN and QC.
var name = 'singapore_mod11a1_061_2025_native_daily';
var region = ee.Geometry.Rectangle([103.54, 1.11, 104.16, 1.54], null, false);
var inputs = ee.ImageCollection('MODIS/061/MOD11A1')
  .filterDate('2025-01-01', '2026-01-01')
  .select(['LST_Day_1km', 'QC_Day'])
  .sort('system:time_start');
var dated = inputs.map(function (image) {
  return image.set('date_utc', image.date().format('yyyyMMdd', 'UTC'));
});

dated.aggregate_array('date_utc').evaluate(function (dates, error) {
  if (error || !dates) {
    print('Could not read the Earth Engine calendar', error);
    return;
  }
  var expected = [];
  for (var i = 0; i < 365; i++) {
    expected.push(new Date(Date.UTC(2025, 0, 1 + i)).toISOString().slice(0, 10).replace(/-/g, ''));
  }
  if (!dates.length || dates.some(function (date, index) {
    return expected.indexOf(date) < 0 || (index > 0 && dates[index - 1] >= date);
  })) {
    print('Export stopped because dates are absent, repeated, unordered, or outside 2025', dates);
    return;
  }
  var missing = expected.filter(function (date) { return dates.indexOf(date) < 0; });
  var first = ee.Image(dated.first());
  var projection = first.select('LST_Day_1km').projection().getInfo();
  var exportImages = dated.map(function (image) {
    var date = ee.String(image.get('date_utc'));
    var dn = image.select('LST_Day_1km').toInt32().unmask(-9999, false).rename('lst_dn');
    // The product uses QC=0 both as a storage fill code and as a best-quality flag.
    // Keep it as zero when paired with nonzero LST. The server audits overlapping mirror days.
    var qc = image.select('QC_Day').toInt32().unmask(0, false)
      .where(dn.lte(0), -9999).rename('qc');
    // A date band makes the order verifiable even if the GeoTIFF omits band names.
    return ee.Image.constant(ee.Number.parse(date)).rename('date')
      .addBands(dn).addBands(qc)
      .toInt32()
      .set('system:index', date);
  }).sort('system:index');
  Export.image.toDrive({
    image: exportImages.toBands(),
    description: name,
    fileNamePrefix: name,
    folder: 'SG_LST_2025',
    region: region,
    crs: projection.crs,
    crsTransform: projection.transform,
    maxPixels: 10000000,
    fileFormat: 'GeoTIFF',
    formatOptions: {noData: -9999}
  });
  print('Available dates, missing dates, and source projection', dates.length, missing, projection);
  print('Open Tasks and run ' + name + '.');
});
