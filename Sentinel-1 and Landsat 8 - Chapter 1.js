//Code: Sentinel-1 and Landsat 8 Vegetation Indices
//Inputs/Arguments: Polygons of study areas; date (month) of analysis
//Author: Erli Pinto dos Santos
//        Agronomy Engineer
//        PhD candidate in Agricultural Engineering
//        Universidade Federal de Viçosa (Viçosa, Minas Gerais, Brazil)
//        E-mail: erlispinto@outlook.com or erli.santos@ufv.br

// Code associated with the paper at International Journal of Remote Sensing:
// "Vegetation cover monitoring in tropical regions using SAR-C dual-polarization
// index: seasonal and spatial influences".
// Check it out at: < https://doi.org/10.1080/01431161.2021.1959955 >. If
// the code was useful, please, cite this paper.

//Last update: March 27, 2021

/* ANALYSIS MONTH: */
// The input vars 'start' and 'finish' controls the imagery acquisition.
// The date must initiate in the first day of month (start) in analysis
// and finish at the second day of the subsequent month (finish),
// totalizing 30 days. The exception is to February month, that the subsequent
// date must finish in March 4th, to complete 30 days.
var start = ee.Date('2015-10-01');
var finish = ee.Date('2015-11-02');

/** STUDY AREA POLYGONS: */
// The geometry roi is especially needs to when apply the algorithm over entire Doce river Watershed,
// to avoid select excess imagery
var roi = ee.FeatureCollection('users/erlisantos/Poligono_buffer_BH_RioDoce'
                               //'users/erlisantos/Poligono_BH_77992'
                               //'users/erlisantos/Poligono_BH_78623'
                               //'users/erlisantos/Poligono_BH_779926'
                               //'users/erlisantos/Poligono_BH_786232'
                               ).geometry();
// Note that from roi to aoi only the first geometry changes (to 'users/erlisantos/Poligono_BH_UPH_RioDoce'):
var aoi = ee.FeatureCollection('users/erlisantos/Poligono_BH_UPH_RioDoce'
                               //'users/erlisantos/Poligono_BH_77992'
                               //'users/erlisantos/Poligono_BH_78623'
                               //'users/erlisantos/Poligono_BH_779926'
                               //'users/erlisantos/Poligono_BH_786232'
                               ).geometry();

/** LANDSAT 8 DATA PROCESSING: */
// The goal of this algorithm is take Landsat 8 Surface Reflectance scenes over the AOI's and compute
// both NDVI and EVI spectral indices for all area and free of cloud cover and cloud shadown. Thus,
// the output of this code is a table with pixel samples (NDVI) collected trought the selected AOI.
// Note that all pixels with NDVI < 0.3 are masked, because this algorithm is designed to vegetation analysis.
// Enjoy it! 

  /** SRTM IMAGE IN THE SCRIPT */
  // Is known that the topographical features of an area infer in the backscattering behavior of radar data
  // collected by SAR (Synthetic Aperture Radar) sensors. To study the topographical influences about
  // Sentinel-1 SAR data, we include SRTM information in this script, with the goal avoid samples radar
  // pixels under effect of shadow or layover mechanims.

    /** LOAD THE SRTM IMAGE: */
    var srtm = ee.Image('USGS/SRTMGL1_003')
                 .clip(aoi);

    /** APPLY AN ALGORITHM TO COMPUTE TERRAIN SLOPE: */
    var slope = ee.Terrain.slope(srtm);

    /** MASK SLOPE IMAGE OMMITING SLOPE PIXELS > 10°: */
    var Slope10 = slope.select('slope').lte(10);
    var slopeNew = slope.updateMask(Slope10);

    /** FUNCTION TO ADD SLOPE BAND TO EACH LANDSAT 8 SCENE: */
    function Combine(image){
      return image.addBands(slopeNew.select('slope'));
    }

  /** FUNCTION TO MASK INVALID BANDS VALUES */
  // This function will mask values out of valid range to the band B2, B4 and B5 necessary to this approach
  // according to Valid Range, noticed in Landsat 8 Surface Reflectance Code (LaSRC),
  // section Surface Reflectance Specifications. Available for consultation here:
  //      https://www.usgs.gov/media/files/land-surface-reflectance-code-lasrc-product-guide
  function maskInvalid(image){
    var validPixels = image.gte(0.0);
    return image.updateMask(validPixels);
  }

  /** FUNCTION TO MASK CLOUDS BASED ON THE pixel_qa BAND OF LANDSAT 8 SR DATA */
  function maskL8sr(image) {
    // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    // Get the pixel QA band.
    var qa = image.select('pixel_qa');
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    return image.updateMask(mask);
  }

  /** FUNCTION TO COMPUTE NDVI AND EVI */
  function VI(image){
    var NDVI = image.expression(
      '(B5 - B4) / (B5 + B4)', {
        'B5': image.select('B5'),
        'B4': image.select('B4')
      });
    var EVI = image.expression(
      '(G * (B5 - B4))/(L+B5+(B4*C1)+(B2*C2))', {
        'B5': image.select('B5'),
        'B4': image.select('B4'),
        'B2': image.select('B2'),
        'L': ee.Number(1.0),
        'G': ee.Number(2.5),
        'C1': ee.Number(6.0),
        'C2': ee.Number(7.5)
      });
    return image.addBands([NDVI.rename("NDVI"), EVI.rename("EVI")]);
  }

  /** FUNCTIONS TO MASK NDVI < 0.3 and EVI < 0.15 (Non vegetated areas): */
  function maskNDVI(image){
    var NDVIlessThan = image.select('NDVI').gte(0.3);
    return image.updateMask(NDVIlessThan);
  }

  /** FUNCTION TO SUBSET ALL IMAGES FROM A COLLECTION: */
  function subset(image) {
    return image.clip(aoi);
  }

  /** CREATING THE VARIABLE AS COLLECTION TO STORE THE FILTERED AND PROCESSED SCENES: */
  // Defining the range period:
  var landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                   .filterBounds(aoi)
                   .filterDate(start, finish)
                   .map(maskInvalid)
                   .map(maskL8sr)
                   .map(VI)
                   .map(subset)
                   .map(maskNDVI)
                   .map(Combine);
  //print(landsat8, 'Landsat 8 Surface Reflectance without clouds');

  /** MOSAICKING LANDSAT 8 SCENES BY DATE (ONE MOSAIC BY MONTH) */
  // Difference in days between start and finish
  var diff = finish.difference(start, 'month');
  //print(diff, "Difference of months");

  // Make a list of all dates
  var range = ee.List.sequence(0, diff.subtract(1)).map(function(day){return start.advance(day,'month')});
  //print(range, "List of all months");

  // Funtion for iteraton over the range of dates
  var day_mosaics = function(date, newlist) {
    // Cast
    date = ee.Date(date);
    newlist = ee.List(newlist);
  
    // Filter collection between date and the next day
    var filtered = landsat8.filterDate(date, date.advance(1,'month'));

    // get date as YEARMONTHDAY. For example, for January 8th 2010
    // would be: 20100108
    var date_formatted = ee.Number.parse(date.format('YYYYMMdd'));
  
    // Make the mosaic
    var image = ee.Image(filtered.mosaic()).set('Mosaic', date_formatted);

    // Add the mosaic to a list only if the collection has images
    return ee.List(ee.Algorithms.If(filtered.size(), newlist.add(image), newlist));
  };

  // Iterate over the range to make a new list, and then cast the list to an imagecollection
  var newcoll8 = ee.ImageCollection(ee.List(range.iterate(day_mosaics, ee.List([]))));
  //print(newcoll8, "Collection with Landsat 8 mosaic images");

/** SENTINEL-1 DATA PROCESSING: */
// The goal of this algorithm is process SAR images from Sentinel-1. The code does:

// 1) Creates a function to apply an speckle filter: in this case by Focal Median method;

// 2) Creates a function that will, based on band histogram (from VV and VH FLOAT bands),
// remove an "heavy-tail" from their histograms: this part of distribution that is understood as
// outliers;

// 3) Applying all functions above, the code creates a variable containing a image
// collection for each band. The result are the bands Sigma0_VV and Sigma0_VH;

// 5) Creates and apply functions to compute all the SAR vegetation indices of interest:
// DPSVI original (Periasamy, 2018) and DPSVI modified;

// 6) Finally, the code create mosaics for the area, combining images from the data period
// selectioned;

  /** Function to apply speckle noise removal filter by Focal Median: */
  // The window (the kernel) used to apply the filter is in pixel number, that is odd number
  function MedianVV(image){
    return image.addBands(image.focal_median({radius: 5,
                                              kernelType: 'square',
                                              units: 'pixels'})
                               .rename(['Sigma0_VV']));
  }
  function MedianVH(image){
    return image.addBands(image.focal_median({radius: 5,
                                              kernelType: 'square',
                                              units: 'pixels'})
                               .rename(['Sigma0_VH']));
  }

  /** Function to mask outliers in VV band based on histogram analysis */
  function maskVV(image){
    var quantil90 = image.reduceRegion({reducer:
                            ee.Reducer.percentile({percentiles: ee.List([95]),
                                                   outputNames: ee.List(['p95'])}),
                                        scale: 10,
                                        bestEffort: true});
    var Outliers = image.select('VV').lte(ee.Number(quantil90.get('VV')));
    return image.updateMask(Outliers);}
  function maskVH(image){
    var quantil90 = image.reduceRegion({reducer:
                            ee.Reducer.percentile({percentiles: ee.List([95]),
                                                   outputNames: ee.List(['p95'])}),
                                        scale: 10,
                                        bestEffort: true});
    var Outliers = image.select('VH').lte(ee.Number(quantil90.get('VH')));
    return image.updateMask(Outliers);}
  
  /** Vars that will receive and pre-process the Sentinel-1 SAR choiced images: */
  var sentinel1VV = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
                      .filter(ee.Filter.eq('instrumentMode', 'IW'))
                      .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                      //.filter(ee.Filter.eq('relativeOrbitNumber_start', 82)) //Filter necessary to basin 779926
                      .filter(ee.Filter.eq('resolution_meters', 10))
                      .filterBounds(aoi) //roi: Necessary to Doce River Basin, to avoid select excess images
                      .filterDate(start, finish)
                      .select('VV')
                      .map(maskVV)
                      .map(MedianVV)
                      .select('Sigma0_VV');
                  
  // VH band
  var sentinel1VH = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
                      .filter(ee.Filter.eq('instrumentMode', 'IW'))
                      .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                      //.filter(ee.Filter.eq('relativeOrbitNumber_start', 82)) //Filter necessary to basin 779926
                      .filter(ee.Filter.eq('resolution_meters', 10))
                      .filterBounds(aoi) //roi: Necessary to Doce River Basin, to avoid select excess images
                      .filterDate(start, finish)
                      .select('VH')
                      .map(maskVH)
                      .map(MedianVH)
                      .select('Sigma0_VH');

  /** Combining both image collections (Sigma0_VV and Sigma0_VH): */
  var sentinel1 = sentinel1VV.combine(sentinel1VH);
  //print(sentinel1, "Sentinel-1 filtered images");

  /** Computing vegetation index over each image in collection: */
  // The outputs of this function are:
      // The original Dual Polarization SAR Vegetation Index (DPSVI; Perisamy, 2018);
      // and the modified DPSVI (proposed here):
  function indices(image){
    var max = image.reduceRegion({reducer: ee.Reducer.max(), scale: 10,
                                  geometry: aoi, bestEffort: true});
    var DPSVIoriginal = image.expression(
      '(((VVmax - VV)+VH)/1.414213562) * ((VV+VH)/VV) * VH', {
        'VH': image.select('Sigma0_VH'),
        'VV': image.select('Sigma0_VV'),
        'VVmax': ee.Number(max.get('Sigma0_VV'))
      });
    var DPSVIm = image.expression(
      '(VV*VV+VV*VH)/1.414213562',{
        'VH': image.select('Sigma0_VH'),
        'VV': image.select('Sigma0_VV')
      });
    
    return image.addBands([DPSVIm.rename("DPSVIm"),
                           DPSVIoriginal.rename("DPSVIo")]);}

  // The follow function apply the DPSVI normalization procedure:
  function DPSVInorm(image){
    var max = image.reduceRegion({reducer: ee.Reducer.max(),
                                  scale: 10,
                                  geometry: aoi,
                                  bestEffort: true});
    var min = image.reduceRegion({reducer: ee.Reducer.min(),
                                  scale: 10,
                                  geometry: aoi,
                                  bestEffort: true});
    var DPSVI = image.expression(
      '(DPSVI - DPSVImin) /(DPSVImax - DPSVImin)',{
        'DPSVI': image.select('DPSVIm'),
        'DPSVImax': ee.Number(max.get('DPSVIm')),
        'DPSVImin': ee.Number(min.get('DPSVIm'))
      });
    return image.addBands(DPSVI.rename('DPSVI'));
  }

  // Now the functions of this section are applied over Image Collection:
  var S1 = sentinel1.map(indices);
  var S1new = S1.map(DPSVInorm)
                .select(['Sigma0_VV', 'Sigma0_VH', 'DPSVIo', 'DPSVI'])
                .map(subset);
  //print(S1new, "Pos processing images S1");
  
  /** Mosaicking Sentinel-1 scenes by date */
  // Difference in days between start and finish
  var diff = finish.difference(start, 'month');
  //print(diff, "Difference of months");

  // Make a list of all dates
  var range = ee.List.sequence(0, diff.subtract(1)).map(function(day){return start.advance(day,'month')});
  print(range, "List of all months");

  // Funtion for iteraton over the range of dates
  var day_mosaics = function(date, newlist) {
    // Cast
    date = ee.Date(date);
    newlist = ee.List(newlist);
  
  // Filter collection between date and the next day
  var filtered = S1new.filterDate(date, date.advance(1,'month'));
 
  // get date as YEARMONTHDAY. For example, for January 8th 2010
  // would be: 20100108
  var date_formatted = ee.Number.parse(date.format('YYYYMMdd'));
  
  // Make the mosaic
  var image = ee.Image(filtered.mosaic()).set('Mosaic', date_formatted);

  // Add the mosaic to a list only if the collection has images
  return ee.List(ee.Algorithms.If(filtered.size(), newlist.add(image), newlist));
};

  // Iterate over the range to make a new list, and then cast the list to an imagecollection
  var newcols1 = ee.ImageCollection(ee.List(range.iterate(day_mosaics, ee.List([]))));
  //print(newcols1, "Collection with Sentinel-1 mosaic images");

/** COMBINING LANDSAT 8 AND SENTINEL-1 MOSAICS IN ONE IMAGE COLLECTION: */
var images = newcoll8.select(['NDVI','EVI','slope'])
                     .combine(newcols1);
//print(images, "Image collection containing L8 and S1 data");

/** VISUALIZING ON THE MAP THE AOI AND SCENES PROCESSED: */
Map.centerObject(aoi);
//Map.addLayer(aoi);
//var NDVIvisParams = {min: 1, max: 33};
//Map.addLayer(images.first(), {bands: ['EVI'], visParams: NDVIvisParams}, "EVI");
//Map.addLayer(images.first(), {bands: ['DPSVI'], visParams: NDVIvisParams}, "DPSVI");

/** GETTING RANDOMLY SAMPLES: */
// Sample mosaic images
//var listOfImages = images.toList(images.size());
//print(listOfImages, "Lista de imagens selecionadas e processadas.");

//var image = ee.Image(listOfImages.get(0));
var image = images.first();
//print(image, "Product in analysis");

var sample = image.sample({region: aoi,
                           scale: 30,
                           tileScale: 16,
                           numPixels: 3000 //Note that numPixels is the approximated number of samples to collect
});
//print(sample, "Amostras coletadas");

// Export sample tables obtained from Sentinel-1 and Landsat-8 mosaics at same position:
Export.table.toDrive({collection: sample,
                      folder: 'Chapter1_Indices_baciaDoce',
                      description: 'DPSVI-NDVI-EVI_2015-10'});

// Visualize the behavior of spectral and microwave indices:
var chart = ui.Chart.feature.byFeature(sample, 'NDVI', ['DPSVI'])
              .setChartType('ScatterChart')
              .setOptions({pointSize: 2,
                           pointColor: 'red',
                           title: "Bacia Doce",
                           titleY: 'DPSVI (dl)',
                           titleX: 'NDVI (dl)',
                           trendlines:{0: {
                                          type: 'linear',
                                          color: 'green',
                                          showR2: true,
                                          visibleInLegend: true}}});
print(chart);
var chart1 = ui.Chart.feature.byFeature(sample, 'NDVI', ['DPSVIo'])
              .setChartType('ScatterChart')
              .setOptions({pointSize: 2,
                           pointColor: 'red',
                           title: "Bacia Doce",
                           titleY: 'DPSVI',
                           titleX: 'NDVI (dl)',
                           trendlines:{0: {
                                          type: 'linear',
                                          color: 'green',
                                          showR2: true,
                                          visibleInLegend: true}}});
//print(chart1);

/** ANALYSING THE IMAGE USING A BAND SCATTERPLOT */
//var result = image.reduceRegion({reducer: ee.Reducer.toList(),
//                                 geometry: aoi,
//                                 bestEffort: true,
//                                 scale: 1000});
// Convert the band data to plot on the y-axis to arrays.
//var yValues = ee.Array(result.get('DPSVI'));
// The band data to plot on the x-axis is a List.
//var xValues = result.get('EVI');
// Make a band correlation chart.
//var chart2 = ui.Chart.array.values(yValues, 0, xValues)
//               .setOptions({title: '2015-10-01',
//                            hAxis: {'title': 'EVI (dl)'},
//                            vAxis: {'title': 'DPSVI (dl)'},
//                            pointSize: 3});
//Print the chart.
//print(chart2, "Band scatter plot");