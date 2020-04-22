#!/usr/bin/env python

from osgeo import ogr
from osgeo import gdal
import osgeo.osr as osr

def main(shapefile):
    print(shapefile)

    #making the shapefile as an object.
    input_shp = ogr.Open(shapefile)

    #getting layer information of shapefile.
    shp_layer = input_shp.GetLayer()

    #pixel_size determines the size of the new raster.
    #pixel_size is proportional to size of shapefile.
    pixel_size = 180

    #get extent values to set size of output raster.
    x_min, x_max, y_min, y_max = shp_layer.GetExtent()

    #calculate size/resolution of the raster.
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)

    #get GeoTiff driver by 
    image_type = 'GTiff'
    driver = gdal.GetDriverByName(image_type)

    #passing the filename, x and y direction resolution, no. of bands, new raster.
    new_raster = driver.Create(output_raster, x_res, y_res, 1, gdal.GDT_Byte)

    #transforms between pixel raster space to projection coordinate space.
    new_raster.SetGeoTransform((x_min, pixel_size, 0, y_min, 0, pixel_size))

    #get required raster band.
    band = new_raster.GetRasterBand(1)

    #assign no data value to empty cells.
    no_data_value = -3.4028234663852886e38
    band.SetNoDataValue(no_data_value)
    band.FlushCache()

    #main conversion method
    gdal.RasterizeLayer(new_raster, [1], shp_layer, burn_values=[255], options=['ALL_TOUCHED=FALSE'])

    #adding a spatial reference
    new_rasterSRS = osr.SpatialReference()
    new_rasterSRS.ImportFromEPSG(4326)
    new_raster.SetProjection(new_rasterSRS.ExportToWkt())
    return gdal.Open(output_raster)
    
    shapefile = '20181109_mountalex_evac_ffdi100d_grid.shp'
    
    main(shapefile)