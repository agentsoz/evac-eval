#!/usr/bin/env python

from osgeo import ogr
from osgeo import gdal

shape_file = r'../scenarios/mount-alexander-shire/castlemaine-region/20181109_mountalex_evac_ffdi100d_grid.shp'
output_raster = r'../scenarios/mount-alexander-shire/castlemaine-region/gssem/20181109_mountalex_evac_ffdi100d_grid.tif'

input_shp = ogr.Open(shape_file)
shp_layer = input_shp.GetLayer()

pixel_size = 180
xmin, xmax, ymin, ymax = shp_layer.GetExtent()

#get GeoTiff driver by 
image_type = 'GTiff'
driver = gdal.GetDriverByName(image_type)


ds = gdal.Rasterize(output_raster, shape_file, xRes=pixel_size, yRes=pixel_size, 
                    burnValues=255, outputBounds=[xmin, ymin, xmax, ymax], 
                    outputType=gdal.GDT_Byte)
ds = None


