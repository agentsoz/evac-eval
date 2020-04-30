#!/usr/bin/env python

from osgeo import ogr, osr
from osgeo import gdal
import os

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
                    
#os.sys('gdalwarp ds output_raster -t_srs "+proj=longlat +ellps=WGS84"')
# Reproject to EPSG: 3111
dest_srs = osr.SpatialReference()
dest_srs.ImportFromEPSG(3111)
dst_wkt = dest_srs.ExportToWkt()
ds = gdal.AutoCreateWarpedVRT(ds,
                              None, # src_wkt : left to default value --> will use the one from source
                              dst_wkt)

                    

ds = None                




