#import libraries
import geostack
from geostack.raster import Raster
from geostack.vector import Vector
from geostack.io import shapefileToVector
from geostack.io import geoJsonToVector
from geostack.io import vectorToShapefile
from geostack.gs_enums import GeometryType
from geostack.core import convert, ProjectionParameters
from geostack.utils import get_epsg
import numpy as np

# Input & output file path for fire file and population file
shape_file = r'../scenarios/mount-alexander-shire/castlemaine-region/20181109_mountalex_evac_ffdi100d_grid.shp'
output_raster_firefie = r'../scenarios/mount-alexander-shire/castlemaine-region/gssem/20181109_mountalex_evac_ffdi100d_grid_geostack.tif'
output_raster_popfile = r'../scenarios/mount-alexander-shire/castlemaine-region/gssem/population.tif'


vec = shapefileToVector(shape_file)
#bounds_ff = vec.getBounds()
#bounds_ff.extend(0.1)

#Creating a raster 
#rasterise_ff = Raster(name="rasterise_ff", data_type=np.uint32)
#rasterise_ff.init_with_bbox(bounds_ff, 0.1)

#rasterise_ff.setAllCellValues(0)

# Setting up the projection
#proj_EPSG4326 = ProjectionParameters.from_proj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#vec.setProjectionParameters(proj_EPSG4326)


raster_ff = vec.rasterise(1.0, "output = x < 505755 ? index : noData_REAL;")
#output_raster = r'./20181109_mountalex_evac_ffdi100d_grid_geostack_ndt.tif'
raster_ff.write(output_raster_firefie);

# Reading Population GeoJson
geofile = "castlemaine-area_pop.json"
#geofile = "Population_file.geojson"
geo = geoJsonToVector(geofile)

#geo.setProjectionParameters(get_epsg(4326))
#vectorToShapefile(geo,  geofile.replace(".json", ".shp"), geom_type=GeometryType.Point)

bounds = geo.getBounds()
bounds.extend(0.1)

#Creating a raster 
rasterise = Raster(name="rasterise", data_type=np.uint32)
rasterise.init_with_bbox(bounds, 0.1)

rasterise.setAllCellValues(0)

# Setting up the projection
proj_EPSG4326 = ProjectionParameters.from_proj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

rasterise.setProjectionParameters(proj_EPSG4326)

rasterise.rasterise(geo, "atomic_inc();")
rasterise.write(output_raster_popfile);
