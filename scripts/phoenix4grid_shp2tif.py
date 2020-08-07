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

# TODO: Make this a script parameter
pheonix4GridShapeFile = '../data/20181109_mountalex_evac_ffdi100d_grid.shp'

print("Loading shapefile as vector " + pheonix4GridShapeFile)
vector = shapefileToVector(pheonix4GridShapeFile)

#print("Setting CRS to EPSG4326")
#proj_EPSG4326 = ProjectionParameters.from_proj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#vector.setProjectionParameters(proj_EPSG4326)

print("Rasterising")
raster = vector.rasterise(100.0, "output = HOUR_BURNT;")

outfile = pheonix4GridShapeFile + ".tif"
print("Writing " + outfile)
raster.write(outfile);
