#import libraries
import geostack
from geostack.raster import Raster
from geostack.vector import Vector
from geostack.utils import get_epsg

# TODO: Make this a script parameter
shp = '../data/20181109_mountalex_evac_ffdi100d_grid.shp'

print("Loading shapefile " + shp)
vector = Vector.from_shapefile(shp)
print("Setting CRS to EPSG:28355")
vector = vector.convert(get_epsg(28355))
print("Rasterising")
raster = vector.rasterise(10, "output = HOUR_BURNT;")
outfile = shp + ".tif"
print("Writing " + shp.replace("shp", "tif"))
raster.write(outfile)
