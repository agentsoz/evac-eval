import subprocess
import csv
from geostack.io import geoJsonToVector
from geostack.raster import Raster
from geostack.core import ProjectionParameters
import numpy as np

# Create the input CSV using this bash script or some other means
#subprocess.run(["./createMATSimPlansHomeCsv.sh"])

filename = '../data/population-archetypes.home.csv'
geojson = '{"features": ['
with open(filename, "r") as f:
    epsg = next(f).rstrip('\n') # header row is the epsg for the coords
    csv_reader = csv.reader(f, delimiter=',') # comma separated X,Y
    for r in csv_reader:
        s = '\n{"geometry": {"coordinates": ['+r[0]+', '+r[1]+'], "type": "Point"}, "properties": {}, "type": "Feature"},'
        geojson += s
geojson = geojson[:-1] # remove the last comma
geojson += '\n], "type": "FeatureCollection"}'

vec = geoJsonToVector(geojson)

# Get bounds of vector
bounds = vec.getBounds()
bounds.extend(0.05)
# Create raster
count = Raster(name="count", data_type=np.uint32)
count.init_with_bbox(bounds, 0.001)
# Set raster projection
proj_EPSG4326 = ProjectionParameters.from_proj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
count.setProjectionParameters(proj_EPSG4326)
# Set raster to zero
count.setAllCellValues(0)
# Rasterise vector
count.rasterise(vec, "atomic_inc();")
# Write it out
outfile=filename.replace("csv", "tif")
print("Writing " + outfile)
count.write(outfile)
