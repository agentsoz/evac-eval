import gzip
import xml.etree.ElementTree as ET
import numpy as np
from geostack.raster import Raster
from geostack.vector import vector
from geostack.utils import get_epsg

#TODO: make these command line arguments
popnfile='../scenarios/mount-alexander-shire/castlemaine-region/population-archetypes.xml.gz'
popnepsg='28355'
activityTypeToFilter='home'
outfile='../data/population-archetypes.tif'

# Get the EPSG
epsg = get_epsg(int(popnepsg))

# Vectorise
print("Reading MATSim population from " + popnfile)
input = gzip.open(popnfile, 'r')
tree = ET.parse(input)
root = tree.getroot()
print("Parsing home activity coordinates into a vector")
vec = vector.Vector()
for act in root.iter('activity'):
    if(act.attrib['type']==activityTypeToFilter):
        c = vector.Coordinate.from_list([float(act.attrib['x']), float(act.attrib['y'])])
        pointIdx = vec.addPoint(c)
        vec.setProperty(pointIdx, "newproperty", "newstr")
vec.setProjectionParameters(epsg)
input.close()

# Rasterise
print("Rasterising")
bounds = vec.getBounds()
bounds.extend(0.05)
count = Raster(name="count", data_type=np.uint32)
count.init_with_bbox(bounds, 50)
count.setProjectionParameters(epsg)
count.setAllCellValues(0)
count.rasterise(vec, "atomic_inc();")

# Write it out
print("Writing " + outfile)
count.write(outfile)
