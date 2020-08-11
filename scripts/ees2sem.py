import os
import gzip
import xml.etree.ElementTree as ET
import numpy as np
from geostack.raster import Raster
from geostack.vector import vector
from geostack.vector import Vector
from geostack import gs_enums
from geostack.utils import get_epsg

#TODO: make these command line arguments
outdir='../data/sem/castlemaine-region/'

popnFile='../scenarios/mount-alexander-shire/castlemaine-region/population-archetypes.xml.gz'
popnEpsg='28355'
popnActivityFilter='home'
popnOutfile= outdir + 'population-archetypes.tif'
popnRasterCellSize='10'

firePhoenix4GridShp = '../data/20181109_mountalex_evac_ffdi100d_grid.shp'
fireOutfile= outdir + '20181109_mountalex_evac_ffdi100d_grid.tif'
fireRasterCellSize='10'

networkFile='../scenarios/mount-alexander-shire/mount_alexander_shire_network_2018.xml.gz'
networkEpsg='28355'
networkOutfilePrefix= outdir + 'mount_alexander_shire_network_2018'


# Parse args
popnEpsg = get_epsg(int(popnEpsg))
popnRasterCellSize = int(popnRasterCellSize)
networkEpsg = get_epsg(int(networkEpsg))

# Create the output dir
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Vectorise the MATSim network
print("Reading MATSim network from " + networkFile)
input = gzip.open(networkFile, 'r')
tree = ET.parse(input)
root = tree.getroot()
print("Parsing network nodes into a vector")
nodes=dict()
vec = vector.Vector()
for elem in root.iter('node'):
    xy=[float(elem.attrib['x']), float(elem.attrib['y'])]
    c = vector.Coordinate.from_list(xy)
    pointIdx = vec.addPoint(c)
    vec.setProperty(pointIdx, "newproperty", "newstr")
    nodes[elem.attrib['id']] = xy # record the mapping of node id to xy coords
vec.setProjectionParameters(networkEpsg)
input.close()
outfile = networkOutfilePrefix + ".nodes.shp"
print("Writing " + outfile)
vec.to_shapefile(outfile, gs_enums.GeometryType.Point)

print("Parsing network links into a vector")
vec = vector.Vector()
for elem in root.iter('link'):
    line = [nodes[elem.attrib['from']], nodes[elem.attrib['to']]]
    # FIXME: 'line' gets somethig like:
    # [[259715.66711096585, 5926741.257714603], [259722.21628692746, 5926750.115245353]] but does not parse below
    #lineIdx = vec.addLineString(line)
    #vec.setProperty(lineIdx, "newproperty", "newstr")
vec.setProjectionParameters(networkEpsg)
input.close()
outfile = networkOutfilePrefix + ".links.shp"
print("Writing " + outfile)
vec.to_shapefile(outfile, gs_enums.GeometryType.LineString)

# Rasterise the fire shapefile
print("Loading Phoenix4 fire grid shapefile " + firePhoenix4GridShp)
vec = Vector.from_shapefile(firePhoenix4GridShp)
print("Rasterising the fire vector")
raster = vec.rasterise(fireRasterCellSize, "output = HOUR_BURNT;")
print("Writing " + fireOutfile)
raster.write(fireOutfile)

# Rasterise the MATSim population home locations
print("Reading MATSim population from " + popnFile)
input = gzip.open(popnFile, 'r')
tree = ET.parse(input)
root = tree.getroot()
print("Parsing home activity coordinates into a vector")
vec = vector.Vector()
for act in root.iter('activity'):
    if(act.attrib['type']==popnActivityFilter):
        c = vector.Coordinate.from_list([float(act.attrib['x']), float(act.attrib['y'])])
        pointIdx = vec.addPoint(c)
        vec.setProperty(pointIdx, "newproperty", "newstr")
vec.setProjectionParameters(popnEpsg)
input.close()
print("Rasterising the home activity locations vector ( CRS size cells)")
bounds = vec.getBounds()
bounds.extend(0.05)
count = Raster(name="count", data_type=np.uint32)
count.init_with_bbox(bounds, popnRasterCellSize)
count.setProjectionParameters(popnEpsg)
count.setAllCellValues(0)
count.rasterise(vec, "atomic_inc();")
print("Writing " + popnOutfile)
count.write(popnOutfile)
