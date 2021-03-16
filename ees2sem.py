import os
import sys
import gzip
import urllib.request
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import numpy as np
from geostack.raster import Raster
from geostack.vector import vector
from geostack.vector import Vector
from geostack import gs_enums
from geostack.utils import get_epsg
from geostack.runner import runScript

def main(outdir):
    # #TODO: make these command line arguments
    # outdir='./data/castlemaine-region/'

    popnFileUrl='https://github.com/agentsoz/ees/raw/dbfab224daaeb02294b5dabb62f55b5f8755b6ce/ees/scenarios/mount-alexander-shire/castlemaine-region-archetypes/population-archetypes.xml.gz'
    popnFile= outdir + 'abm/population-archetypes.xml.gz'
    popnEpsg='28355'
    popnActivityFilter='home'
    popnOutfile= outdir + 'sem/population-archetypes.tif'
    popnRasterCellSize='100'

    firePhoenix4GridShpUrl='https://github.com/agentsoz/ees-data/raw/master/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100d/20181109_mountalex_evac_ffdi100d_grid.shp'
    firePhoenix4GridShp = outdir + 'abm/20181109_mountalex_evac_ffdi100d_grid.shp'
    fireOutfile= outdir + 'sem/20181109_mountalex_evac_ffdi100d_grid.tif'
    fireRasterCellSize='10'

    networkFileUrl='https://github.com/agentsoz/ees/raw/dbfab224daaeb02294b5dabb62f55b5f8755b6ce/ees/scenarios/mount-alexander-shire/mount_alexander_shire_network_2018.xml.gz'
    networkFile= outdir + 'abm/mount_alexander_shire_network_2018.xml.gz'
    networkEpsg='28355'
    networkOutfilePrefix= outdir + 'sem/mount_alexander_shire_network_2018'

    # Parse args
    popnEpsg = get_epsg(int(popnEpsg))
    popnRasterCellSize = int(popnRasterCellSize)
    networkEpsg = get_epsg(int(networkEpsg))
    epsg4326 = get_epsg(4326)

    # Create the output dir
    os.makedirs(outdir + "/abm", exist_ok=True)
    os.makedirs(outdir + "/sem", exist_ok=True)

    # Download all files
    if not os.path.exists(popnFile):
        print("Downloading " + popnFileUrl + " to " + popnFile)
        urllib.request.urlretrieve(popnFileUrl, popnFile)
    if not os.path.exists(networkFile):
        print("Downloading " + networkFileUrl + " to " + networkFile)
        urllib.request.urlretrieve(networkFileUrl, networkFile)
    for ext in ['shp', 'prj', 'dbf', 'cpg', 'shx']:
        src = firePhoenix4GridShpUrl.replace('shp', ext)
        dst = firePhoenix4GridShp.replace('shp', ext)
        if not os.path.exists(dst):
            print("Downloading " + src + " to " + dst)
            urllib.request.urlretrieve(src, dst)

    # Vectorise the MATSim network
    print("Reading MATSim network from " + networkFile)
    input = gzip.open(networkFile, 'r')
    tree = ET.parse(input)
    root = tree.getroot()
    print("Parsing network nodes into a vector")
    nodes=dict()
    vec = vector.Vector()
    for elem in root.iter('node'):
        c = vector.Coordinate(float(elem.attrib['x']), float(elem.attrib['y']))
        pointIdx = vec.addPoint(c)
        vec.setProperty(pointIdx, "newproperty", "newstr")
        nodes[elem.attrib['id']] = c[:2] # record the mapping of node id to xy coords
    vec.setProjectionParameters(networkEpsg)
    input.close()
    outfile = networkOutfilePrefix + ".nodes.shp"
    print("Writing " + outfile)
    vec.to_shapefile(outfile, gs_enums.GeometryType.Point)
    outfile = networkOutfilePrefix + ".nodes.geojson"
    print("Writing " + outfile)
    vec.convert(epsg4326)
    vec.to_geojson(outfile)

    print("Parsing network links into a vector")
    # vec = vector.Vector()
    for elem in root.iter('link'):
        line = [nodes[elem.attrib['from']], nodes[elem.attrib['to']]]
        lineIdx = vec.addLineString(line)
        vec.setProperty(lineIdx, "matsim_linkID", elem.attrib['id'])
        vec.setProperty(lineIdx, "diameter", 1.0)
    vec.setProjectionParameters(networkEpsg)
    input.close()
    outfile = networkOutfilePrefix + ".links.shp"
    print("Writing " + outfile)
    vec.to_shapefile(outfile, gs_enums.GeometryType.LineString)
    outfile = networkOutfilePrefix + ".links.geojson"
    print("Writing " + outfile)
    vec.convert(epsg4326)
    vec.to_geojson(outfile)

    # Rasterise the fire shapefile
    print("Loading Phoenix4 fire grid shapefile " + firePhoenix4GridShp)
    vec = Vector.from_shapefile(firePhoenix4GridShp)
    print("Rasterising the fire vector")
    raster = vec.rasterise(fireRasterCellSize, "output = HOUR_BURNT;")
    # convert hours to seconds
    runScript("rasterised = rasterised * 3600.0;", [raster])
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
            c = vector.Coordinate(float(act.attrib['x']), float(act.attrib['y']))
            pointIdx = vec.addPoint(c)
            vec.setProperty(pointIdx, "newproperty", "newstr")
    vec.setProjectionParameters(popnEpsg)
    input.close()
    print("Rasterising the home activity locations vector")
    bounds = vec.getBounds()
    bounds.extend(popnRasterCellSize)
    count = Raster(name="count", data_type=np.uint32)
    count.init_with_bbox(bounds, popnRasterCellSize)
    count.setProjectionParameters(popnEpsg)
    count.setAllCellValues(0)
    count.rasterise(vec, "atomic_inc();")
    print("Writing " + popnOutfile)
    count.write(popnOutfile)

if __name__ == "__main__":
    parser = ArgumentParser(description="script to convert vector and raster files using geostack")
    parser.add_argument("--outdir", dest="outdir", type=str, help="path to the output directory")
    args = parser.parse_args()

    if args.outdir is None:
        parser.print_help()
    else:
        main(args.outdir)
