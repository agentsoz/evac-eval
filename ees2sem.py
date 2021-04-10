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

# Run this program inside directory evac-eval/.

def main(scenario, outdir):
    # TODO: make these parameters into command-line arguments
    # outdir = './data/castlemaine-region/'

    if scenario in ('Mount-Alexander-Shire', 'Loddon-Mallee-Northern-Cluster-Shires'):  # TODO: is it correct to use these (population and fire) files for the Loddon-Mallee scenario?
        popnFileUrl = 'https://github.com/agentsoz/ees/raw/dbfab224daaeb02294b5dabb62f55b5f8755b6ce/ees/scenarios/mount-alexander-shire/castlemaine-region-archetypes/population-archetypes.xml.gz'
        popnFile = outdir + 'abm/population-archetypes.xml.gz'
        popnEpsg = '28355'
        popnActivityFilter = 'home'
        popnOutfile = outdir + 'sem/population-archetypes.tif'
        popnRasterCellSize = '100'
    
        firePhoenix4GridShpUrl = 'https://github.com/agentsoz/ees-data/raw/master/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100d/20181109_mountalex_evac_ffdi100d_grid.shp'
        firePhoenix4GridShp = outdir + 'abm/20181109_mountalex_evac_ffdi100d_grid.shp'

    else:
        popnFileUrl = 'Not applicable'
        popnFile = 'Not applicable'
        popnEpsg = 'Not applicable'
        popnActivityFilter = 'Not applicable'
        popnOutfile = 'Not applicable'
        popnRasterCellSize = 'Not applicable'
    
        firePhoenix4GridShpUrl = 'Not applicable'
        firePhoenix4GridShp = 'Not applicable'

    if scenario == 'Mount-Alexander-Shire':  # run with, e.g., python ees2sem.py --scenario Mount-Alexander-Shire --outdir modelinputsMtAlexander/
        fireOutfile = outdir + 'sem/20181109_mountalex_evac_ffdi100d_grid.tif'
        fireRasterCellSize = '10'
    
        networkFileUrl = 'https://github.com/agentsoz/ees/raw/dbfab224daaeb02294b5dabb62f55b5f8755b6ce/ees/scenarios/mount-alexander-shire/mount_alexander_shire_network_2018.xml.gz'
        networkFile = outdir + 'abm/mount_alexander_shire_network_2018.xml.gz'
        networkEpsg = '28355'
        networkOutfilePrefix = outdir + 'sem/mount_alexander_shire_network_2018'

    elif scenario == 'Loddon-Mallee-Northern-Cluster-Shires':  # run with, e.g., python ees2sem.py --scenario loddon-mallee-northern-cluster-shires --outdir modelinputsLoddonMallee/
        fireOutfile = outdir + 'sem/20181109_loddonMalleeNorthernCluster_evac_ffdi100d_grid.tif'
        fireRasterCellSize = '10'
    
        networkFileUrl = 'https://github.com/agentsoz/ees/raw/dbfab224daaeb02294b5dabb62f55b5f8755b6ce/ees/scenarios/loddon-mallee-northern-cluster-shires/loddon_mallee_northern_cluster_shires_network.xml.gz'
        networkFile = outdir + 'abm/loddon_mallee_northern_cluster_shires_network.xml.gz'
        networkEpsg = '28355'
        networkOutfilePrefix = outdir + 'sem/loddon_mallee_northern_cluster_shires_network'

    elif scenario == 'cmr_1s1d1r1k':  # run with, e.g., python ees2sem.py --scenario cmr_1s1d1r1k --outdir modelinputsCmr_1s1d1r1k/
        fireOutfile = 'Not applicable'  # for this simple scenario the location 'fireOutfile' is unspecified, because no fire-layer is necessary for this scenario
        fireRasterCellSize = 'Not applicable'

        networkFileUrl = 'https://github.com/agentsoz/evac-eval/blob/b57d9f72fc72d470020bb84c30def116f7637f3b/scenarios/mount-alexander-shire/cmr_1s1d1r1k/cmr_1s1d1r_network.xml.gz'  # ST: couldn't successfully download from this URL, getting a file of format HTML rather than the desired GZIP, so I cheated by downloading from this URL manually, i.e. in a browser
        networkFile = outdir + 'abm/cmr_1s1d1r_network.xml.gz'
        networkEpsg = '28355'
        networkOutfilePrefix = outdir + 'sem/cmr_1s1d1r_network'

    else:
      print(f"Scenario '{scenario}' is unrecognised; terminating program.")
      sys.exit()

    # Parse args
    if popnEpsg != 'Not applicable':
        popnEpsg = get_epsg(int(popnEpsg))
        popnRasterCellSize = int(popnRasterCellSize)
    networkEpsg = get_epsg(int(networkEpsg))
    epsg4326 = get_epsg(4326)

    # Create the output dir
    os.makedirs(outdir + "/abm", exist_ok=True)
    os.makedirs(outdir + "/sem", exist_ok=True)

    # Download all files
    if popnFile != 'Not applicable':
        if not os.path.exists(popnFile):
            print("Downloading " + popnFileUrl + " to " + popnFile)
            urllib.request.urlretrieve(popnFileUrl, popnFile)
    if not os.path.exists(networkFile):
        print("Downloading " + networkFileUrl + " to " + networkFile)
        urllib.request.urlretrieve(networkFileUrl, networkFile)
    if firePhoenix4GridShpUrl != 'Not applicable':
        for ext in ['shp', 'prj', 'dbf', 'cpg', 'shx']:
            src = firePhoenix4GridShpUrl.replace('shp', ext)
            dst = firePhoenix4GridShp.replace('shp', ext)
            if not os.path.exists(dst):
                print("Downloading " + src + " to " + dst)
                urllib.request.urlretrieve(src, dst)

    # Vectorise the MATSim network
    print("Reading MATSim network from " + networkFile)
    input = gzip.open(networkFile, 'r')
#    try:
#        input = gzip.open(networkFile, 'r')
#    except:
#        input = open(networkFile, 'r')
    tree = ET.parse(input)
    root = tree.getroot()
    print("Parsing network nodes into a vector")
    nodes = dict()
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
        vec.setProperty(lineIdx, "diameter", 1.0)  # TODO: "diameter" applies only to the monotonic-flow version of SEM, but for that version it shouldn't always be 1.0
        vec.setProperty(lineIdx, "matsim_linkID", elem.attrib['id'])
        vec.setProperty(lineIdx, "capacity", elem.attrib['capacity'])  # for version of SEM with non-monotonic flow, need capacity, free speed, number of lanes, and length of each link
        vec.setProperty(lineIdx, "freespeed", elem.attrib['freespeed'])
        vec.setProperty(lineIdx, "permlanes", elem.attrib['permlanes'])
        vec.setProperty(lineIdx, "diameter", elem.attrib['permlanes'])
        vec.setProperty(lineIdx, "length", elem.attrib['length'])
        vec.setProperty(lineIdx, "oneway", elem.attrib['oneway'])  # might not need 'oneway' property
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
    if fireOutfile != 'Not applicable':  # ST: write the fire raster-layer to the location in 'fireOutfile' if one was specified; for the simple scenario cmr_1s1d1r1k the location is 'Not applicable', because no fire-layer is used in the simple scenario
        print("Loading Phoenix4 fire grid shapefile " + firePhoenix4GridShp)
        vec = Vector.from_shapefile(firePhoenix4GridShp)
        print("Rasterising the fire vector")
        raster = vec.rasterise(fireRasterCellSize, "output = HOUR_BURNT;")
        # convert hours to seconds
        runScript("rasterised = rasterised * 3600.0;", [raster])
        print("Writing " + fireOutfile)
        raster.write(fireOutfile)

    # Rasterise the MATSim population home locations
    if popnFile != 'Not applicable':
        print("Reading MATSim population from " + popnFile)
        input = gzip.open(popnFile, 'r')
        tree = ET.parse(input)
        root = tree.getroot()
        print("Parsing home activity coordinates into a vector")
        vec = vector.Vector()
        for act in root.iter('activity'):
            if (act.attrib['type'] == popnActivityFilter):
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
    parser = ArgumentParser(description="script to convert vector and raster files using Geostack")
    parser.add_argument("--outdir", dest="outdir", type=str, help="path to the output directory - must include a forward slash ('/') after its name")
#    parser.add_argument("--scenario", dest="scenario", type=str, default='Scenario not supplied', required=True, help="name of the scenario (choose from 'mount-alexander-shire' and 'loddon-mallee-northern-cluster-shires'")
    parser.add_argument("--scenario", dest="scenario", type=str, required=True, help="name of the scenario (choose from 'Mount-Alexander-Shire', 'Loddon-Mallee-Northern-Cluster-Shires', or 'cmr_1s1d1r1k'")
    args = parser.parse_args()

#    if args.outdir is None:
    if args.scenario is None or args.outdir is None:
        parser.print_help()
    else:
        main(args.scenario, args.outdir)
