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
import shutil, time, copy, json  # ST
from geostack.vector import BoundingBox

# Run this program inside directory evac-eval/.

def main(scenario, outdir):
  # TODO: make these parameters into command-line arguments
  # outdir = './data/castlemaine-region/'

  if scenario in ('Mount-Alexander-Shire_100a', 'Mount-Alexander-Shire_100b', 'Mount-Alexander-Shire_100c', 'Mount-Alexander-Shire_100d', 'Loddon-Mallee-Northern-Cluster-Shires'):  # TODO: is it correct to use these (population and fire) files for the Loddon-Mallee scenario?
    popnFileUrl = 'https://github.com/agentsoz/ees/raw/dbfab224daaeb02294b5dabb62f55b5f8755b6ce/ees/scenarios/mount-alexander-shire/castlemaine-region-archetypes/population-archetypes.xml.gz'
    popnFile = outdir + 'abm/population-archetypes.xml.gz'
    popnEpsg = '28355'
    popnActivityFilter = 'home'
##    popnRasterCellSize = '100'
#    popnRasterCellSize = '1000'  # coarser-grained cells to facilitate definition of "population nodes" or injection-nodes: this includes each node of the road-network that lies near a population of decent size, which means the population raster must be coarse
#    popnRasterCellSize = '10000'  # very coarse-grained cells, for debugging
    popnRasterCellSizes = (100, 200, 500, 1000, 2000, 5000, 10000)
#    popnOutfile = outdir + f'sem/populationArchetypes_{int(popnRasterCellSize)}.tif'
  
    if scenario == 'Mount-Alexander-Shire_100a':
      firePhoenix4GridShpUrl = 'https://github.com/agentsoz/ees-data/raw/master/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100a/20181109_mountalex_evac_ffdi100a_grid.shp'
      firePhoenix4GridShp = outdir + 'abm/20181109_mountalex_evac_ffdi100a_grid.shp'
    elif scenario == 'Mount-Alexander-Shire_100b':
      firePhoenix4GridShpUrl = 'https://github.com/agentsoz/ees-data/raw/master/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100b/20181109_mountalex_evac_ffdi100b_grid.shp'
      firePhoenix4GridShp = outdir + 'abm/20181109_mountalex_evac_ffdi100b_grid.shp'
    elif scenario == 'Mount-Alexander-Shire_100c':
      firePhoenix4GridShpUrl = 'https://github.com/agentsoz/ees-data/raw/master/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100c/20181109_mountalex_evac_ffdi100c_grid.shp'
      firePhoenix4GridShp = outdir + 'abm/20181109_mountalex_evac_ffdi100c_grid.shp'
    elif scenario == 'Mount-Alexander-Shire_100d':
      firePhoenix4GridShpUrl = 'https://github.com/agentsoz/ees-data/raw/master/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100d/20181109_mountalex_evac_ffdi100d_grid.shp'
      firePhoenix4GridShp = outdir + 'abm/20181109_mountalex_evac_ffdi100d_grid.shp'
    elif scenario == 'Loddon-Mallee-Northern-Cluster-Shires':  # TODO: is it correct to use these (population and fire) files for the Loddon-Mallee scenario?
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


#  populationArchetypesCSVfile = 'Not applicable'  # will be redefined below if a population-archetypes CSV file exists
  populationArchetypesXMLfile = 'Not applicable'  # will be redefined below if a population-archetypes XML file exists

  if scenario == 'Mount-Alexander-Shire_100a':  # run with, e.g., python ees2sem.py --scenario Mount-Alexander-Shire_100a --outdir modelinputsMtAlexander_100a/
    fireOutfile = outdir + 'sem/20181109_mountalex_evac_ffdi100a_grid.tif'
  elif scenario == 'Mount-Alexander-Shire_100b':  # run with, e.g., python ees2sem.py --scenario Mount-Alexander-Shire_100b --outdir modelinputsMtAlexander_100b/
    fireOutfile = outdir + 'sem/20181109_mountalex_evac_ffdi100b_grid.tif'
  elif scenario == 'Mount-Alexander-Shire_100c':  # run with, e.g., python ees2sem.py --scenario Mount-Alexander-Shire_100c --outdir modelinputsMtAlexander_100c/
    fireOutfile = outdir + 'sem/20181109_mountalex_evac_ffdi100c_grid.tif'
  elif scenario == 'Mount-Alexander-Shire_100d':  # run with, e.g., python ees2sem.py --scenario Mount-Alexander-Shire_100d --outdir modelinputsMtAlexander_100d/
    fireOutfile = outdir + 'sem/20181109_mountalex_evac_ffdi100d_grid.tif'

  if scenario in ('Mount-Alexander-Shire_100a', 'Mount-Alexander-Shire_100b', 'Mount-Alexander-Shire_100c', 'Mount-Alexander-Shire_100d'):
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

  elif scenario == 'cmr_3s2d3k':  # run with, e.g., python ees2sem.py --scenario cmr_3s2d3k --outdir modelinputsCmr_3s2d3k/
    fireOutfile = 'Not applicable'  # for this scenario the location 'fireOutfile' is unspecified, because no fire-layer is necessary for this scenario
    fireRasterCellSize = 'Not applicable'
    networkFileUrl = 'scenarios/mount-alexander-shire/cmr_3s2d3k/cmr_3s2d3k_network.xml.gz'  # relative path to a local file  #(N.B. isn't a zipped-up file)
    networkFile = outdir + 'abm/cmr_3s2d3k_network.xml.gz'
#    networkFile = outdir + 'abm/cmr_3s2d3k_network.xml'  # not zipped up
    networkEpsg = '28355'
    networkOutfilePrefix = outdir + 'sem/cmr_3s2d3k_network'

  elif scenario == 'cmr1full':  # run with, e.g., python ees2sem.py --scenario cmr1full --outdir modelinputsCmr1full/
    fireOutfile = 'Not applicable'  # for this scenario the location 'fireOutfile' is unspecified, because no fire-layer is necessary for this scenario
    fireRasterCellSize = 'Not applicable'
    networkFileUrl = 'scenarios/mount-alexander-shire/cmr1full/loddon_mallee_northern_cluster_shires_network.xml.gz'  # relative path to a local file
    networkFile = outdir + 'abm/loddon_mallee_northern_cluster_shires_network.xml.gz'
    networkEpsg = '28355'
    networkOutfilePrefix = outdir + 'sem/cmr1full_network'
    assignedSubflowsOutfile = outdir + 'sem/cmr1full_assignedSubflows.geojson'

  else:
    print(f"Scenario '{scenario}' is unrecognised; terminating program.")
    sys.exit()

  # Parse args
  if popnEpsg != 'Not applicable':
    popnEpsg = get_epsg(int(popnEpsg))
#  if popnRasterCellSize != 'Not applicable':
#    popnRasterCellSize = int(popnRasterCellSize)
  networkEpsg = get_epsg(int(networkEpsg))
  epsg4326 = get_epsg(4326)

  # Create the output dir
  os.makedirs(outdir + "/abm", exist_ok=True)
  os.makedirs(outdir + "/sem", exist_ok=True)

  # Download all files
  if popnFile != 'Not applicable':
    if not os.path.exists(popnFile):
      print("Downloading " + popnFileUrl + " to " + popnFile + " ...")
      urllib.request.urlretrieve(popnFileUrl, popnFile)
    else:
      print(f"Population-file {popnFile} already exists: not downloading it.")
  if not os.path.exists(networkFile):
    if networkFileUrl[:8] == 'https://':
      print("Downloading " + networkFileUrl + " to " + networkFile + " ...")
      urllib.request.urlretrieve(networkFileUrl, networkFile)
    else:  # assume networkFileUrl gives a relative path to a local (gzipped) file
      print("Copying local file with relative path " + networkFileUrl + " to " + networkFile + " ...")
      shutil.copy(networkFileUrl, networkFile)
  else:
    print(f"Network-file {networkFile} already exists: not downloading it.")
  if firePhoenix4GridShpUrl != 'Not applicable':
    for ext in ['shp', 'prj', 'dbf', 'cpg', 'shx']:
      src = firePhoenix4GridShpUrl.replace('shp', ext)
      dst = firePhoenix4GridShp.replace('shp', ext)
      if not os.path.exists(dst):
        print("Downloading " + src + " to " + dst + " ...")
        urllib.request.urlretrieve(src, dst)
      else:
        print(f"Fire Shape-file {dst} already exists: not downloading it.")

  if scenario == 'cmr1full':  # TODO: might need to add these lines in general for all scenarios
    # For THIS scenario (perhaps also for later ones), also need to 'cp -pi ~/staticevacuation/condaenvGeostack/evac-eval/scenarios/mount-alexander-shire/cmr1full/population-archetypes.csv.gz ~/staticevacuation/condaenvGeostack/evac-eval/modelinputsCmr1full/abm/': the below four lines automate this
    populationArchetypesCSVfileUrl = 'scenarios/mount-alexander-shire/cmr1full/population-archetypes.csv.gz'  # relative path to a local file
#    populationArchetypesCSVfile = '/home/tay373/staticevacuation/condaenvGeostack/evac-eval/modelinputsCmr1full/sem/population-archetypes.csv.gz'
    populationArchetypesCSVfile = outdir + 'sem/population-archetypes.csv.gz'
    print("Copying local file with relative path " + populationArchetypesCSVfileUrl + " to " + populationArchetypesCSVfile + " ...")
    shutil.copy(populationArchetypesCSVfileUrl, populationArchetypesCSVfile)
    # For THIS scenario (perhaps also for later ones), also need to 'cp -pi ~/staticevacuation/condaenvGeostack/evac-eval/scenarios/mount-alexander-shire/cmr1full/population-archetypes.xml.gz ~/staticevacuation/condaenvGeostack/evac-eval/modelinputsCmr1full/abm/': the below four lines automate this
    populationArchetypesXMLfileUrl = 'scenarios/mount-alexander-shire/cmr1full/population-archetypes.xml.gz'  # relative path to a local file
#    popnFile = outdir + 'sem/population-archetypes.xml.gz'
    populationArchetypesXMLfile = outdir + 'sem/population-archetypes.xml.gz'
#    print("Copying local file with relative path " + populationArchetypesXMLfileUrl + " to " + popnFile + " ...")
    print("Copying local file with relative path " + populationArchetypesXMLfileUrl + " to " + populationArchetypesXMLfile + " ...")
#    shutil.copy(populationArchetypesXMLfileUrl, popnFile)
    shutil.copy(populationArchetypesXMLfileUrl, populationArchetypesXMLfile)

  # Vectorise the MATSim network
  print("Reading MATSim network from " + networkFile + " ...")
  input = gzip.open(networkFile, 'r')  # 'networkFile' must be a gzipped file
#  try:
#    input = gzip.open(networkFile, 'r')
#  except:
#    input = open(networkFile, 'r')
  tree = ET.parse(input)
  root = tree.getroot()
#  print(f"root is {root}")
  print("Parsing network nodes into a vector ...")
  nodes = dict()
  pointIDxByNodeID = dict()
  vec = vector.Vector()
  for elem in root.iter('node'):
    c = vector.Coordinate(float(elem.attrib['x']), float(elem.attrib['y']))
    pointIDx = vec.addPoint(c)
#    print(f"pointIDx is {pointIDx}, nodeID is {elem.attrib['id']}")
    vec.setProperty(pointIDx, "pointID", pointIDx)  # needed for mapping of matsim_nodeID to the vector's point-ID
    vec.setProperty(pointIDx, "matsim_nodeID", elem.attrib['id'])
    nodes[elem.attrib['id']] = c[:2]  # record the mapping of node id to xy coords
    pointIDxByNodeID[elem.attrib['id']] = pointIDx
    vec.setProperty(pointIDx, "maxCapacityOut", -99.9)  # will record the maximum capacity of any link leaving this node, i.e. any link whose 'from' property equals this node's ID
  vec.setProjectionParameters(networkEpsg)
  input.close()
  outfile = networkOutfilePrefix + ".nodes.shp"
  print("Writing " + outfile)
  vec.to_shapefile(outfile, gs_enums.GeometryType.Point)
  outfile = networkOutfilePrefix + ".nodes.geojson"
  print("Writing " + outfile)
  vec.convert(epsg4326)
  vec.to_geojson(outfile)

#  nodeIndices = vec.getPointIndexes()  # indices of all nodes
#  print(f"len(nodeIndices) is {len(nodeIndices)}")
#  for i, nodeIndex in enumerate(nodeIndices):  # check that coordinates look sensible
#    nodeCoordinates = vec.getPointCoordinate(nodeIndex).to_list()
#    print(f"nodeCoordinates from 'vec' is {nodeCoordinates}")
#    if i == 4: break  # print out only the first five nodes

#  vec_nodes = vector.Vector(vec)  # copy 'vec' to 'vec_nodes'
#  vec_nodes = copy.deepcopy(vec)  # copy 'vec' to 'vec_nodes'
  vec_nodes = vec  # ST: copy 'vec' to 'vec_nodes'

  print("Parsing network links into a vector ...")
  starttime_parseNetworkLinksToVector  = time.time()
  # vec = vector.Vector()
  for elem in root.iter('link'):
    line = [nodes[elem.attrib['from']], nodes[elem.attrib['to']]]
    lineIDx = vec.addLineString(line)
    vec.setProperty(lineIDx, "lineStringID", lineIDx)  # needed for mapping of matsim_linkID to the vector's lineString-ID
#    vec.setProperty(lineIDx, "diameter", 1.0)  # (of very low priority, as the monotonic-flow version SEM1 is not useful): "diameter" applies only to the monotonic-flow version of SEM, and even for that version it shouldn't always be 1.0
    vec.setProperty(lineIDx, "matsim_linkID", elem.attrib['id'])
    # For every version of SEM with non-monotonic flow, need capacity, free speed, number of lanes, length of each link, and whether the link is one-way:
##    print(f"elem.attrib['capacity'] is {repr(elem.attrib['capacity'])}")
#    vec.setProperty(lineIDx, "capacity", elem.attrib['capacity'])
    vec.setProperty(lineIDx, "capacity", float(elem.attrib['capacity']))  # must convert capacity-value from string to float, so that findEvacuationRisks.py can calculate 'maxlinkcapacity' raster, in line 'maxlinkcapacity = network.rasterise(rasterCellSize_large, "output = max(output, capacity);", GeometryType.LineString)' - TODO: in fact the float() here results in integer-valued capacities being saved to the vector-property "capacity", which is Geostack behaviour and probably won't cause problems but is unexpected
    vec.setProperty(lineIDx, "freespeed", elem.attrib['freespeed'])
    vec.setProperty(lineIDx, "permlanes", elem.attrib['permlanes'])
#    vec.setProperty(lineIDx, "diameter", elem.attrib['permlanes'])
    vec.setProperty(lineIDx, "length", elem.attrib['length'])
    vec.setProperty(lineIDx, "oneway", elem.attrib['oneway'])  # might not need 'oneway' property
    pointIDx = pointIDxByNodeID[ elem.attrib['from'] ]  # ASSUMPTION: every link in networkfile is one-way, i.e. every two-way link is specified as two one-way links; otherwise the 'from'-node is ill-defined
    if float(elem.attrib['capacity']) > vec_nodes.getProperty(pointIDx, "maxCapacityOut", float):
      vec_nodes.setProperty(pointIDx, "maxCapacityOut", float(elem.attrib['capacity']))
#      print(f"Setting 'maxCapacityOut' property of node with ID {elem.attrib['from']} to {float(elem.attrib['capacity'])}")
  timeToParseNetworkLinksToVector = time.time() - starttime_parseNetworkLinksToVector
  print(f"Parsing of network links from networkFile {networkFile} took {timeToParseNetworkLinksToVector:.2f} seconds.")
  vec.setProjectionParameters(networkEpsg)
#  input.close()
  outfile = networkOutfilePrefix + ".links.shp"
  print("Writing " + outfile)
  vec.to_shapefile(outfile, gs_enums.GeometryType.LineString)
  outfile = networkOutfilePrefix + ".links.geojson"
  print("Writing " + outfile)
  vec.convert(epsg4326)
  vec.to_geojson(outfile)

  # Convert the MATSim population home-locations into nearest injection-nodes, and 'EvacLocationPreference' into nearest exit-nodes:
##  if popnFile != 'Not applicable':
#  if populationArchetypesCSVfile != 'Not applicable':
  if populationArchetypesXMLfile != 'Not applicable':
    starttime_extractAssignedSubflowsAndExitNodesFromXML  = time.time()
##    (subflowsToAssignedExitNodesByInjectionNodeID, exitnodes) = extractAssignedSubflowsAndExitNodesFromCSV( populationArchetypesCSVfile )
#    print("Reading MATSim population's home-coordinates and evacuation location-preferences from " + populationArchetypesCSVfile + " ...")
    print("Reading MATSim population's home-coordinates and evacuation location-preferences from " + populationArchetypesXMLfile + " ...")
#    input = gzip.open(populationArchetypesCSVfile, 'r')
    input = gzip.open(populationArchetypesXMLfile, 'r')
    tree = ET.parse(input)
    root = tree.getroot()
    print("Parsing home-coordinates (i.e. geographical coordinates) into a vector ...")
    vec_nearestNodes = vector.Vector()
    print(f"vec_nodes projection is {vec_nodes.getProjectionParameters().to_proj4()}")
    nodeIndices = vec_nodes.getPointIndexes()  # indices of all nodes
    print(f"len(nodeIndices) is {len(nodeIndices)}")
    for i, nodeIndex in enumerate(nodeIndices):  # check that coordinates look sensible
      nodeCoordinates = vec_nodes.getPointCoordinate(nodeIndex).to_list()
      print(f"nodeCoordinates from 'vec_nodes' is {nodeCoordinates}")
      if i == 4: break  # print out only the first five nodes

    vec_homeLocations = vector.Vector()
    vec_homeLocations.setProjectionParameters(epsg4326)
    for homeCoords in root.iter('person'):
#      print(f"homeCoords is {homeCoords}")
      for attr in homeCoords.iter('attribute'):
        if (attr.attrib['name'] == 'Geographical.Coordinate'):
          print(f"attr.attrib is {attr.attrib}, attr.tag is '{attr.tag}', attr.text is '{attr.text}'")
          clist = json.loads(attr.text)  # convert string representing list of coords, e.g. '[144.2085, -37.0773]', into the equivalent list
#          print(f"clist is {clist}")
          c = vector.Coordinate(float(clist[0]), float(clist[1]))
          homeLocation_ID = vec_homeLocations.addPoint(c)
    vec_homeLocations = vec_homeLocations.convert(networkEpsg)

    homeLocationIndices = vec_homeLocations.getPointIndexes()
    print(f"len(homeLocationIndices) is {len(homeLocationIndices)}")
    for i, HLindex in enumerate(homeLocationIndices):  # check that coordinates look sensible
      HLcoords = vec_homeLocations.getPointCoordinate(HLindex).to_list()
      print(f"HLcoords from 'vec_homeLocations' is {HLcoords}")
      if i == 4: break  # print out only the first five homeLocations

#    for homeCoords in root.iter('person'):
#      for attr in homeCoords.iter('attribute'):
#        if (attr.attrib['name'] == 'Geographical.Coordinate'):
#          print(f"attr.attrib is {attr.attrib}, attr.tag is '{attr.tag}', attr.text is '{attr.text}'")
#          clist = json.loads(attr.text)  # convert string representing list of coords, e.g. '[144.2085, -37.0773]', into the equivalent list
##          print(f"clist is {clist}")
#          c = vector.Coordinate(float(clist[0]), float(clist[1]))
#          homeLocation_ID = vec_homeLocations.addPoint(c)
#          bb = BoundingBox.from_list([clist, clist])  # Geostack BoundingBox that reduces to a single point
#          nearestNode = vec_nodes.nearest(bb)  # find the nearest road-network node(s) (there may be more than one) to the current home-location
#          nearestNodeIndices = nearestNode.getPointIndexes()  # indices of all nodes that are nearest to c (there may be more than one)
#          print(f"nearestNodeIndices is {list(nearestNodeIndices)}")
#          for nodeIndex in nearestNodeIndices:
#            nodeCoordinates = nearestNode.getPointCoordinate(nodeIndex).to_list()
#            print(f"nearestNode to c is {nodeCoordinates}")
#            break  # need only one of the nearest nodes
#          if nodeCoordinates not in coordsAlreadyAddedToNearestNodes:
#            coordsAlreadyAddedToNearestNodes.append(nodeCoordinates)
#            nearestNode_ID = vec_nearestNodes.addPoint(nodeCoordinates)
##            print(f"nearestNode to c is {vec_nearestNodes.getPointCoordinate(nearestNode_ID).to_list()}")
###            vec.setProperty(lineIDx, "oneway", elem.attrib['oneway'])
###            vec.setProperty(pointIDx, "matsim_nodeID", elem.attrib['id'])
#            numHomes = vec_nearestNodes.getProperty(nearestNode_ID, "homesCount", int)
#            numHomes += 1
#            print(f"nearestNode_ID is {nearestNode_ID}")
#            vec_nearestNodes.setProperty(nearestNode_ID, "homesCount", numHomes)
#            print(f"numHomes is {numHomes}")
#      sys.exit()

#    coordsAlreadyAddedToNearestNodes = []
    nodeIndicesAlreadyAddedToNearestNodes = []
    numHLcoordsWithNearestNode = 0
    numHLcoordsWithNoNearestNode = 0
    for i, HLindex in enumerate(homeLocationIndices):
      HLcoords = vec_homeLocations.getPointCoordinate(HLindex).to_list()
      bb = BoundingBox.from_list([HLcoords, HLcoords])  # Geostack BoundingBox that reduces to a single point
      nearestNode = vec_nodes.nearest(bb)  # find the nearest road-network node(s) (there may be more than one) to the current home-location
      nearestNodeIndices = nearestNode.getPointIndexes()  # indices of all nodes that are nearest to c (there may be more than one)
      if list(nearestNodeIndices) == []:
#        print(f"HLcoords {HLcoords} has no nearest node.")
        numHLcoordsWithNoNearestNode += 1
#        nodeCoordinates = None
        nodeIndex = None
      else:
        numHLcoordsWithNearestNode += 1
        print(f"HLcoords {HLcoords} has nearestNodeIndices {list(nearestNodeIndices)}")
      for nodeIndex in nearestNodeIndices:
        nodeCoordinates = nearestNode.getPointCoordinate(nodeIndex).to_list()
#        print(f"nearestNode to c is {nodeCoordinates}")
        print(f"nodeIndex is {nodeIndex}, nodeCoordinates is {nodeCoordinates}")
        break  # need only one of the nearest nodes
##      assert len(nearestNodeIndices) in range(7)  # the largest number of nearest nodes so far found is six
#      if nodeCoordinates:
      if nodeIndex:
#        if nodeCoordinates not in coordsAlreadyAddedToNearestNodes:
        if nodeIndex not in nodeIndicesAlreadyAddedToNearestNodes:
#          coordsAlreadyAddedToNearestNodes.append(nodeCoordinates)
          nodeIndicesAlreadyAddedToNearestNodes.append( nodeIndex )
#          nearestNode_ID = vec_nearestNodes.addPoint(nodeCoordinates)
#          print(f"nearestNode to c is {vec_nearestNodes.getPointCoordinate(nearestNode_ID).to_list()}")
##          vec.setProperty(lineIDx, "oneway", elem.attrib['oneway'])
##          vec.setProperty(pointIDx, "matsim_nodeID", elem.attrib['id'])
          numHomes = 1
#          vec_nearestNodes.setProperty(nearestNode_ID, "homesCount", numHomes)
          vec_nodes.setProperty(nodeIndex, "homesCount", numHomes)
        else:
#          numHomes = vec_nearestNodes.getProperty(nearestNode_ID, "homesCount", int)
          numHomes = vec_nodes.getProperty(nodeIndex, "homesCount", int)
          numHomes += 1
#          vec_nearestNodes.setProperty(nearestNode_ID, "homesCount", numHomes)
          vec_nodes.setProperty(nodeIndex, "homesCount", numHomes)
#        print(f"i is {i}, nearestNode_ID is {nearestNode_ID}, numHomes is {numHomes}")
        print(f"i is {i}, nodeIndex is {nodeIndex}, numHomes is {numHomes}")
#        assert numHomes == 1  # not generally True
##      if i == 5:
#      if nodeCoordinates:
#        print(f"i={i} gives first non-None nodeCoordinates.")
#        sys.exit()

####    vec.setProjectionParameters(popnEpsg)
###    vec.setProjectionParameters(networkEpsg)
##    print(f"coordsAlreadyAddedToNearestNodes is {coordsAlreadyAddedToNearestNodes}")
#    print(f"nodeIndicesAlreadyAddedToNearestNodes is {nodeIndicesAlreadyAddedToNearestNodes}")
    print(f"len(nodeIndicesAlreadyAddedToNearestNodes) is {len(nodeIndicesAlreadyAddedToNearestNodes)}")
#    nearestNodesIndices = vec_nearestNodes.getPointIndexes()
#    print(f"len(nearestNodesIndices) is {len(nearestNodesIndices)}")

#    for i, NNindex in enumerate(nearestNodesIndices):  # check that coordinates look sensible
#      NNcoords = vec_nearestNodes.getPointCoordinate(NNindex).to_list()
#      numHomes = vec_nearestNodes.getProperty(NNindex, "homesCount", int)
#      print(f"i={i}: NNcoords from 'vec_nearestNodes' is {NNcoords}, with numHomes {numHomes}")
#      if i == 9: break  # print out only the first ten nearestNodes
    for i, NNindex in enumerate(nodeIndicesAlreadyAddedToNearestNodes):  # check that coordinates look sensible
      NNcoords = vec_nodes.getPointCoordinate(NNindex).to_list()
      numHomes = vec_nodes.getProperty(NNindex, "homesCount", int)
      print(f"i={i}: NNindex {NNindex} from 'vec_nodes' has NNcoords {NNcoords}, with numHomes {numHomes}")
#      if i == 9: break  # print out only the first ten nearestNodes

    print(f"numHLcoordsWithNearestNode is {numHLcoordsWithNearestNode}")
    print(f"numHLcoordsWithNoNearestNode is {numHLcoordsWithNoNearestNode}")

    input.close()
#    print("Rasterising the home activity locations vector ...")
#    bounds = vec.getBounds()
#    bounds.extend(popnRasterCellSize)
#    count = Raster(name="count", data_type=np.uint32)
#    count.init_with_bbox(bounds, popnRasterCellSize)
#    count.setProjectionParameters(popnEpsg)
#    count.setAllCellValues(0)
#    count.rasterise(vec, "atomic_inc();")
#    print("Writing " + popnOutfile)
    print("Writing " + assignedSubflowsOutfile)
#    count.write(popnOutfile)
#    print("Writing " + outfile)
    vec.convert(epsg4326)
    vec.to_geojson( assignedSubflowsOutfile )
#    endtime_extractAssignedSubflowsAndExitNodesFromXML = time.time()
#    timeToExtractAssignedSubflowsAndExitNodesFromXML = endtime_extractAssignedSubflowsAndExitNodesFromXML - starttime_extractAssignedSubflowsAndExitNodesFromXML
    timeToExtractAssignedSubflowsAndExitNodesFromXML = time.time() - starttime_extractAssignedSubflowsAndExitNodesFromXML
    print(f"Extraction of assigned subflows and exit-nodes from populationArchetypesXMLfile {populationArchetypesXMLfile} took {timeToExtractAssignedSubflowsAndExitNodesFromXML:.5f} seconds.")

  # Rasterise the fire shapefile
  if fireOutfile != 'Not applicable':  # ST: write the fire raster-layer to the location in 'fireOutfile' if one was specified; for the simple scenario cmr_1s1d1r1k the location is 'Not applicable', because no fire-layer is used in the simple scenario
    print("Loading Phoenix4 fire grid shapefile " + firePhoenix4GridShp + " ...")
    vec = Vector.from_shapefile(firePhoenix4GridShp)
    print("Rasterising the fire vector ...")
    raster = vec.rasterise(fireRasterCellSize, "output = HOUR_BURNT;")
    # convert hours to seconds
    runScript("rasterised = rasterised * 3600.0;", [raster])
    print("Writing " + fireOutfile)
    raster.write(fireOutfile)

  # Rasterise the MATSim population home locations
  if popnFile != 'Not applicable':
    print("Reading MATSim population from " + popnFile + " ...")
    input = gzip.open(popnFile, 'r')
    tree = ET.parse(input)
    root = tree.getroot()
    print("Parsing home activity coordinates into a vector ...")
    vec = vector.Vector()
    for act in root.iter('activity'):
      if (act.attrib['type'] == popnActivityFilter):
        c = vector.Coordinate(float(act.attrib['x']), float(act.attrib['y']))
        pointIDx = vec.addPoint(c)
#        vec.setProperty(pointIDx, "newproperty", "newstr")
    vec.setProjectionParameters(popnEpsg)
    input.close()
    print("Rasterising the home activity locations vector ...")
    bounds = vec.getBounds()
    for popnRasterCellSize in popnRasterCellSizes:
      bounds.extend(popnRasterCellSize)
      count = Raster(name="count", data_type=np.uint32)
      count.init_with_bbox(bounds, popnRasterCellSize)
      count.setProjectionParameters(popnEpsg)
      count.setAllCellValues(0)
      count.rasterise(vec, "atomic_inc();")
      popnOutfile = outdir + f'sem/populationArchetypes_{int(popnRasterCellSize)}.tif'
      print("Writing " + popnOutfile)
      count.write(popnOutfile)



if __name__ == "__main__":
  parser = ArgumentParser(description="script to convert vector and raster files using Geostack")
  parser.add_argument("--outdir", dest="outdir", type=str, help="path to the output directory - must include a forward slash ('/') after its name")
#  parser.add_argument("--scenario", dest="scenario", type=str, default='Scenario not supplied', required=True, help="name of the scenario (choose from 'mount-alexander-shire' and 'loddon-mallee-northern-cluster-shires'")
  parser.add_argument("--scenario", dest="scenario", type=str, required=True, help="name of the scenario (choose from 'Mount-Alexander-Shire_100[a-d]', 'Loddon-Mallee-Northern-Cluster-Shires', 'cmr_1s1d1r1k', 'cmr_3s2d3k', or 'cmr1full'")
  args = parser.parse_args()

#  if args.outdir is None:
  if args.scenario is None or args.outdir is None:
    parser.print_help()
  else:
    main(args.scenario, args.outdir)
