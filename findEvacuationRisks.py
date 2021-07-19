import os
import os.path as pth
import sys
#sys.path.append("/home/gar305/Documents/geostack/py_build/python/src")
import json
from time import time
#from datetime import datetime
import numpy as np
from argparse import ArgumentParser
#from geostack.solvers import NetworkFlowSolver
from geostack.raster import Raster, RasterFile
from geostack.vector import BoundingBox#, Vector
#from geostack.core import ProjectionParameters
from geostack.io import geoJsonToVector, vectorToGeoJson
from geostack.runner import runScript
from geostack.gs_enums import GeometryType#, ReductionType
import solversSEMversions
import scipy.linalg

# Check parameter types and values
def getParam(cfgJson, param, type_options=None, checkGreater=None, default=None):

  # Return default if not found
  if param not in cfgJson:
    return default

  # Check type
  if type_options is not None:
    if not isinstance(type_options, tuple):
      raise TypeError("type_options should be tuple of datatype")

    for item in type_options:
      if item not in [int, float, str, dict, list]:
        raise TypeError(f"type_option {item.__name__} is not recognised")

    if not isinstance(cfgJson[param], type_options):
      data_type = ','.join([f"{item.__name__}" for item in type_options])
      raise TypeError(f"{param} should be {data_type}")

  # Check value
  if checkGreater is not None:
    if cfgJson[param] <= checkGreater:
      raise ValueError(f"{param} should be greater than {checkGreater}")

  return cfgJson[param]

def findRisks(cfgJson):

#  print("===Geostack evacuation-risks identifier===")
  startTime = time()

  # Check json entries
  try:
#    evacuationLayer = getParam(cfgJson, "evacuationLayer", (str,))
    populationLayer = getParam(cfgJson, "populationLayer", (str,))
#    networkGeoJSON = getParam(cfgJson, "networkGeoJSON", (str,))
    JSONnetworkfilename = getParam(cfgJson, "networkGeoJSON", (str,))
    outputGeoJSONfilename = getParam(cfgJson, "outputGeoJSON", (str,))
#    networkConfig = getParam(cfgJson, "networkConfig", (dict,))
    FIRE_BOUNDS_EXTENSION_m = getParam(cfgJson, "FIRE_BOUNDS_EXTENSION_m", (float,))  # distance in metres by which the bounds of one fire-layer are extended in order to encompass all relevant fire-layers; TODO: find a better general method
    print(f"FIRE_BOUNDS_EXTENSION_m is {FIRE_BOUNDS_EXTENSION_m} metres")
    CELLSIZEFINE_m = getParam(cfgJson, "CELLSIZEFINE_m", (float,))  # in metres
#    CELLSIZESMALL_m = getParam(cfgJson, "CELLSIZESMALL_m", (float,))  # in metres
    CELLSIZELARGE_m = getParam(cfgJson, "CELLSIZELARGE_m", (float,))  # in metres
    INNERDIST_SAFEBUFFER_m = getParam(cfgJson, "INNERDIST_SAFEBUFFER_m", (float,))  # distance in metres from fire-perimeter to inner boundary of the "safe" buffer in which exit-nodes might be placed
    OUTERDIST_SAFEBUFFER_m = getParam(cfgJson, "OUTERDIST_SAFEBUFFER_m", (float,))  # distance in metres from fire-perimeter to outer boundary of the "safe" buffer in which exit-nodes might be placed

#    CELLSIZEFINE_m = 80.0

#    CELLSIZESMALL_m = 100.0

#    CELLSIZELARGE_m = 1000.0
##    CELLSIZELARGE_m = 10000.0
#    print(f"CELLSIZEFINE_m is {CELLSIZEFINE_m}, CELLSIZESMALL_m is {CELLSIZESMALL_m}, CELLSIZELARGE_m is {CELLSIZELARGE_m}")
    print(f"CELLSIZEFINE_m is {CELLSIZEFINE_m} metres, CELLSIZELARGE_m is {CELLSIZELARGE_m} metres")

#    rasterCellSize_small = CELLSIZESMALL_m
#    rasterCellSize_large = CELLSIZELARGE_m

#    INNERDIST_SAFEBUFFER_m = 1000.0  # distance in metres from fire-perimeter to inner boundary of the "safe" buffer in which exit-nodes might be placed
#    OUTERDIST_SAFEBUFFER_m = 2000.0  # distance in metres from fire-perimeter to outer boundary of the "safe" buffer
    print(f"INNERDIST_SAFEBUFFER_m is {INNERDIST_SAFEBUFFER_m} metres, OUTERDIST_SAFEBUFFER_m is {OUTERDIST_SAFEBUFFER_m} metres")
#    sys.exit()

  except Exception as e:

    # Return error if any parameters are incorrect
    print("Parameter error: " + str(e))
    return -1

  outputDir = pth.dirname(pth.abspath(outputGeoJSONfilename))
  print(f"outputDir is {outputDir}")

  # Read population input-layer:
  print(f"Reading population_unaligned data from {populationLayer}...")
  if pth.splitext(populationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:
    population_unaligned = Raster(name="population_unaligned")
    population_unaligned.read(populationLayer)
  else:
    population_unaligned = RasterFile(name = "population_unaligned", filePath = populationLayer, backend = 'gdal')
    population_unaligned.read()
#  population_unaligned_bounds = population_unaligned.getBounds()
  dims = population_unaligned.getDimensions()
  print(f"population_unaligned dimensions is {dims}")
  population_unaligned.write(pth.join(outputDir, f"population_unaligned.tif"))

  # Read just one fire-layer to obtain its projection, and extend its bounds to define bounds encompassing all fires:
#  evacuationLayer = 'modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif'
  evacuationLayer = 'modelinputsMtAlexander_100c/sem/20181109_mountalex_evac_ffdi100c_grid.tif'
  print(f"To obtain projection, reading fire data from {evacuationLayer}...")
  if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:  # use native raster-readers
    fire = Raster(name="fire")
    fire.read(evacuationLayer)
  else:  # use gdal-based raster-reader
    fire = RasterFile(name = "fire", filePath = evacuationLayer, backend = 'gdal')
    fire.read()
  allFires_proj = fire.getProjectionParameters()
#  print(f"allFires_proj is {allFires_proj}")
  allFires_bounds = fire.getBounds()
#  allFires_bounds.extend(20000.0)  # extend bounds of the fire layer by 20 km
  allFires_bounds.extend(FIRE_BOUNDS_EXTENSION_m)  # extend bounds of the fire layer; TODO: the purpose of this extension is to cover all fires; in general, it might be better to calculate a union of the fire-rasters

  # Read road-network vector-layer:
  print(f"Reading road-network data from {JSONnetworkfilename} ...")
  roadNetwork = geoJsonToVector(JSONnetworkfilename)
  network_proj = roadNetwork.getProjectionParameters()
  # Restrict the network's bounds to a smaller area around the fire, then convert these bounds to the network projection:
#  network_bounds = fire.getBounds()
##  network_bounds = allFires_bounds  # bug: allFires_bounds will be changed by 'extend' in the next line
#  network_bounds.extend(30000.0)
  network_bounds = allFires_bounds
#  network_bounds.extend(10000.0)  # extend bounds of the layer by 10 km
  print(f"allFires_bounds is {allFires_bounds}")
  print(f"network_bounds is {network_bounds}")
  assert network_bounds == allFires_bounds

#  networkBounds0 = roadNetwork.getBounds()
#  print(f"roadNetwork.getBounds() is {networkBounds0}")

  network_bounds = network_bounds.convert(network_proj, allFires_proj)
  print(f"network_bounds after conversion to network_proj is {network_bounds}")

  roadNetwork = roadNetwork.region(network_bounds)  # clip roadNetwork to allFires_bounds
  networkBounds1 = roadNetwork.getBounds()
  print(f"roadNetwork.getBounds() after clipping to roadNetwork.region(network_bounds) is {networkBounds1}")

  roadNetwork = roadNetwork.convert(allFires_proj)  # convert roadNetwork to same projection as the 'fire' raster, for point-sampling
#  networkBounds2 = roadNetwork.getBounds()
#  print(f"roadNetwork.getBounds() after conversion to allFires_proj is {networkBounds2}")  # roadNetwork bounds are now wide again - a bug in Geostack?

  print("Creating 'numFiresInside' layer ...")  # to count the number of fire-perimeters each cell is inside
#  numFiresInside = Raster(name="numFiresInside", data_type=np.uint32)
#  numFiresInside = Raster(name="numFiresInside", data_type=np.int32)
  numFiresInside = Raster(name="numFiresInside")
#  ifp_bounds = fire.getBounds()
##  ifp_bounds.extend(30000.0)  # extend bounds of the layer by 30 km
#  ifp_bounds.extend(20000.0)  # extend bounds of the layer by 20 km
##  print(f"allFires_bounds is {allFires_bounds}")
#  ifp_bounds = allFires_bounds
##  print(f"ifp_bounds is {ifp_bounds}")
#  numFiresInside.init_with_bbox(population_unaligned_bounds, CELLSIZESMALL_m)
#  numFiresInside.init_with_bbox(network_bounds, CELLSIZESMALL_m)
#  numFiresInside.init_with_bbox(ifp_bounds, CELLSIZESMALL_m)
  numFiresInside.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
  numFiresInside.setProjectionParameters(allFires_proj)
  numFiresInside.setAllCellValues(0)

  # Read all fire-layers to create 'numFiresInside' layer:
##  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif')
#  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif', 'modelinputsMtAlexander_100c/sem/20181109_mountalex_evac_ffdi100c_grid.tif', 'modelinputsMtAlexander_100d/sem/20181109_mountalex_evac_ffdi100d_grid.tif')
  allFireNames = ('a', 'b', 'c', 'd')  # fires ffdi100a, ffdi100b, ...
#  evacuationLayers = [f'modelinputsMtAlexander_100{es}/sem/20181109_mountalex_evac_ffdi100{es}_grid.tif' for es in allFireNames]
#  print(f"evacuationLayers is {evacuationLayers}")
#  for evacuationLayer in evacuationLayers:
  for fireName in allFireNames:
    evacuationLayer = f'modelinputsMtAlexander_100{fireName}/sem/20181109_mountalex_evac_ffdi100{fireName}_grid.tif'
    print(f"Reading fire data from {evacuationLayer}...")
    if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:
      # use native raster-readers:
##      currentFire = Raster(name=f"fire_{evacuationLayer}")
#      currentFire = Raster(name=f"fire")
#      currentFire = Raster(name=f"fire_fireName")
      currentFire = Raster(name=f"currentFire")
      currentFire.read(evacuationLayer)
    else:
      # use gdal based raster reader:
#      currentFire = RasterFile(name = f"fire_{evacuationLayer}", filePath = evacuationLayer, backend = 'gdal')
#      currentFire = RasterFile(name = f"fire", filePath = evacuationLayer, backend = 'gdal')
#      currentFire = RasterFile(name = f"fire_fireName", filePath = evacuationLayer, backend = 'gdal')
      currentFire = RasterFile(name = f"currentFire", filePath = evacuationLayer, backend = 'gdal')
      currentFire.read()
##    fire_proj = currentFire.getProjectionParameters()
##    print(f"fire_proj is {fire_proj}")
#    fire_max_time = np.nanmax(currentFire.data)  # time that defines perimeter of fire's maximum extent
#    currentFire.setVariableData("fireMaxTime", fire_max_time)  # set fireMaxTime in 'currentFire' layer
#    runScript('numFiresInside += currentFire < currentFire::fireMaxTime ? 1 : 0;', [numFiresInside, currentFire] )  # add 1 to value of cell if inside the fire-perimeter, 0 otherwise
    runScript('numFiresInside += currentFire > 0 ? 1 : 0;', [numFiresInside, currentFire] )  # add 1 to value of cell if inside the fire-perimeter, 0 otherwise

  dims = numFiresInside.getDimensions()
  print(f"numFiresInside dimensions is {dims}")
  numFiresInside.write(pth.join(outputDir, f"numFiresInside_{int(CELLSIZELARGE_m)}.tif"))  # write raster-layer 'numFiresInside' to a file

#  numFiresPopulationIsInside = Raster(name="numFiresPopulationIsInside", data_type=np.int32)
  numFiresPopulationIsInside = Raster(name="numFiresPopulationIsInside")
  numFiresPopulationIsInside.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
  numFiresPopulationIsInside.setProjectionParameters(allFires_proj)
#  runScript("numFiresPopulationIsInside = population_unaligned > 0 ? numFiresInside : noData_UINT;", [numFiresPopulationIsInside, population_unaligned, numFiresInside])  # results in very large raster-cell values, for some reason
  runScript("numFiresPopulationIsInside = population_unaligned > 0 ? numFiresInside : -1;", [numFiresPopulationIsInside, population_unaligned, numFiresInside])  # -1 indicates population_unaligned is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to; numFiresPopulationIsInside is 'nan' where population_unaligned > 0 and numFiresInside is 'nan'


# Using lines like the following with 'output' seems to result in 'numFiresPopulationIsInside' raster having a number of bands equal to 'CELLSIZELARGE_m' - why?
###  numFiresPopulationIsInside = runScript("output = population_unaligned * numFiresInside;", [population_unaligned, numFiresInside])
###  numFiresPopulationIsInside = runScript("output = population_unaligned && numFiresInside;", [population_unaligned, numFiresInside])
##  numFiresPopulationIsInside = runScript("output = population_unaligned > 0 ? numFiresInside : 0;", [population_unaligned, numFiresInside])
#  numFiresPopulationIsInside = runScript("output = population_unaligned > 0 ? numFiresInside : -1;", [population_unaligned, numFiresInside])  # -1 indicates population_unaligned is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to

  dims_numFiresPopulationInside = numFiresPopulationIsInside.getDimensions()
#  print(f"dims_numFiresPopulationInside is {dims_numFiresPopulationInside}")
  numFiresPopulationIsInside.write(pth.join(outputDir, f"numFiresPopulationIsInside_{int(CELLSIZELARGE_m)}.tif"))  # write raster 'numFiresPopulationIsInside' to a file

  population = runScript("output = numFiresPopulationIsInside;", [numFiresPopulationIsInside])  # make new raster with same dimensions as numFiresPopulationIsInside, to store population values; copies the contents of numFiresPopulationIsInside rather than simply making a link to the raster, so that runScript on population will NOT change the contents of numFiresPopulationIsInside
  population.name = 'population'
  runScript("population = population_unaligned;", [population, population_unaligned])  # set population to a raster aligned with numFiresPopulationIsInside
  dims_population = population.getDimensions()
  print(f"dims_population is {dims_population}")
  assert dims_population == dims_numFiresPopulationInside
  population.write(pth.join(outputDir, f"populationArchetypes_{int(CELLSIZELARGE_m)}.tif"))  # write 'population' to a file

#  sys.exit()

#  roadNetwork.pointSample(numFiresInside)  # sample 'numFiresInside' raster at each point and write resulting value to the 'roadNetwork' vector-layer
#  with open(pth.join(outputDir, f"clippedNetwork_{int(CELLSIZELARGE_m)}.geojson"), "w") as f:
#    f.write(vectorToGeoJson(roadNetwork))


  # For each fire, define injection- and exit-nodes, then calculate maximum flow from injection- to exit-nodes:
  meanTimeToEvacuateByPopulationNode = dict()  # mean time, in hours, taken to evacuate population-node, over all fires that contain the node within their perimeters
  numPopulationNodesByFire = dict()  # number of positive-population population-nodes within fire's perimeter
  numCriticalLinksByFire = dict()  # number of critical links within fire's perimeter
  numEvacuationsLinkIsCriticalTo = dict()
  # Read fire-raster input-layers:
##  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', )  # just one fire, for simplicity
##  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif')
##  for evacuationLayer in evacuationLayers:
#  evacuationScenarios = ('a', 'b', 'c', 'd')  # fires ffdi100a, ffdi100b, ...
  evacuationScenarios = ('a', )  # the single fire ffdi100a
  for fireName in evacuationScenarios:
    evacuationLayer = f'modelinputsMtAlexander_100{fireName}/sem/20181109_mountalex_evac_ffdi100{fireName}_grid.tif'
    print(f"################# Reading currentFire data from {evacuationLayer}...")
    if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:
      # use native raster-readers:
      currentFire = Raster(name="currentFire")
      currentFire.read(evacuationLayer)
    else:
      # use gdal-based raster-reader:
      currentFire = RasterFile(name = "currentFire", filePath = evacuationLayer, backend = 'gdal')
      currentFire.read()
#    fire_proj = currentFire.getProjectionParameters()


#    # Get fire-arrival-time bounds:
#    fire_bounds = currentFire.getBounds()
##    fire_bounds.extend(10000.0)  # extend bounds of the currentFire layer by 10 km

#    # Create initialFireExtent layer, to hold the area of the fire's epicentre:
##    print("Creating 'initialFireExtent' layer ...")
#    initialFireExtent = Raster(name="initialFireExtent")
#    initialFireExtent.init_with_bbox(allFires_bounds, CELLSIZESMALL_m)
#    initialFireExtent.setProjectionParameters(allFires_proj)
#    fire_time_epicentre = 1.0  # time in seconds used to define fire's epicentre
#    print(f"Finding fire's epicentre by running for time (fire-duration) {fire_time_epicentre}...")
#    currentFire.setVariableData("time", fire_time_epicentre)  # set current time in 'currentFire' layer
#    runScript("initialFireExtent = currentFire < currentFire::time ? 1 : 0;", [initialFireExtent, currentFire])
##    initialFireExtent.write(pth.join(outputDir, f"firePerimeter_{int(fire_time_epicentre):06d}.tif"))

##    print("Creating 'contourOfInitialFireExtent' layer ...")
#    contourOfInitialFireExtent = initialFireExtent.vectorise( [ 0.5 ] )
#    with open(pth.join(outputDir, f"contourOfInitialFireExtent_{int(fire_time_epicentre):06d}.geojson"), "w") as contourfile:
#      contourfile.write(vectorToGeoJson(contourOfInitialFireExtent))

#    perimeterBounds = contourOfInitialFireExtent.getBounds()
#    print(f"perimeterBounds[0] is {repr(perimeterBounds[0])}")  # returns nothing useful
##    bb = BoundingBox.from_list([[0.5, 0.5], [0.5, 0.5]])  # create search-point
##    nearestNodes = roadNetwork.nearest(bb)  # find node nearest to epicentre
#    nearestNodes = roadNetwork.nearest(perimeterBounds)  # find node nearest to epicentre
#    nearestNodeIndices = nearestNodes.getPointIndexes()  # indices of all nodes that are nearest to epicentre (there may be more than one)
#    print(f"Fire epicentre has len(nearestNodeIndices)={len(nearestNodeIndices)}")
#    print(f"nearestNodeIndices is {list(nearestNodeIndices)}")
#    propertyNames = nearestNodes.getProperties().getPropertyNames()
#    print(f"nearestNodes' propertyNames is {propertyNames}")
##    nearestNodeGeometryIndices = nearestNodes.getGeometryIndexes()  # includes links as well as nodes, so not needed: use nearestNodeIndices instead
##    print(f"nearestNodeGeometryIndices is {list(nearestNodeGeometryIndices)}")
##    for i in nearestNodeGeometryIndices:
#    for i in nearestNodeIndices:
#      print(f"{i} {dict([ (p, nearestNodes.getProperty(i, p)) for p in propertyNames ])}")
#    if nearestNodeIndices:
#      nodeIndexEpicentreSnapsTo = nearestNodeIndices[0]
#      print(f"nodeIndexEpicentreSnapsTo is {nodeIndexEpicentreSnapsTo}")
#      MATsimID_nodeEpicentreSnapsTo = nearestNodes.getProperty( nodeIndexEpicentreSnapsTo, 'matsim_nodeID' )
#      print(f"MATsimID_nodeEpicentreSnapsTo is {MATsimID_nodeEpicentreSnapsTo}")
#      with open(pth.join(outputDir, f"nearestNodesToFireEpicentre_{int(fire_time_epicentre):06d}.geojson"), "w") as f:
#        f.write(vectorToGeoJson(nearestNodes))
#    if list(nearestNodeIndices) == []:
##      print(f"HLcoords {HLcoords} has no nearest node.")
###      nodeCoordinates = None
#      nodeIndex = None
#    else:
###      print(f"Fire epicentre, with HLcoords {HLcoords} has nearestNodeIndices {list(nearestNodeIndices)}")
#      pass
#    for nodeIndex in nearestNodeIndices:
#      nodeCoordinates = nearestNodes.getPointCoordinate(nodeIndex).to_list()
##      print(f"nearestNodes to c is {nodeCoordinates}")
#      print(f"Nearest node to currentFire epicentre: nodeIndex is {nodeIndex}, nodeCoordinates is {nodeCoordinates}")
#      break  # need only one of the nearest nodes
##    sys.exit()


#    if evacuationLayer == 'modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif':
#      inflowsByNodeID = {'619133802': 500, '618559403': 800, '357004729': 300}  # these values for populations within raster-cells are estimates
#      print(f"inflowsByNodeID is {inflowsByNodeID}")
#      exitnodes = {'3511596886', '357001685'}  # for now, selected by hand using QGIS; TODO: ought to lie inside the band lying between one and two kilometres from the fire's maximum-extent perimeter
#    elif evacuationLayer == 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif':
##      inflowsByNodeID = {'294313111': 14}  # value of 14 obtained by using QGIS's 'Identify Features' on 'population-archetypes' layer at node '294313111'
#      inflowsByNodeID = {'294313111': 800, '402305951': 600}  # these values for populations within raster-cells are estimates
#      print(f"inflowsByNodeID is {inflowsByNodeID}")
#      exitnodes = {'316917533'}  # for now, selected by hand using QGIS, from the band lying between one and two kilometres from the fire's maximum-extent perimeter


#    fire_max_time = np.nanmax(currentFire.data)  # time that defines perimeter of currentFire's maximum extent
##    print(f"Running for time (fire-duration) {fire_max_time}...")
#    print(f"fire_max_time is {fire_max_time}")

#    # Create inflowPerNodeInsidePerimeter layer:
##    print("Creating 'inflowPerNodeInsidePerimeter' layer ...")
#    inflowPerNodeInsidePerimeter = Raster(name="inflowPerNodeInsidePerimeter")
#    inflowPerNodeInsidePerimeter.init_with_bbox(allFires_bounds, CELLSIZESMALL_m)
#    inflowPerNodeInsidePerimeter.setProjectionParameters(allFires_proj)
#    currentFire.setVariableData("fireMaxTime", fire_max_time)  # set fireMaxTime in 'currentFire' layer
#    runScript("inflowPerNodeInsidePerimeter = currentFire < currentFire::fireMaxTime ? 1 : 0;", [inflowPerNodeInsidePerimeter, currentFire])
##    inflowPerNodeInsidePerimeter.write(pth.join(outputDir, f"inflowPerNodeInsidePerimeter_{int(fire_max_time):06d}.tif"))

    # Create isInsideFirePerimeter layer:
#    print("Creating 'isInsideFirePerimeter' layer ...")
    isInsideFirePerimeter = Raster(name="isInsideFirePerimeter")
#    isInsideFirePerimeter.init_with_bbox(allFires_bounds, CELLSIZESMALL_m)
    isInsideFirePerimeter.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
    isInsideFirePerimeter.setProjectionParameters(allFires_proj)
#    currentFire.setVariableData("fireMaxTime", fire_max_time)  # set fireMaxTime in 'currentFire' layer
#    runScript("isInsideFirePerimeter = currentFire < currentFire::fireMaxTime ? 1 : 0;", [isInsideFirePerimeter, currentFire])
    runScript("isInsideFirePerimeter = currentFire > 0 ? 1 : 0;", [isInsideFirePerimeter, currentFire])
    isInsideFirePerimeter.write(pth.join(outputDir, f"isInsideFirePerimeter_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.tif"))  # write raster 'isInsideFirePerimeter' to a file

#    sys.exit()

#    print("Creating 'maxlinkcapacity_unaligned' layer ...")  # get maximum link-capacity within each cell
###    maxlinkcapacity_unaligned = Raster(name="maxlinkcapacity_unaligned")
##    maxlinkcapacity_unaligned = Raster(name="maxlinkcapacity_unaligned", data_type=np.uint32)
##    maxlinkcapacity_unaligned.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
##    maxlinkcapacity_unaligned.setProjectionParameters(allFires_proj)
##    maxlinkcapacity_unaligned.setAllCellValues(0)
#    # James Hilton, 17th June 2021:
#    # > The atomic_max function only works for integers in OpenCL (I have no idea why). If it did work I think you could just do: r.rasterise(v, "atomic_max(capacity);") where r is your count raster, v the vector and capacity the field name. For this capacity would have to be an integer. I was thinking you could do this with a vector script but I think the same issue will occur.
##    maxlinkcapacity_unaligned.rasterise(roadNetwork, "atomic_max('capacity');", GeometryType.LineString)  # find maximum link-capacity within each cell; doesn't work, for some reason, and runs even if 'capacity' is replaced with a typo such as 'cpacity'
#    # 'roadNetwork' vector has much wider bounds than 'numFiresInside', so clip 'roadNetwork' to the right size:
#    roadNetwork = roadNetwork.region(allFires_bounds)  # clip roadNetwork again, but to allFires_bounds, which are now in the same projection; this makes the bounds similar to those in allFires_bounds, so that rasterising the roadNetwork to create maxlinkcapacity_unaligned will give the raster-layer suitable bounds
#    networkBounds3 = roadNetwork.getBounds()
#    print(f"roadNetwork.getBounds() after clipping to roadNetwork.region(allFires_bounds) is {networkBounds3}")
##    maxlinkcapacity_unaligned = roadNetwork.rasterise(CELLSIZELARGE_m, "output = max(output, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = roadNetwork.rasterise(CELLSIZELARGE_m, "output = max(output, capacity);")  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = roadNetwork.rasterise(CELLSIZELARGE_m, "output = max(capacity, output);", 2)  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = roadNetwork.rasterise(CELLSIZELARGE_m, "output = min(5000, capacity);", GeometryType.LineString)  # gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = roadNetwork.rasterise(CELLSIZELARGE_m, "output = max(0, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run; this line doesn't seem to work unless the raster-cell size is small (it works certainly for a cell-size between 80 and 640)
#    maxlinkcapacity_unaligned = roadNetwork.rasterise(80, "output = max(0, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell: this works, at its fine-grained scale of 80
##    maxlinkcapacity_unaligned = roadNetwork.rasterise(640, "output = max(0, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell; as the scale is increased this line works less well: the higher-capacity links get "overwritten" by lower-capacity ones
#    maxlinkcapacity_unaligned.write(pth.join(outputDir, f"maxlinkcapacity_unaligned_{int(CELLSIZELARGE_m)}.tif"))
#    # TODO: 'maxlinkcapacity_unaligned' raster has different dimensions from, and is not aligned with, 'numFiresInside'; must find a way to align it with 'numFiresInside'
#    maxlinkcapacity_unaligned.name = 'maxlinkcapacity_unaligned'
###    maxlinkcapacity = numFiresPopulationIsInside  # make new raster with same dimensions as numFiresPopulationIsInside, to store maxlinkcapacity values; this line is buggy, though, because it makes the new raster simply a link to the existing raster, and so runScript on maxlinkcapacity will change the contents of numFiresPopulationIsInside
###    maxlinkcapacity = numFiresPopulationIsInside.clip(numFiresPopulationIsInside)  # this line is buggy, because it makes the new raster simply a link to the existing raster
#    maxlinkcapacity = runScript("output = numFiresPopulationIsInside;", [numFiresPopulationIsInside])  # make new raster with same dimensions as numFiresPopulationIsInside, to store maxlinkcapacity values; copies the contents of numFiresPopulationIsInside rather than simply making a link to the raster, so that runScript on maxlinkcapacity will NOT change the contents of numFiresPopulationIsInside
#    maxlinkcapacity.name = 'maxlinkcapacity'
#    runScript("maxlinkcapacity = maxlinkcapacity_unaligned;", [maxlinkcapacity, maxlinkcapacity_unaligned])  # set maxlinkcapacity to a raster aligned with numFiresPopulationIsInside that contains the maximum link-capacity within each cell
#    dims_maxlinkcapacity = maxlinkcapacity.getDimensions()
#    print(f"dims_maxlinkcapacity is {dims_maxlinkcapacity}")
#    assert dims_maxlinkcapacity == dims_numFiresPopulationInside
#    maxlinkcapacity.write(pth.join(outputDir, f"maxlinkcapacity_{int(CELLSIZELARGE_m)}.tif"))  # write 'maxlinkcapacity', actually containing maximum link-capacity per cell, to a file
#    print(f"maxlinkcapacity.max() is {maxlinkcapacity.max()}")
##    sys.exit()
#
#    maxlinkcapacity_proj = maxlinkcapacity.getProjectionParameters()
#    assert maxlinkcapacity_proj == allFires_proj
#
###    network_bounds = network_bounds.convert(network_proj, allFires_proj)
###    maxlinkcapacity = maxlinkcapacity.region(network_bounds)  # clip raster to the bounds: doesn't work, because seems to re-initialise the raster
##    maxlinkcapacity = maxlinkcapacity.region(allFires_bounds)  # clip raster to the bounds: doesn't work, because seems to re-initialise the raster
###    maxlinkcapacity = maxlinkcapacity.convert(allFires_proj)  # convert maxlinkcapacity to same projection as the 'currentFire' raster; doesn't work, as "'Raster' object has no attribute 'convert'"
##    dims = maxlinkcapacity.getDimensions()
##    print(f"maxlinkcapacity dimensions is {dims}")
##    maxlinkcapacity.write(pth.join(outputDir, f"maxlinkcapacity_{int(CELLSIZELARGE_m)}.tif"))

    populationIsInsideFire = Raster(name="populationIsInsideFire")
    populationIsInsideFire.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
    populationIsInsideFire.setProjectionParameters(allFires_proj)
#    runScript("populationIsInsideFire = population_unaligned > 0 ? isInsideFirePerimeter : -1;", [populationIsInsideFire, population_unaligned, isInsideFirePerimeter])  # -1 indicates population_unaligned is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to; populationIsInsideFire is 'nan' where population_unaligned > 0 and isInsideFirePerimeter is 'nan'
    runScript("populationIsInsideFire = population > 0 ? isInsideFirePerimeter : -1;", [populationIsInsideFire, population, isInsideFirePerimeter])  # -1 indicates population is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to; populationIsInsideFire is 'nan' where population > 0 and isInsideFirePerimeter is 'nan'
    populationIsInsideFire.write(pth.join(outputDir, f"populationIsInsideFire_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.tif"))  # write raster 'populationIsInsideFire' to a file

    # Define the injection-nodes: there is at most one in each raster-cell that both has positive population and lies within currentFire's perimeter: for each such raster-cell, the injection-node is the road-network node, if any, that has the largest maximum out-capacity:
    # We make the assumption that every network-node having the largest maximum out-capacity is equally suitable as an injection-node; in reality this is not always true, for example in the case of a highway - a highway has nearby nodes on links that go in opposite directions, which means a cell's population might "snap to" the highway node that takes traffic in the wrong direction, away from the exit-nodes, and so will have to double back.
    print("######## Defining the injection-nodes ...")
##    roadNetwork.addProperty("numFirePerimetersInside")
#    dims = numFiresPopulationIsInside.getDimensions()
#    print(f"numFiresPopulationIsInside dimensions is {dims}")
    dims = populationIsInsideFire.getDimensions()
    print(f"populationIsInsideFire dimensions is {dims}")
#    NEGATIVEVALUE = -99.9
    inflowsByNodeID = dict()
    positivePopulationInsideFireByNodeID = dict()
    inflowCoordsByNodeID = dict()
    # A raster seems to be indexed by [j,i] rather than [i,j]:
    for j in range(0, dims.ny):
      for i in range(0, dims.nx):
##        print(f"i {i}, j {j}, numFiresPopulationIsInside[{i},{j}] {numFiresPopulationIsInside[i,j]}")  # values are 'nan', -1, 0, 1, 2
#        print(f"j {j}, i {i}, numFiresPopulationIsInside[{j},{i}] {numFiresPopulationIsInside[j,i]}")  # values are -1, 0, 1, 2
#        print(f"numFiresPopulationIsInside[{j},{i}] is {numFiresPopulationIsInside[j,i]}")
##        if numFiresPopulationIsInside[i,j] > -1:
##        if numFiresPopulationIsInside[j,i] > -1:  # -1 indicates population is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to
#        if numFiresPopulationIsInside[j,i] > 0:  # -1 indicates population is zero in that raster-cell, and 0 indicates the cell is inside no fire-perimeters; in either case, no population-node need be found within that cell for the population to snap to
        if populationIsInsideFire[j,i] > 0:  # -1 indicates population is zero in that raster-cell, and 0 indicates the cell is not inside currentFire's perimeter; in either case, no population-node need be found within that cell for the population to snap to
##          print(f"i {i}, j {j}, numFiresPopulationIsInside[{i},{j}] {numFiresPopulationIsInside[i,j]}")  # values are 0, 1, 2
#          print(f"j {j}, i {i}, numFiresPopulationIsInside[{j},{i}] {numFiresPopulationIsInside[j,i]}")  # values are 0, 1, 2
          print(f"j {j}, i {i}, populationIsInsideFire[{j},{i}] {populationIsInsideFire[j,i]}")  # values are -1, 0, 1
          b = BoundingBox.from_list([ [dims.ox+dims.hx*i, dims.oy+dims.hy*j], [dims.ox+dims.hx*(i+1), dims.oy+dims.hy*(j+1)] ])  # create search box for this cell
#          b = BoundingBox.from_list([ [dims.oy+dims.hy*j, dims.ox+dims.hx*i], [dims.oy+dims.hy*(j+1), dims.ox+dims.hx*(i+1)] ])  # create search box for this cell
          network_currentCell = roadNetwork.region(b)  # find geometry within bounding-box
##          maxlinkcapacity_currentCell = maxlinkcapacity.currentCell(b)  # find raster within bounding-box: doesn't work, because (so far) maxlinkcapacity raster has incorrect values to begin with
          maxlinkcapacity_currentCell = network_currentCell.rasterise(CELLSIZEFINE_m, "output = max(0, capacity);", GeometryType.LineString)  # create finer-grained raster containing within each cell the maximum link-capacity: with CELLSIZEFINE_m = 80, this works
#          maxlinkcapacity_currentCell = network_currentCell.rasterise(640, "output = max(0, capacity);", GeometryType.LineString)  # this works, and slightly faster than with cell-size of 80, but the difference is only about 1%
##          maxlinkcapacity_currentCell = network_currentCell.rasterise(1280, "output = max(0, capacity);", GeometryType.LineString)  # cell-size of 1280 FAILS
          largest_maxlinkcapacity_currentCell = maxlinkcapacity_currentCell.max()
          print(f"largest_maxlinkcapacity_currentCell is {largest_maxlinkcapacity_currentCell}")
##          largestMaxOutCapacityFound_cell = NEGATIVEVALUE  # will record, over all nodes/Points in the raster-cell, the node's maximum out-capacity that is the largest
##          INITIALLOWERBOUND_CURRENTCELL = maxlinkcapacity[j,i] - 1.0
##          INITIALLOWERBOUND_CURRENTCELL = 600.0 - 1.0  # problem-dependent, unless there is a way to establish that all links have capacity of at least 600.0?
#          INITIALLOWERBOUND_CURRENTCELL = largest_maxlinkcapacity_currentCell - 1.0
#          alpha = INITIALLOWERBOUND_CURRENTCELL  # lower bound, over all nodes/Points in the raster-cell, on the node's maximum out-capacity that is the largest
          # (Search stops if finds a link with capacity at least beta; a lower beta means a quicker search, but risks not finding the highest capacity)
##          beta = maxlinkcapacity[j,i]  # upper bound, over all nodes/Points in the raster-cell, on the node's maximum out-capacity that is the largest
##          beta = 1500.0  # runs in 38 to 40 seconds
##          beta = 1000.0  # runs in 16 seconds
#          beta = largest_maxlinkcapacity_currentCell
#          print(f"alpha is {alpha}, beta is {beta}")
          print(f"len(network_currentCell.getPointIndexes()) is {len(network_currentCell.getPointIndexes())}")
          print(f"len(network_currentCell.getLineStringIndexes()) is {len(network_currentCell.getLineStringIndexes())}")
          largestMaxOutCapacityFound_cell = None
          idx_largestMaxOutCapacityFound_cell = None
          for idx in network_currentCell.getPointIndexes():  # set property on all Points within raster-cell; James Hilton: "One thing to note is that a property set on a derived Vector affects the parent Vector. So just setting the properties on each of the Vectors returned from region sets them on the parent Vector - these don't have to be copied back."
#            network_currentCell.setProperty(idx, "numFirePerimetersInside", numFiresPopulationIsInside[j,i] )
            maxOutCapacity_currentPoint = network_currentCell.getProperty(idx, "maxCapacityOut", float) 
#            if maxOutCapacity_currentPoint > alpha:
#              alpha = largestMaxOutCapacityFound_cell = maxOutCapacity_currentPoint
#              idx_largestMaxOutCapacityFound_cell = idx
#            if maxOutCapacity_currentPoint >= beta:
##              beta = maxOutCapacity_currentPoint
##              idx_largestMaxOutCapacityFound_cell = idx
#              break
#            elif maxOutCapacity_currentPoint > beta:
#              print(f"WARNING: maxOutCapacity_currentPoint {maxOutCapacity_currentPoint} exceeds beta {beta}.")
#              break
            if maxOutCapacity_currentPoint == largest_maxlinkcapacity_currentCell:
              largestMaxOutCapacityFound_cell = maxOutCapacity_currentPoint
              idx_largestMaxOutCapacityFound_cell = idx
              break
            elif maxOutCapacity_currentPoint > largest_maxlinkcapacity_currentCell:
              print(f"WARNING: maxOutCapacity_currentPoint {maxOutCapacity_currentPoint} exceeds largest_maxlinkcapacity_currentCell {largest_maxlinkcapacity_currentCell}.")
              break
          if idx_largestMaxOutCapacityFound_cell:
            MATsimID_node = network_currentCell.getProperty( idx_largestMaxOutCapacityFound_cell, 'matsim_nodeID' )
            print(f"population[{j},{i}] is {population[j,i]}")
            inflowsByNodeID[MATsimID_node] = population[j,i]
            positivePopulationInsideFireByNodeID[MATsimID_node] = population[j,i]
            inflowCoordsByNodeID[MATsimID_node] = np.array(network_currentCell.getPointCoordinate(idx_largestMaxOutCapacityFound_cell).to_list()[:2])  # take the first two coordinates only
#            print(f"inflowCoordsByNodeID[{MATsimID_node}] is {inflowCoordsByNodeID[MATsimID_node]}")
          else:
            MATsimID_node = None
          print(f"Property 'largestMaxOutCapacityFound_cell' is {largestMaxOutCapacityFound_cell}, at node with idx {idx_largestMaxOutCapacityFound_cell} and MATsimID_node {MATsimID_node}")
##          if largestMaxOutCapacityFound_cell == NEGATIVEVALUE:
#          if largestMaxOutCapacityFound_cell == INITIALLOWERBOUND_CURRENTCELL:
          if largestMaxOutCapacityFound_cell == None:
            print(f"WARNING: this raster-cell contains no node with a road-link leaving it.")
          network_currentCell.setProperty(idx_largestMaxOutCapacityFound_cell, "largestMaxOutCapacityInCell", largestMaxOutCapacityFound_cell)  # on the node inside the raster-cell that has the largest maximum out-capacity, set the property "largestMaxOutCapacityInCell" to this largest value
    print(f"inflowsByNodeID is {inflowsByNodeID}")
    print(f"len(inflowsByNodeID) is {len(inflowsByNodeID)}")
    print(f"positivePopulationInsideFireByNodeID is {positivePopulationInsideFireByNodeID}")
    with open(pth.join(outputDir, f"networkWithPropertiesByCell_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.geojson"), "w") as f:
      f.write(vectorToGeoJson(roadNetwork))
    elapsedTime = time()
    print(f"findEvacuationRisks.py has run for {elapsedTime - startTime} seconds so far.")
#    sys.exit()

    # Find currentFire's contour:
    print("Creating 'contourVector' layer ...")
    contourVector = isInsideFirePerimeter.vectorise( [ 0.5 ] )  # 0.5 is between the possible values, 0 and 1, of isInsideFirePerimeter
##    contourVector.write(pth.join(outputDir, f"contourVector_{int(fire_max_time):06d}.tif"))  # returns "'AttributeError: 'Vector' object has no attribute 'write'"
#    with open(pth.join(outputDir, f"contourVector_{int(fire_max_time):06d}.geojson"), "w") as contourfile:
#      contourfile.write(vectorToGeoJson(contourVector))

#    runScript('numFiresInside += currentFire < currentFire::fireMaxTime ? 1 : 0;', [numFiresInside, currentFire] )  # value of raster-cell is 1 if inside the fire-perimeter, 0 otherwise

    # Find distance from currentFire's contour:
    print("Creating 'dist' layer ...")
    dist = contourVector.mapDistance(resolution=180.0, geom_type=GeometryType.LineString, bounds=allFires_bounds)
    dist.name = 'dist'
    dist.setProjectionParameters(allFires_proj)
#    runScript('dist = currentFire < currentFire::fireMaxTime ? -dist : dist;', [dist, currentFire] )
    runScript('dist = currentFire > 0 ? -dist : dist;', [dist, currentFire] )  # ensure 'dist' is a signed-distance raster, where points inside the fire perimeter have a negative distance: this is needed so that exit-nodes will be selected only from the nodes that lie outside the perimeter and not inside it
#    dist.write(pth.join(outputDir, f"dist_{int(fire_max_time):06d}.tif"))  # write raster-layer 'dist' to a file
    dist.write(pth.join(outputDir, f"dist_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.tif"))  # write raster-layer 'dist' to a file

#    sys.exit()

#    # Create count layer:
#    print("Creating 'count' layer ...")
#    count = Raster(name="count", data_type=np.uint32)
##    count.init_with_bbox(allFires_bounds, CELLSIZESMALL_m)
#    count.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
#    count.setProjectionParameters(allFires_proj)
#    # Count points in cell:
#    count.setAllCellValues(0)
#    count.rasterise(roadNetwork, "atomic_inc();", GeometryType.Point)  # atomic_inc() increments the value in any cell containing vector-nodes from 'roadNetwork' layer
##    runScript('count = count == 0 ? nodata : count;', [count] )  # set each cell containing no vector-nodes to 'nodata'
#    runScript('count = count < 1 ? nodata : count;', [count] )  # set each cell containing no vector-nodes to 'nodata'
#    count.write(pth.join(outputDir, f"count_{int(fire_max_time):06d}.tif"))

##    # Create flow based on population:
##    print(f"Creating 'flowbasedonpopandcount_{int(fire_max_time):06d}.tif' layer ...")
##    runScript('''
##      flow = nodata;
##      if (dist < 10000.0 && count > 0) {
##        flow = population/(REAL)count;
##      }
##      ''', [flow, dist, count, population])
#    # Create injection-flows based on population:
#    print(f"Creating 'inflowbasedonpopandcount_{int(fire_max_time):06d}.tif' layer ...")
#    runScript('''
#      inflowPerNodeInsidePerimeter = nodata;
#      if (dist <= 0.0 && count > 0) {
#        inflowPerNodeInsidePerimeter = population/(REAL)count;
#      }
#      ''', [inflowPerNodeInsidePerimeter, dist, count, population])
#    inflowPerNodeInsidePerimeter.write(pth.join(outputDir, f"inflowbasedonpopandcount_{int(fire_max_time):06d}.tif"))


    # Create exit-nodes based on safe distance from currentFire's perimeter:
    atSafeDistFromPerimeter = Raster(name="atSafeDistFromPerimeter")
#    atSafeDistFromPerimeter.init_with_bbox(allFires_bounds, CELLSIZESMALL_m)
    atSafeDistFromPerimeter.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
    atSafeDistFromPerimeter.setProjectionParameters(allFires_proj)
#    print(f"Creating 'atSafeDistFromPerimeter_{int(fire_max_time):06d}.tif' layer ...")
    print(f"Creating 'atSafeDistFromPerimeter_ffdi100{fireName}.tif' layer ...")
#    runScript("atSafeDistFromPerimeter = dist >= 1000.0 && dist < 2000.0 ? 1 : 0;", [atSafeDistFromPerimeter, dist])
    dist.setVariableData("INNERDIST_SAFEBUFFER_m", INNERDIST_SAFEBUFFER_m)
    dist.setVariableData("OUTERDIST_SAFEBUFFER_m", OUTERDIST_SAFEBUFFER_m)
#    runScript("atSafeDistFromPerimeter = dist >= INNERDIST_SAFEBUFFER_m && dist < OUTERDIST_SAFEBUFFER_m ? 1 : 0;", [atSafeDistFromPerimeter, dist])
    runScript("atSafeDistFromPerimeter = dist >= dist::INNERDIST_SAFEBUFFER_m && dist < dist::OUTERDIST_SAFEBUFFER_m ? 1 : 0;", [atSafeDistFromPerimeter, dist])
##    atSafeDistFromPerimeter.write(pth.join(outputDir, f"atSafeDistFromPerimeter_{int(fire_max_time):06d}.tif"))
#    atSafeDistFromPerimeter.write(pth.join(outputDir, f"atSafeDistFromPerimeter__ffdi100{fireName}.tif"))
    atSafeDistFromPerimeter.write(pth.join(outputDir, f"atSafeDistFromPerimeter__{int(CELLSIZELARGE_m)}_ffdi100{fireName}.tif"))

    # Define the exit-nodes, at most one per raster-cell, as lying inside the "safe" buffer beyond currentFire's perimeter and being, for each raster-cell, the road-network node, if any, that has the largest maximum in-capacity (or more simply, the first node found inside the raster-cell):
    print("######## Defining the exit-nodes ...")
    dims = atSafeDistFromPerimeter.getDimensions()
    print(f"atSafeDistFromPerimeter dimensions is {dims}")
#    exitnodes = list()
    exitnodes = set()
#    EXITNODES_MUSTBECLOSETO_INJECTIONNODES = True  # remove each candidate exit-node that is not any injection-node's (equal-) closest candidate exit-node. The intention of this is to remove exit-nodes that are far from the injection-nodes. One problem with this, however, is that an exit-node on a highway might be removed simply because it's slightly farther from an injection-node than is another exit-node, even though the second exit-node is on a minor road and is therefore a less suitable exit-point.
    EXITNODES_MUSTBECLOSETO_INJECTIONNODES = False  # A problem with this: if the number of exit-nodes is not reduced in some way then for fire ffdi100a it takes about 59 seconds to find shortest paths from all injection-nodes to all exit-nodes. If shortest paths aren't found, though, maximum-flow will find some very long paths from injection-nodes to exit-nodes. TODO: a solution to this dilemma might be to allow an exit-node only if it lies on a high-capacity road-link, such as a highway.
#    candidateExitNodes = dict()
    candidateExitIndexByNodeID = dict()
    candidateExitNodeCoords = dict()
    # A raster seems to be indexed by [j,i] rather than [i,j]:
    for j in range(0, dims.ny):
      for i in range(0, dims.nx):
#        if atSafeDistFromPerimeter[j,i] > 0:
        if atSafeDistFromPerimeter[j,i] == 1:  # 0 indicates the raster-cell is outside the "safe" buffer, in which case no exit-node need be found within that cell
          print(f"j {j}, i {i}, atSafeDistFromPerimeter[{j},{i}] {atSafeDistFromPerimeter[j,i]}")  # values are 0, 1
          b = BoundingBox.from_list([ [dims.ox+dims.hx*i, dims.oy+dims.hy*j], [dims.ox+dims.hx*(i+1), dims.oy+dims.hy*(j+1)] ])  # create search box for this cell
          network_currentCell = roadNetwork.region(b)  # find geometry within bounding-box
          print(f"len(network_currentCell.getPointIndexes()) is {len(network_currentCell.getPointIndexes())}")
          print(f"len(network_currentCell.getLineStringIndexes()) is {len(network_currentCell.getLineStringIndexes())}")
#          largestMaxInCapacityFound_cell = None  # will record, over all nodes/Points in the raster-cell, the node's maximum in-capacity that is the largest
          idx_firstNodeFound_cell = None
          for idx in network_currentCell.getPointIndexes():  # set property on all Points within raster-cell; James Hilton: "One thing to note is that a property set on a derived Vector affects the parent Vector. So just setting the properties on each of the Vectors returned from region sets them on the parent Vector - these don't have to be copied back."
            idx_firstNodeFound_cell = idx
            break
          if idx_firstNodeFound_cell:
            MATsimID_node = network_currentCell.getProperty( idx_firstNodeFound_cell, 'matsim_nodeID' )
            if EXITNODES_MUSTBECLOSETO_INJECTIONNODES:
#              candidateExitNodes[MATsimID_node] = True
              pass
            else:
#              exitnodes.append(MATsimID_node)
              exitnodes.add(MATsimID_node)
          else:
            MATsimID_node = None
          if MATsimID_node:
            print(f"At node with idx {idx_firstNodeFound_cell}, MATsimID_node is {MATsimID_node}")
            if EXITNODES_MUSTBECLOSETO_INJECTIONNODES:
              candidateExitIndexByNodeID[MATsimID_node] = idx_firstNodeFound_cell
              candidateExitNodeCoords[idx_firstNodeFound_cell] = np.array(network_currentCell.getPointCoordinate(idx_firstNodeFound_cell).to_list()[:2])  # take the first two coordinates only
            else:
              network_currentCell.setProperty(idx_firstNodeFound_cell, "isExitNode", 1)  # on the first node found inside the raster-cell, set the property "isExitNode" to 1
          else:
            print(f"WARNING: this raster-cell contains no node of the road-network.")
    if EXITNODES_MUSTBECLOSETO_INJECTIONNODES:
      INFINITEDISTANCE = 9e9
      print(f"INFINITEDISTANCE is {INFINITEDISTANCE}.")
      # Eliminate each candidate exit-node that is not any injection-node's (equal-) closest candidate exit-node. (The intention of this is to remove exit-nodes that are far from the injection-nodes.)
      print(f"Candidate exit-nodes is {list(candidateExitIndexByNodeID.keys())}")
      print(f"len(candidateExitIndexByNodeID) is {len(candidateExitIndexByNodeID)}")
      for inflowNodeMATsimID in inflowCoordsByNodeID:
        injectionNodeCoords = inflowCoordsByNodeID[inflowNodeMATsimID]
        mindist_current = INFINITEDISTANCE
#        for candidateExitMATsimID in candidateExitNodes:
        for candidateExitMATsimID in candidateExitIndexByNodeID:
          idx_candidateExit = candidateExitIndexByNodeID[candidateExitMATsimID]
          candidateExitCoords = candidateExitNodeCoords[idx_candidateExit]
          dist_twonodes = scipy.linalg.norm(injectionNodeCoords - candidateExitCoords)
          assert dist_twonodes < INFINITEDISTANCE
          if dist_twonodes < mindist_current:
            mindist_current = dist_twonodes
            idx_closestCandidateExit = idx_candidateExit
            MATsimID_closestCandidateExit = candidateExitMATsimID
        network_currentCell.setProperty(idx_closestCandidateExit, "isExitNode", 1)  # on this injection-node's closest exit-node, set the property "isExitNode" to 1
        print(f"Adding to exitnodes the candidate exit-node {MATsimID_closestCandidateExit}, which is at distance {mindist_current} from injection-node {inflowNodeMATsimID}.")
        exitnodes.add(MATsimID_closestCandidateExit)
    exitnodes = list(exitnodes)  # convert from set to list
    print(f"exitnodes is {exitnodes}")
    print(f"len(exitnodes) is {len(exitnodes)}")
    exitNodesAtSafeDistfilename = pth.join(outputDir, f"exitNodesAtSafeDist_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.geojson")
#    with open(pth.join(outputDir, f"exitNodesAtSafeDist_{int(CELLSIZELARGE_m)}.geojson"), "w") as f:
    with open(exitNodesAtSafeDistfilename, "w") as f:
      f.write(vectorToGeoJson(roadNetwork))
#    networkBounds_clippedToAllFires = roadNetwork.getBounds()  # is in wrong projection
#    print(f"networkBounds_clippedToAllFires is {networkBounds_clippedToAllFires}")

    elapsedTime = time()
    print(f"findEvacuationRisks.py has run for {elapsedTime - startTime} seconds so far.")
#    sys.exit()  # for debugging, exit before running SEM

#    # Map population per node to network-node inflow:
##    roadNetwork.pointSample(inflowPerNodeInsidePerimeter)  # sample 'inflowPerNodeInsidePerimeter' raster at each point and write resulting value to the 'roadNetwork' vector-layer
#    # Map whether sufficiently far from currentFire's perimeter to network-node safety:
#    roadNetwork.pointSample(atSafeDistFromPerimeter)  # sample 'atSafeDistFromPerimeter' raster at each point and write resulting value to the 'roadNetwork' vector-layer
###    with open(pth.join(outputDir, f"networkNodesAtSafeDist_{int(fire_max_time):06d}.geojson"), "w") as networkfile:
##    with open(pth.join(outputDir, f"networkNodesAtSafeDist_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.geojson"), "w") as networkfile:
#    with open(pth.join(outputDir, f"networkNodesAtSafeDist_ffdi100{fireName}.geojson"), "w") as networkfile:
#      networkfile.write(vectorToGeoJson(roadNetwork))


    # Run flow-solver:
#    networkFlowSolver.run()
    SEMversion = 'SEM5'  # static-flow model based on maximum-flow method#, with each specified subflow from an injection-node to an exit-node being added to the arcs on that path
    print(f"SEMversion is {SEMversion}.")
##    inflowsByNodeID = {'matsimnode0': 700, 'matsimnode1': 500}
#    inflowsByNodeID = getInjectionNodesandInflows()  # injection-nodes must be inside currentFire's perimeter; for now, this is defined above by hand
    inflowsByNodeID = None  # not needed for SEM5
    inflowsandflowperiodsByNodeID = None  # not needed for SEM5

    positivePopulationInsideFireByNodeID = {'259608426': 6.0}  # one population-ode, for debugging only
##    exitnodes = {'matsimnode5', 'matsimnode8'}  # two exit-nodes, for debugging only
    exitnodes = {'1681513618'}  # one exit-node, for debugging only

    subflowsToAssignedExitNodesByInjectionNodeID = simulDurationInHours = None
    print(f"JSONnetworkfilename is {JSONnetworkfilename}")
#    SEMoutputGeoJSON = solversSEMversions.runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, subflowsToAssignedExitNodesByInjectionNodeID, simulDurationInHours)
    print(f"exitNodesAtSafeDistfilename is {exitNodesAtSafeDistfilename}")
#    (SEMoutputGeoJSON, criticalLinks) = solversSEMversions.runSEM( SEMversion, exitNodesAtSafeDistfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, positivePopulationInsideFireByNodeID, subflowsToAssignedExitNodesByInjectionNodeID, simulDurationInHours)  # run SEM version on the roadNetwork as clipped to 'allFires' bounds
    currentFireBounds = currentFire.getBounds()
    print(f"currentFireBounds is {currentFireBounds}.")
    currentFireBounds = currentFireBounds.convert(network_proj, allFires_proj)  # convert currentFireBounds to the roadNetwork projection
    print(f"currentFireBounds converted to roadNetwork projection is {currentFireBounds}.")
    (maxFlowSolnGeoJSON, criticalLinksInsideFireBoundingBox, evacuationTimeByNodeID) = solversSEMversions.runSEM( SEMversion, exitNodesAtSafeDistfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, positivePopulationInsideFireByNodeID, currentFireBounds, subflowsToAssignedExitNodesByInjectionNodeID, simulDurationInHours)  # run SEM version on the roadNetwork 'exitNodesAtSafeDistfilename', as this has been clipped to the 'allFires' bounds

#    print(f"Reading roadNetwork data from {maxFlowSolnGeoJSON} into a vector ...")
#    network_maxFlowSoln = geoJsonToVector(maxFlowSolnGeoJSON)
#    maxFlowSoln_proj = network_maxFlowSoln.getProjectionParameters()
#    # Restrict network_maxFlowSoln's bounds to a smaller area around currentFire, then convert these bounds to the network_maxFlowSoln projection:
#    maxFlowSoln_bounds = currentFire.getBounds()
#    print(f"maxFlowSoln_bounds is {maxFlowSoln_bounds}")
##    network_maxFlowSolnBounds0 = network_maxFlowSoln.getBounds()
##    print(f"network_maxFlowSoln.getBounds() is {network_maxFlowSolnBounds0}")
#    maxFlowSoln_bounds = maxFlowSoln_bounds.convert(maxFlowSoln_proj, allFires_proj)
#    print(f"maxFlowSoln_bounds after conversion to maxFlowSoln_proj is {maxFlowSoln_bounds}")
#    network_maxFlowSoln = network_maxFlowSoln.region(maxFlowSoln_bounds)  # clip network_maxFlowSoln to currentFire's bounds
##    network_maxFlowSolnBounds1 = network_maxFlowSoln.getBounds()
##    print(f"network_maxFlowSoln.getBounds() after clipping to network_maxFlowSoln.region(maxFlowSoln_bounds) is {network_maxFlowSolnBounds1}")
#    network_maxFlowSoln = network_maxFlowSoln.convert(allFires_proj)  # convert network_maxFlowSoln to same projection as 'currentFire' raster, for point-sampling
#    networkMaxFlowSolnfilename = pth.join(outputDir, f"networkMaxFlowSoln_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.geojson")
#    with open(networkMaxFlowSolnfilename, "w") as f:
#      f.write(vectorToGeoJson(network_maxFlowSoln))

#    # Mark a link as not critical if its endpoints' distances are more than OUTERDIST_SAFEBUFFER_m from currentFire's perimeter:
#    # Mark a link as not critical if both its endpoints fall outside currentFire's perimeter:
##    insideFireOrSafeBuffer = Raster(name="insideFireOrSafeBuffer")
##    insideFireOrSafeBuffer.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
##    insideFireOrSafeBuffer.setProjectionParameters(allFires_proj)
##    print(f"Creating 'insideFireOrSafeBuffer_ffdi100{fireName}.tif' layer ...")
##    dist.setVariableData("INNERDIST_SAFEBUFFER_m", INNERDIST_SAFEBUFFER_m)
##    dist.setVariableData("OUTERDIST_SAFEBUFFER_m", OUTERDIST_SAFEBUFFER_m)
###    runScript("insideFireOrSafeBuffer = dist >= dist::INNERDIST_SAFEBUFFER_m && dist < dist::OUTERDIST_SAFEBUFFER_m ? 1 : 0;", [insideFireOrSafeBuffer, dist])
##    runScript("insideFireOrSafeBuffer = dist < dist::OUTERDIST_SAFEBUFFER_m ? 1 : 0;", [insideFireOrSafeBuffer, dist])
##    insideFireOrSafeBuffer.write(pth.join(outputDir, f"insideFireOrSafeBuffer__{int(CELLSIZELARGE_m)}_ffdi100{fireName}.tif"))
#
#    print("######## Marking as not critical all links that fall entirely outside currentFire's perimeter ...")
#    dims = isInsideFirePerimeter.getDimensions()
#    print(f"isInsideFirePerimeter dimensions is {dims}")
#    criticalLinksOutsideFire = set()
#    # A raster seems to be indexed by [j,i] rather than [i,j]:
#    for j in range(0, dims.ny):
#      for i in range(0, dims.nx):
#        if isInsideFirePerimeter[j,i] == 0:  # 1 indicates the raster-cell is inside currentFire's perimeter
#          print(f"j {j}, i {i}, isInsideFirePerimeter[{j},{i}] {isInsideFirePerimeter[j,i]}")  # values are 0, 1
#          b = BoundingBox.from_list([ [dims.ox+dims.hx*i, dims.oy+dims.hy*j], [dims.ox+dims.hx*(i+1), dims.oy+dims.hy*(j+1)] ])  # create search box for this cell
#          network_currentCell = network_maxFlowSoln.region(b)  # find geometry within bounding-box
#          print(f"len(network_currentCell.getPointIndexes()) is {len(network_currentCell.getPointIndexes())}")
#          print(f"len(network_currentCell.getLineStringIndexes()) is {len(network_currentCell.getLineStringIndexes())}")
#          MATsimID_link = None
#          for idLS in network_currentCell.getLineStringIndexes():  # set property on all LineStrings within raster-cell; James Hilton: "One thing to note is that a property set on a derived Vector affects the parent Vector. So just setting the properties on each of the Vectors returned from region sets them on the parent Vector - these don't have to be copied back."
#            MATsimID_link = network_currentCell.getProperty( idLS, 'matsim_linkID' )
#            linkIsCritical = network_currentCell.getProperty( idLS, 'isCritical' )
#            print(f"At link with idLS {idLS}, MATsimID_link is {MATsimID_link}, linkIsCritical is {linkIsCritical}")
#            if linkIsCritical:
#              network_currentCell.setProperty(idLS, "isCritical", 0.0)  # on each link found inside the raster-cell, set the property "isCritical" to mark the link as being not critical
#              criticalLinksOutsideFire.add(MATsimID_link)
#          if not MATsimID_link:
#            print(f"WARNING: this raster-cell contains no link of the road-network.")
#    criticalLinksOutsideFire = list(criticalLinksOutsideFire)  # convert from set to list
#    print(f"criticalLinksOutsideFire is {criticalLinksOutsideFire}")
#    print(f"len(criticalLinksOutsideFire) is {len(criticalLinksOutsideFire)}")
#    networkWithCriticalLinksOnlyInsideFirefilename = pth.join(outputDir, f"networkWithCriticalLinksOnlyInsideFire_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.geojson")
#    with open(networkWithCriticalLinksOnlyInsideFirefilename, "w") as f:
#      f.write(vectorToGeoJson(network_maxFlowSoln))


    print(f"len(evacuationTimeByNodeID) is {len(evacuationTimeByNodeID)}")
    for (popNode, pointID) in evacuationTimeByNodeID:
      if (popNode, pointID) in meanTimeToEvacuateByPopulationNode:
#        sumEvacuationTimesSoFar = meanTimeToEvacuateByPopulationNode[popNode]['sumEvacTimes']
#        numFiresInsideSoFar = meanTimeToEvacuateByPopulationNode[popNode]['numFiresInside']
        meanTimeToEvacuateByPopulationNode[popNode, pointID]['sumEvacTimes'] += evacuationTimeByNodeID[ popNode, pointID ]
        meanTimeToEvacuateByPopulationNode[popNode, pointID]['numFiresInside'] += 1
      else:
        meanTimeToEvacuateByPopulationNode[popNode, pointID] = dict()
        meanTimeToEvacuateByPopulationNode[popNode, pointID]['sumEvacTimes'] = evacuationTimeByNodeID[ popNode, pointID ]
        meanTimeToEvacuateByPopulationNode[popNode, pointID]['numFiresInside'] = 1
    print(f"len(positivePopulationInsideFireByNodeID) is {len(positivePopulationInsideFireByNodeID)}")
    numPopulationNodesByFire[ fireName ] = len(positivePopulationInsideFireByNodeID)
##    print(f"criticalLinks is {criticalLinks}")
##    print(f"len(criticalLinks) is {len(criticalLinks)}")
#    print(f"criticalLinksInsideFireBoundingBox is {criticalLinksInsideFireBoundingBox}")
    print(f"len(criticalLinksInsideFireBoundingBox) is {len(criticalLinksInsideFireBoundingBox)}")
    numCriticalLinksByFire[ fireName ] = len(criticalLinksInsideFireBoundingBox)
    for (MATsimID_link, lineStringID) in criticalLinksInsideFireBoundingBox:
      if (MATsimID_link, lineStringID) in numEvacuationsLinkIsCriticalTo:
        numEvacuationsLinkIsCriticalTo[ MATsimID_link, lineStringID ] += 1
      else:
        numEvacuationsLinkIsCriticalTo[ MATsimID_link, lineStringID ] = 1


    # Write data:
#    outName, outExt = pth.splitext(outputGeoJSONfilename)
#    print(f"outName is {outName}, outExt is {outExt}")
##    outVector = networkFlowSolver.getNetwork().convert(network_proj)
#    with open(f"{outName}_{int(fire_max_time):06d}_{outExt}", "w") as outfile:

    assert outputGeoJSONfilename[-8:] == ".geojson"
#    with open(outputGeoJSONfilename, "w") as outfile:
    with open(outputGeoJSONfilename[:-8] + f"_{int(CELLSIZELARGE_m)}_ffdi100{fireName}.geojson", "w") as outfile:
###      outfile.write(vectorToGeoJson(outVector))
##      outfile.write(SEMoutputGeoJSON)
#      json.dump(SEMoutputGeoJSON, outfile)
      json.dump(maxFlowSolnGeoJSON, outfile)


  MAXEVACTIME_FINITE = 100.0  # (in hours) maximum time that any population-node should take to evacuate if it isn't infinite
  for (popNode, pointID) in meanTimeToEvacuateByPopulationNode:
    meanTimeToEvacuateByPopulationNode[popNode, pointID]['meanTimeToEvacuate'] = meanTimeToEvacuateByPopulationNode[popNode, pointID]['sumEvacTimes'] / meanTimeToEvacuateByPopulationNode[popNode, pointID]['numFiresInside']
    if meanTimeToEvacuateByPopulationNode[popNode, pointID]['meanTimeToEvacuate'] > MAXEVACTIME_FINITE:
      meanTimeToEvacuateByPopulationNode[popNode, pointID]['meanTimeToEvacuate'] = MAXEVACTIME_FINITE
  print(f"meanTimeToEvacuateByPopulationNode is {meanTimeToEvacuateByPopulationNode}")
  print(f"numPopulationNodesByFire is {numPopulationNodesByFire}")
  print(f"numCriticalLinksByFire is {numCriticalLinksByFire}")
  print(f"numEvacuationsLinkIsCriticalTo is {numEvacuationsLinkIsCriticalTo}")
  print(f"len(numEvacuationsLinkIsCriticalTo) is {len(numEvacuationsLinkIsCriticalTo)}")


  print("Creating 'ignitionRisk' layer ...")  # to record the ignition-risk of each raster-cell, i.e. the total number of population-nodes within all fires the cell is inside plus the total number of critical links for all fires whose perimeters the cell is inside
  ignitionRisk = Raster(name="ignitionRisk")
  ignitionRisk.init_with_bbox(allFires_bounds, CELLSIZELARGE_m)
  ignitionRisk.setProjectionParameters(allFires_proj)
  ignitionRisk.setAllCellValues(0)
  # Read all fire-layers to create 'ignitionRisk' layer:
#  for fireName in allFireNames:
  for fireName in evacuationScenarios:
    evacuationLayer = f'modelinputsMtAlexander_100{fireName}/sem/20181109_mountalex_evac_ffdi100{fireName}_grid.tif'
    print(f"Reading fire data from {evacuationLayer}...")
    if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:  # use native raster-readers:
      currentFire = Raster(name=f"currentFire")
      currentFire.read(evacuationLayer)
    else:  # use gdal based raster reader:
      currentFire = RasterFile(name = f"currentFire", filePath = evacuationLayer, backend = 'gdal')
      currentFire.read()
    currentFire.setVariableData("ignitionRisk", numPopulationNodesByFire[fireName] + numCriticalLinksByFire[fireName])
#    runScript('numFiresInside += currentFire < currentFire::fireMaxTime ? 1 : 0;', [numFiresInside, currentFire] )  # add 1 to value of cell if inside the fire-perimeter, 0 otherwise
    runScript('ignitionRisk += currentFire > 0 ? currentFire::ignitionRisk : 0;', [ignitionRisk, currentFire] )  # add 1 to value of cell if inside the fire-perimeter, 0 otherwise
#  dims = ignitionRisk.getDimensions()
#  print(f"ignitionRisk dimensions is {dims}")
  ignitionRisk.write(pth.join(outputDir, f"ignitionRisk_{int(CELLSIZELARGE_m)}.tif"))  # write raster-layer 'ignitionRisk' to a file

  elapsedTime = time()
  print(f"findEvacuationRisks.py has run for {elapsedTime - startTime} seconds so far.")


  # Record each road-link's point-of-impact risk:
  print("Creating 'roadNetworkWithRiskProperties' layer ...")  # to record the point-of-impact risk for each road-link, i.e. the number of evacuation-flows the link is critical to, as well as the community risk for each population-node, i.e. the mean time taken to evacuate each node over all fires that contain it
  print("Setting each road-link's 'numEvacsLinkCriticalTo' property ...")
#  numLinksProcessed = 0  # for monitoring progress of this time-consuming computation
#  for idLS in roadNetwork.getLineStringIndexes():
#    MATsimID_link = roadNetwork.getProperty( idLS, 'matsim_linkID' )
#    if MATsimID_link in numEvacuationsLinkIsCriticalTo:
###      print(f"idLS {idLS}, MATsimID_link {MATsimID_link}")
##      roadNetwork.setProperty(idLS, "numEvacsLinkCriticalTo_2", numEvacuationsLinkIsCriticalTo[MATsimID_link])
#      roadNetwork.setProperty(idLS, "numEvacsLinkCriticalTo", numEvacuationsLinkIsCriticalTo[MATsimID_link])
#      numLinksProcessed += 1
#      if numLinksProcessed % 100 == 0:
#        print(f"numLinksProcessed is {numLinksProcessed}")
#    else:
##      print(f"idLS {idLS}: MATsimID_link {MATsimID_link} NOT in numEvacuationsLinkIsCriticalTo.")
#      pass
##    roadNetwork.setProperty(idLS, "numEvacsLinkCriticalTo", numEvacuationsLinkIsCriticalTo.get(MATsimID_link,0))
  for (MATsimID_link, lineStringID) in numEvacuationsLinkIsCriticalTo:
    roadNetwork.setProperty(lineStringID, "numEvacsLinkCriticalTo", numEvacuationsLinkIsCriticalTo[MATsimID_link, lineStringID])

  elapsedTime = time()
  print(f"findEvacuationRisks.py has run for {elapsedTime - startTime} seconds so far.")

  # Record each injection-node's community risk:
  print("Setting each injection-node's 'meanTimeToEvacuate' property ...")
#  for idx in roadNetwork.getPointIndexes():
#    MATsimID_node = roadNetwork.getProperty( idx, 'matsim_nodeID' )
#    if MATsimID_node in meanTimeToEvacuateByPopulationNode:
##      roadNetwork.setProperty(idx, "meanTimeToEvacuate_2", meanTimeToEvacuateByPopulationNode[MATsimID_node]['meanTimeToEvacuate'])
#      roadNetwork.setProperty(idx, "meanTimeToEvacuate", meanTimeToEvacuateByPopulationNode[MATsimID_node]['meanTimeToEvacuate'])
  for (MATsimID_node, pointID) in meanTimeToEvacuateByPopulationNode:
#    print(f"MATsimID_node is {repr(MATsimID_node)}, pointID is {repr(pointID)}")
    roadNetwork.setProperty(pointID, "meanTimeToEvacuate", meanTimeToEvacuateByPopulationNode[MATsimID_node, pointID]['meanTimeToEvacuate'])


  elapsedTime = time()
  print(f"findEvacuationRisks.py has run for {elapsedTime - startTime} seconds so far.")

  with open(pth.join(outputDir, f"roadNetworkWithRiskProperties_{int(CELLSIZELARGE_m)}_allFires.geojson"), "w") as f:
    f.write(vectorToGeoJson(roadNetwork))


  endTime = time()
  print(f"findEvacuationRisks.py ran for {endTime - startTime} seconds in total.")
  return 0


if __name__ == "__main__":

  # Parse command-line arguments:
  parser = ArgumentParser(description="Calculate and plot evacuation risks")
  group = parser.add_mutually_exclusive_group()
  group.add_argument("--configJSON", dest="configJSON", type=str, help="Static evacuation json configuration")
  group.add_argument("--configFile", dest="configFile", type=str, help="Path to static evacuation json configuration file")
  args = parser.parse_args()

  if args.configJSON != None:
    try:
      cfgJson = json.loads(args.configJSON)  # open configuration JSON-file
    except ValueError as e:
      raise RuntimeError("JSON configuration-file is not valid")
    # Find risks using configuration-file:
    if findRisks(cfgJson) == 0:
      raise RuntimeError("Unable to run static evacuation model")

  else:
    if args.configFile is None:
      # Error if no config-file provided:
      raise RuntimeError("Spark Configuration file is not provided")
    elif not pth.exists(args.configFile) or not pth.isfile(args.configFile):
      # Error if file cannot be found:
      raise FileNotFoundError(f"Configuration file '{args.configFile}' not found")
    else:
      with open(args.configFile, 'r') as inp:  # open configuration JSON-file
        cfgJson = json.load(inp)  # parse JSON-file
        # Set working directory:
        workingDir = pth.dirname(pth.abspath(args.configFile))
        os.chdir(workingDir)
        print(f"workingDir is {workingDir}")
        findRisks(cfgJson)  # find risks using configuration-file
