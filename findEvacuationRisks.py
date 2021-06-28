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

  RUNMAXFLOW = False

#  print("===Geostack evacuation-risks identifier===")
  startTime = time()

  # Check json entries
  try:
#    evacuationLayer = getParam(cfgJson, "evacuationLayer", (str,))
    populationLayer = getParam(cfgJson, "populationLayer", (str,))
#    networkGeoJSON = getParam(cfgJson, "networkGeoJSON", (str,))
    JSONnetworkfilename = getParam(cfgJson, "networkGeoJSON", (str,))
    outputGeoJSONfilename = getParam(cfgJson, "outputGeoJSON", (str,))
    networkConfig = getParam(cfgJson, "networkConfig", (dict,))

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

  CELLSIZEFINE = 80.0
##  CELLSIZESMALL = 180.0
##  CELLSIZELARGE = 1800.0
  CELLSIZESMALL = 100.0

  CELLSIZELARGE = 1000.0
#  CELLSIZELARGE = 2000.0
#  CELLSIZELARGE = 10000.0
  print(f"CELLSIZEFINE is {CELLSIZEFINE}, CELLSIZESMALL is {CELLSIZESMALL}, CELLSIZELARGE is {CELLSIZELARGE}")

#  rasterCellSize_small = CELLSIZESMALL
#  rasterCellSize_large = CELLSIZELARGE

  # Read just one fire-layer to obtain its projection, and extend its bounds to define bounds encompassing all fires:
#  evacuationLayer = 'modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif'
  evacuationLayer = 'modelinputsMtAlexander_100c/sem/20181109_mountalex_evac_ffdi100c_grid.tif'
  print(f"Reading fire data from {evacuationLayer}...")
  if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:  # use native raster-readers
    fire = Raster(name="fire")
    fire.read(evacuationLayer)
  else:  # use gdal-based raster-reader
    fire = RasterFile(name = "fire", filePath = evacuationLayer, backend = 'gdal')
    fire.read()
  allFires_proj = fire.getProjectionParameters()
#  print(f"allFires_proj is {allFires_proj}")
  allFires_bounds = fire.getBounds()
  allFires_bounds.extend(20000.0)  # extend bounds of the fire layer by 20 km; TODO: the purpose of this extension is to cover all fires; in general, it might be necessary to calculate a union of the fire-rasters

  # Read network:
  print(f"Reading network data from {JSONnetworkfilename} ...")
  network = geoJsonToVector(JSONnetworkfilename)
  network_proj = network.getProjectionParameters()
  # Restrict the network's bounds to a smaller area around the fire, then convert these bounds to the network projection:
#  network_bounds = fire.getBounds()
##  network_bounds = allFires_bounds  # bug: allFires_bounds will be changed by 'extend' in the next line
#  network_bounds.extend(30000.0)
  network_bounds = allFires_bounds
#  network_bounds.extend(10000.0)  # extend bounds of the layer by 10 km
  print(f"allFires_bounds is {allFires_bounds}")
  print(f"network_bounds is {network_bounds}")
  assert network_bounds == allFires_bounds

  networkBounds0 = network.getBounds()
  print(f"network.getBounds() is {networkBounds0}")

  network_bounds = network_bounds.convert(network_proj, allFires_proj)
  print(f"network_bounds after conversion to network_proj is {network_bounds}")

  network = network.region(network_bounds)  # clip network to the bounds
  networkBounds1 = network.getBounds()
  print(f"network.getBounds() after clipping to network.region(network_bounds) is {networkBounds1}")

  network = network.convert(allFires_proj)  # convert network to same projection as the 'fire' raster, for point-sampling
  networkBounds2 = network.getBounds()
  print(f"network.getBounds() after conversion to allFires_proj is {networkBounds2}")  # network bounds are now wide again - a bug in Geostack?

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
#  numFiresInside.init_with_bbox(population_unaligned_bounds, CELLSIZESMALL)
#  numFiresInside.init_with_bbox(network_bounds, CELLSIZESMALL)
#  numFiresInside.init_with_bbox(ifp_bounds, CELLSIZESMALL)
  numFiresInside.init_with_bbox(allFires_bounds, CELLSIZELARGE)
  numFiresInside.setProjectionParameters(allFires_proj)
  numFiresInside.setAllCellValues(0)

  # Read all fire-layers to create 'numFiresInside' layer:
##  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif')
#  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif', 'modelinputsMtAlexander_100c/sem/20181109_mountalex_evac_ffdi100c_grid.tif', 'modelinputsMtAlexander_100d/sem/20181109_mountalex_evac_ffdi100d_grid.tif')
  evacuationScenarios = ('a', 'b', 'c', 'd')  # fires ffdi100a, ffdi100b, ...
#  evacuationLayers = [f'modelinputsMtAlexander_100{es}/sem/20181109_mountalex_evac_ffdi100{es}_grid.tif' for es in evacuationScenarios]
#  print(f"evacuationLayers is {evacuationLayers}")
#  for evacuationLayer in evacuationLayers:
  for fireName in evacuationScenarios:
    evacuationLayer = f'modelinputsMtAlexander_100{fireName}/sem/20181109_mountalex_evac_ffdi100{fireName}_grid.tif'
    print(f"Reading fire data from {evacuationLayer}...")
    if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:
      # use native raster-readers:
#      fire = Raster(name=f"fire_{evacuationLayer}")
      fire = Raster(name=f"fire")
      fire.read(evacuationLayer)
    else:
      # use gdal based raster reader:
#      fire = RasterFile(name = f"fire_{evacuationLayer}", filePath = evacuationLayer, backend = 'gdal')
      fire = RasterFile(name = f"fire", filePath = evacuationLayer, backend = 'gdal')
      fire.read()
#    fire_proj = fire.getProjectionParameters()
#    print(f"fire_proj is {fire_proj}")
    fire_max_time = np.nanmax(fire.data)  # time that defines perimeter of fire's maximum extent
    fire.setVariableData("time", fire_max_time)  # set time in 'fire' layer
    runScript('numFiresInside += fire < fire::time ? 1 : 0;', [numFiresInside, fire] )  # add 1 to value of cell if inside the fire-perimeter, 0 otherwise

  dims = numFiresInside.getDimensions()
  print(f"numFiresInside dimensions is {dims}")
  numFiresInside.write(pth.join(outputDir, f"numFiresInside_{int(CELLSIZELARGE)}.tif"))  # write raster-layer 'numFiresInside' to a file

#  numFiresPopulationIsInside = Raster(name="numFiresPopulationIsInside", data_type=np.int32)
  numFiresPopulationIsInside = Raster(name="numFiresPopulationIsInside")
  numFiresPopulationIsInside.init_with_bbox(allFires_bounds, CELLSIZELARGE)
  numFiresPopulationIsInside.setProjectionParameters(allFires_proj)
#  runScript("numFiresPopulationIsInside = population_unaligned > 0 ? numFiresInside : noData_UINT;", [numFiresPopulationIsInside, population_unaligned, numFiresInside])  # results in very large raster-cell values, for some reason
  runScript("numFiresPopulationIsInside = population_unaligned > 0 ? numFiresInside : -1;", [numFiresPopulationIsInside, population_unaligned, numFiresInside])  # -1 indicates population_unaligned is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to; numFiresPopulationIsInside is 'nan' where population_unaligned > 0 and numFiresInside is 'nan'


# Using lines like the following with 'output' seems to result in 'numFiresPopulationIsInside' raster having a number of bands equal to 'CELLSIZELARGE' - why?
###  numFiresPopulationIsInside = runScript("output = population_unaligned * numFiresInside;", [population_unaligned, numFiresInside])
###  numFiresPopulationIsInside = runScript("output = population_unaligned && numFiresInside;", [population_unaligned, numFiresInside])
##  numFiresPopulationIsInside = runScript("output = population_unaligned > 0 ? numFiresInside : 0;", [population_unaligned, numFiresInside])
#  numFiresPopulationIsInside = runScript("output = population_unaligned > 0 ? numFiresInside : -1;", [population_unaligned, numFiresInside])  # -1 indicates population_unaligned is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to

  dims_numFiresPopulationInside = numFiresPopulationIsInside.getDimensions()
#  print(f"dims_numFiresPopulationInside is {dims_numFiresPopulationInside}")
  numFiresPopulationIsInside.write(pth.join(outputDir, f"numFiresPopulationIsInside_{int(CELLSIZELARGE)}.tif"))  # write raster 'numFiresPopulationIsInside' to a file

  population = runScript("output = numFiresPopulationIsInside;", [numFiresPopulationIsInside])  # make new raster with same dimensions as numFiresPopulationIsInside, to store population values; copies the contents of numFiresPopulationIsInside rather than simply making a link to the raster, so that runScript on population will NOT change the contents of numFiresPopulationIsInside
  population.name = 'population'
  runScript("population = population_unaligned;", [population, population_unaligned])  # set population to a raster aligned with numFiresPopulationIsInside
  dims_population = population.getDimensions()
  print(f"dims_population is {dims_population}")
  assert dims_population == dims_numFiresPopulationInside
  population.write(pth.join(outputDir, f"populationArchetypes_{int(CELLSIZELARGE)}.tif"))  # write 'population' to a file

#  sys.exit()

#  network.pointSample(numFiresInside)  # sample 'numFiresInside' raster at each point and write resulting value to the 'network' vector-layer
#  with open(pth.join(outputDir, f"clippedNetwork_{int(CELLSIZELARGE)}.geojson"), "w") as f:
#    f.write(vectorToGeoJson(network))


  # Read fire arrival-time input-layers:
##  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', )  # just one fire, for simplicity
##  evacuationLayers = ('modelinputsMtAlexander_100a/sem/20181109_mountalex_evac_ffdi100a_grid.tif', 'modelinputsMtAlexander_100b/sem/20181109_mountalex_evac_ffdi100b_grid.tif')
##  for evacuationLayer in evacuationLayers:
#  evacuationScenarios = ('a', 'b', 'c', 'd')  # fires ffdi100a, ffdi100b, ...
  evacuationScenarios = ('a', )  # the single fire ffdi100a
  for fireName in evacuationScenarios:
    evacuationLayer = f'modelinputsMtAlexander_100{fireName}/sem/20181109_mountalex_evac_ffdi100{fireName}_grid.tif'
    print(f"Reading fire data from {evacuationLayer}...")
    if pth.splitext(evacuationLayer)[-1].lower() in ['.tif', '.asc', '.gsr', '.flt']:
      # use native raster-readers:
      fire = Raster(name="fire")
      fire.read(evacuationLayer)
    else:
      # use gdal-based raster-reader:
      fire = RasterFile(name = "fire", filePath = evacuationLayer, backend = 'gdal')
      fire.read()
#    fire_proj = fire.getProjectionParameters()

#    # Get fire-arrival-time bounds:
#    fire_bounds = fire.getBounds()
##    fire_bounds.extend(10000.0)  # extend bounds of the fire layer by 10 km

    # Create initialFireExtent layer, to hold the area of the fire's epicentre:
#    print("Creating 'initialFireExtent' layer ...")
    initialFireExtent = Raster(name="initialFireExtent")
    initialFireExtent.init_with_bbox(allFires_bounds, CELLSIZESMALL)
    initialFireExtent.setProjectionParameters(allFires_proj)

    fire_time_epicentre = 1.0  # time in seconds used to define fire's epicentre
    print(f"Finding fire's epicentre by running for time (fire-duration) {fire_time_epicentre}...")
    fire.setVariableData("time", fire_time_epicentre)  # set current time in 'fire' layer
    runScript("initialFireExtent = fire < fire::time ? 1 : 0;", [initialFireExtent, fire])
#    initialFireExtent.write(pth.join(outputDir, f"firePerimeter_{int(fire_time_epicentre):06d}.tif"))

##    print("Creating 'contourOfInitialFireExtent' layer ...")
#    contourOfInitialFireExtent = initialFireExtent.vectorise( [ 0.5 ] )
#    with open(pth.join(outputDir, f"contourOfInitialFireExtent_{int(fire_time_epicentre):06d}.geojson"), "w") as contourfile:
#      contourfile.write(vectorToGeoJson(contourOfInitialFireExtent))

#    perimeterBounds = contourOfInitialFireExtent.getBounds()
#    print(f"perimeterBounds[0] is {repr(perimeterBounds[0])}")  # returns nothing useful
##    bb = BoundingBox.from_list([[0.5, 0.5], [0.5, 0.5]])  # create search-point
##    nearestNodes = network.nearest(bb)  # find node nearest to epicentre
#    nearestNodes = network.nearest(perimeterBounds)  # find node nearest to epicentre
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
#      print(f"Nearest node to fire epicentre: nodeIndex is {nodeIndex}, nodeCoordinates is {nodeCoordinates}")
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


    fire_max_time = np.nanmax(fire.data)  # time that defines perimeter of fire's maximum extent
#    print(f"Running for time (fire-duration) {fire_max_time}...")
    print(f"fire_max_time is {fire_max_time}")

#    # Create inflowPerNodeInsidePerimeter layer:
##    print("Creating 'inflowPerNodeInsidePerimeter' layer ...")
#    inflowPerNodeInsidePerimeter = Raster(name="inflowPerNodeInsidePerimeter")
#    inflowPerNodeInsidePerimeter.init_with_bbox(allFires_bounds, CELLSIZESMALL)
#    inflowPerNodeInsidePerimeter.setProjectionParameters(allFires_proj)
#    fire.setVariableData("time", fire_max_time)  # set current time in 'fire' layer
#    runScript("inflowPerNodeInsidePerimeter = fire < fire::time ? 1 : 0;", [inflowPerNodeInsidePerimeter, fire])
##    inflowPerNodeInsidePerimeter.write(pth.join(outputDir, f"inflowPerNodeInsidePerimeter_{int(fire_max_time):06d}.tif"))

    # Create isInsideFirePerimeter layer:
#    print("Creating 'isInsideFirePerimeter' layer ...")
    isInsideFirePerimeter = Raster(name="isInsideFirePerimeter")
#    isInsideFirePerimeter.init_with_bbox(allFires_bounds, CELLSIZESMALL)
    isInsideFirePerimeter.init_with_bbox(allFires_bounds, CELLSIZELARGE)
    isInsideFirePerimeter.setProjectionParameters(allFires_proj)
    fire.setVariableData("time", fire_max_time)  # set current time in 'fire' layer
    runScript("isInsideFirePerimeter = fire < fire::time ? 1 : 0;", [isInsideFirePerimeter, fire])
    isInsideFirePerimeter.write(pth.join(outputDir, f"isInsideFirePerimeter_{int(CELLSIZELARGE)}_ffdi100{fireName}.tif"))  # write raster 'isInsideFirePerimeter' to a file

#    print("Creating 'maxlinkcapacity_unaligned' layer ...")  # get maximum link-capacity within each cell
###    maxlinkcapacity_unaligned = Raster(name="maxlinkcapacity_unaligned")
##    maxlinkcapacity_unaligned = Raster(name="maxlinkcapacity_unaligned", data_type=np.uint32)
##    maxlinkcapacity_unaligned.init_with_bbox(allFires_bounds, CELLSIZELARGE)
##    maxlinkcapacity_unaligned.setProjectionParameters(allFires_proj)
##    maxlinkcapacity_unaligned.setAllCellValues(0)
#    # James Hilton, 17th June 2021:
#    # > The atomic_max function only works for integers in OpenCL (I have no idea why). If it did work I think you could just do: r.rasterise(v, "atomic_max(capacity);") where r is your count raster, v the vector and capacity the field name. For this capacity would have to be an integer. I was thinking you could do this with a vector script but I think the same issue will occur.
##    maxlinkcapacity_unaligned.rasterise(network, "atomic_max('capacity');", GeometryType.LineString)  # find maximum link-capacity within each cell; doesn't work, for some reason, and runs even if 'capacity' is replaced with a typo such as 'cpacity'
#    # 'network' vector has much wider bounds than 'numFiresInside', so clip 'network' to the right size:
#    network = network.region(allFires_bounds)  # clip network again, but to allFires_bounds, which are now in the same projection; this makes the bounds similar to those in allFires_bounds, so that rasterising the network to create maxlinkcapacity_unaligned will give the raster-layer suitable bounds
#    networkBounds3 = network.getBounds()
#    print(f"network.getBounds() after clipping to network.region(allFires_bounds) is {networkBounds3}")
##    maxlinkcapacity_unaligned = network.rasterise(CELLSIZELARGE, "output = max(output, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = network.rasterise(CELLSIZELARGE, "output = max(output, capacity);")  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = network.rasterise(CELLSIZELARGE, "output = max(capacity, output);", 2)  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = network.rasterise(CELLSIZELARGE, "output = min(5000, capacity);", GeometryType.LineString)  # gives different cell-values each time it's run
##    maxlinkcapacity_unaligned = network.rasterise(CELLSIZELARGE, "output = max(0, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell: gives different cell-values each time it's run; this line doesn't seem to work unless the raster-cell size is small (it works certainly for a cell-size between 80 and 640)
#    maxlinkcapacity_unaligned = network.rasterise(80, "output = max(0, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell: this works, at its fine-grained scale of 80
##    maxlinkcapacity_unaligned = network.rasterise(640, "output = max(0, capacity);", GeometryType.LineString)  # create raster containing maximum link-capacity within each cell; as the scale is increased this line works less well: the higher-capacity links get "overwritten" by lower-capacity ones
#    maxlinkcapacity_unaligned.write(pth.join(outputDir, f"maxlinkcapacity_unaligned_{int(CELLSIZELARGE)}.tif"))
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
#    maxlinkcapacity.write(pth.join(outputDir, f"maxlinkcapacity_{int(CELLSIZELARGE)}.tif"))  # write 'maxlinkcapacity', actually containing maximum link-capacity per cell, to a file
#    print(f"maxlinkcapacity.max() is {maxlinkcapacity.max()}")
##    sys.exit()
#
#    maxlinkcapacity_proj = maxlinkcapacity.getProjectionParameters()
#    assert maxlinkcapacity_proj == allFires_proj
#
###    network_bounds = network_bounds.convert(network_proj, allFires_proj)
###    maxlinkcapacity = maxlinkcapacity.region(network_bounds)  # clip raster to the bounds: doesn't work, because seems to re-initialise the raster
##    maxlinkcapacity = maxlinkcapacity.region(allFires_bounds)  # clip raster to the bounds: doesn't work, because seems to re-initialise the raster
###    maxlinkcapacity = maxlinkcapacity.convert(allFires_proj)  # convert maxlinkcapacity to same projection as the 'fire' raster; doesn't work, as "'Raster' object has no attribute 'convert'"
##    dims = maxlinkcapacity.getDimensions()
##    print(f"maxlinkcapacity dimensions is {dims}")
##    maxlinkcapacity.write(pth.join(outputDir, f"maxlinkcapacity_{int(CELLSIZELARGE)}.tif"))

    populationIsInsideFire = Raster(name="populationIsInsideFire")
    populationIsInsideFire.init_with_bbox(allFires_bounds, CELLSIZELARGE)
    populationIsInsideFire.setProjectionParameters(allFires_proj)
#    runScript("populationIsInsideFire = population_unaligned > 0 ? isInsideFirePerimeter : -1;", [populationIsInsideFire, population_unaligned, isInsideFirePerimeter])  # -1 indicates population_unaligned is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to; populationIsInsideFire is 'nan' where population_unaligned > 0 and isInsideFirePerimeter is 'nan'
    runScript("populationIsInsideFire = population > 0 ? isInsideFirePerimeter : -1;", [populationIsInsideFire, population, isInsideFirePerimeter])  # -1 indicates population is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to; populationIsInsideFire is 'nan' where population > 0 and isInsideFirePerimeter is 'nan'
    populationIsInsideFire.write(pth.join(outputDir, f"populationIsInsideFire_{int(CELLSIZELARGE)}_ffdi100{fireName}.tif"))  # write raster 'populationIsInsideFire' to a file

##    network.addProperty("numFirePerimetersInside")
#    dims = numFiresPopulationIsInside.getDimensions()
#    print(f"numFiresPopulationIsInside dimensions is {dims}")
    dims = populationIsInsideFire.getDimensions()
    print(f"populationIsInsideFire dimensions is {dims}")
#    NEGATIVEVALUE = -99.9
    inflowsByNodeID = dict()
    # A raster seems to be indexed by [j,i] rather than [i,j]:
    for j in range(0, dims.ny):
      for i in range(0, dims.nx):
##        print(f"i {i}, j {j}, numFiresPopulationIsInside[{i},{j}] {numFiresPopulationIsInside[i,j]}")  # values are 'nan', -1, 0, 1, 2
#        print(f"j {j}, i {i}, numFiresPopulationIsInside[{j},{i}] {numFiresPopulationIsInside[j,i]}")  # values are -1, 0, 1, 2
#        print(f"numFiresPopulationIsInside[{j},{i}] is {numFiresPopulationIsInside[j,i]}")
##        if numFiresPopulationIsInside[i,j] > -1:
##        if numFiresPopulationIsInside[j,i] > -1:  # -1 indicates population is zero in that raster-cell, and therefore no population-node need be found within that cell for the population to snap to
#        if numFiresPopulationIsInside[j,i] > 0:  # -1 indicates population is zero in that raster-cell, and 0 indicates the cell is inside no fire-perimeters; in either case, no population-node need be found within that cell for the population to snap to
        if populationIsInsideFire[j,i] > 0:  # -1 indicates population is zero in that raster-cell, and 0 indicates the cell is not inside the fire-perimeter; in either case, no population-node need be found within that cell for the population to snap to
##          print(f"i {i}, j {j}, numFiresPopulationIsInside[{i},{j}] {numFiresPopulationIsInside[i,j]}")  # values are 0, 1, 2
#          print(f"j {j}, i {i}, numFiresPopulationIsInside[{j},{i}] {numFiresPopulationIsInside[j,i]}")  # values are 0, 1, 2
          print(f"j {j}, i {i}, populationIsInsideFire[{j},{i}] {populationIsInsideFire[j,i]}")  # values are 0, 1, 2
          b = BoundingBox.from_list([ [dims.ox+dims.hx*i, dims.oy+dims.hy*j], [dims.ox+dims.hx*(i+1), dims.oy+dims.hy*(j+1)] ])  # create search box for this cell
#          b = BoundingBox.from_list([ [dims.oy+dims.hy*j, dims.ox+dims.hx*i], [dims.oy+dims.hy*(j+1), dims.ox+dims.hx*(i+1)] ])  # create search box for this cell
          network_currentCell = network.region(b)  # find geometry within bounding-box
##          maxlinkcapacity_currentCell = maxlinkcapacity.currentCell(b)  # find raster within bounding-box: doesn't work, because (so far) maxlinkcapacity raster has incorrect values to begin with
#          maxlinkcapacity_currentCell = network_currentCell.rasterise(80, "output = max(0, capacity);", GeometryType.LineString)  # create finer-grained raster containing within each cell the maximum link-capacity: this works
          maxlinkcapacity_currentCell = network_currentCell.rasterise(CELLSIZEFINE, "output = max(0, capacity);", GeometryType.LineString)  # create finer-grained raster containing within each cell the maximum link-capacity: this works
#          maxlinkcapacity_currentCell = network_currentCell.rasterise(640, "output = max(0, capacity);", GeometryType.LineString)  # this works, and perhaps slightly faster than with cell-size of 80, but the difference is only about 1%
##          maxlinkcapacity_currentCell = network_currentCell.rasterise(1280, "output = max(0, capacity);", GeometryType.LineString)  # cell-size of 1280 FAILS
          largest_maxlinkcapacity_currentCell = maxlinkcapacity_currentCell.max()
          print(f"largest_maxlinkcapacity_currentCell is {largest_maxlinkcapacity_currentCell}")
##          largestMaxOutCapacityFound_cell = NEGATIVEVALUE  # will record, over all nodes/Points in the raster-cell, the node's maximum out-capacity that is the largest
#          INITIALLOWERBOUND_CURRENTCELL = maxlinkcapacity[j,i] - 1.0
#          INITIALLOWERBOUND_CURRENTCELL = 600.0 - 1.0  # problem-dependent, unless there is a way to establish that all links have capacity of at least 600.0?
          INITIALLOWERBOUND_CURRENTCELL = largest_maxlinkcapacity_currentCell - 1.0  # problem-dependent, unless there is a way to establish that all links have capacity of at least 600.0?
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
          else:
            MATsimID_node = None
          print(f"Property 'largestMaxOutCapacityFound_cell' is {largestMaxOutCapacityFound_cell}, at node with idx {idx_largestMaxOutCapacityFound_cell} and MATsimID_node {MATsimID_node}")
##          if largestMaxOutCapacityFound_cell == NEGATIVEVALUE:
#          if largestMaxOutCapacityFound_cell == INITIALLOWERBOUND_CURRENTCELL:
          if largestMaxOutCapacityFound_cell == None:
            print(f"WARNING: this raster-cell contains no node with a road-link leaving it.")
          network_currentCell.setProperty(idx_largestMaxOutCapacityFound_cell, "largestMaxOutCapacityInCell", largestMaxOutCapacityFound_cell)  # on the node inside the raster-cell that has the largest maximum out-capacity, set the property "largestMaxOutCapacityFound_cell" to this largest value
    print(f"inflowsByNodeID is {inflowsByNodeID}")
    with open(pth.join(outputDir, f"networkWithPropertiesByCell_{int(CELLSIZELARGE)}.geojson"), "w") as f:
      f.write(vectorToGeoJson(network))
#    endTime = time()
#    print(f"findEvacuationRisks.py ran for {endTime - startTime} seconds.")
#    sys.exit()

    # Find fire contour:
    print("Creating 'contourVector' layer ...")
    contourVector = isInsideFirePerimeter.vectorise( [ 0.5 ] )
##    contourVector.write(pth.join(outputDir, f"contourVector_{int(fire_max_time):06d}.tif"))  # returns "'AttributeError: 'Vector' object has no attribute 'write'"
#    with open(pth.join(outputDir, f"contourVector_{int(fire_max_time):06d}.geojson"), "w") as contourfile:
#      contourfile.write(vectorToGeoJson(contourVector))

#    runScript('numFiresInside += fire < fire::time ? 1 : 0;', [numFiresInside, fire] )  # value of raster-cell is 1 if inside the fire-perimeter, 0 otherwise

    # Find distance from fire contour:
    print("Creating 'dist' layer ...")
    dist = contourVector.mapDistance(resolution=180.0, geom_type=GeometryType.LineString, bounds=allFires_bounds)
    dist.name = 'dist'
    dist.setProjectionParameters(allFires_proj)
    runScript('dist = fire < fire::time ? -dist : dist;', [dist, fire] )  # ensure 'dist' is a signed-distance raster, where points inside the fire perimeter have a negative distance: this is needed so that exit-nodes will be selected only from the nodes that lie outside the perimeter and not inside it
#    dist.write(pth.join(outputDir, f"dist_{int(fire_max_time):06d}.tif"))  # write raster-layer 'dist' to a file

#    # Create count layer:
#    print("Creating 'count' layer ...")
#    count = Raster(name="count", data_type=np.uint32)
##    count.init_with_bbox(allFires_bounds, CELLSIZESMALL)
#    count.init_with_bbox(allFires_bounds, CELLSIZELARGE)
#    count.setProjectionParameters(allFires_proj)
#    # Count points in cell:
#    count.setAllCellValues(0)
#    count.rasterise(network, "atomic_inc();", GeometryType.Point)  # atomic_inc() increments the value in any cell containing vector-nodes from 'network' layer
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


    # Create exit-nodes based on safe distance from fire-perimeter:
    atSafeDistFromPerimeter = Raster(name="atSafeDistFromPerimeter")
    atSafeDistFromPerimeter.init_with_bbox(allFires_bounds, CELLSIZESMALL)
    atSafeDistFromPerimeter.setProjectionParameters(allFires_proj)
#    print(f"Creating 'atSafeDistFromPerimeter_{int(fire_max_time):06d}.tif' layer ...")
    print(f"Creating 'atSafeDistFromPerimeter_ffdi100{fireName}.tif' layer ...")
#    runScript("atSafeDistFromPerimeter = dist >= 1000.0 && dist < 2000.0 ? 1 : 0;", [atSafeDistFromPerimeter, dist])
    INNERDIST_SAFEBUFFER = 1000.0
    OUTERDIST_SAFEBUFFER = 2000.0
    runScript("atSafeDistFromPerimeter = dist >= INNERDIST_SAFEBUFFER && dist < OUTERDIST_SAFEBUFFER ? 1 : 0;", [atSafeDistFromPerimeter, dist])
###    atSafeDistFromPerimeter.write(pth.join(outputDir, f"exitnodesatsafedist_{int(fire_max_time):06d}.tif"))
##    atSafeDistFromPerimeter.write(pth.join(outputDir, f"atSafeDistFromPerimeter_{int(fire_max_time):06d}.tif"))
#    atSafeDistFromPerimeter.write(pth.join(outputDir, f"atSafeDistFromPerimeter__{int(rasterCellSize)}_ffdi100{fireName}.tif"))
    atSafeDistFromPerimeter.write(pth.join(outputDir, f"atSafeDistFromPerimeter__ffdi100{fireName}.tif"))

#    # Map flow to network:
    # Map population per node to network-node inflow:
#    network.pointSample(inflowPerNodeInsidePerimeter)  # sample 'inflowPerNodeInsidePerimeter' raster at each point and write resulting value to the 'network' vector-layer
    # Map whether sufficiently far from fire-perimeter to network-node safety:
    network.pointSample(atSafeDistFromPerimeter)  # sample 'atSafeDistFromPerimeter' raster at each point and write resulting value to the 'network' vector-layer
##    with open(pth.join(outputDir, f"networkNodesAtSafeDist_{int(fire_max_time):06d}.geojson"), "w") as networkfile:
#    with open(pth.join(outputDir, f"networkNodesAtSafeDist_{int(rasterCellSize)}_ffdi100{fireName}.geojson"), "w") as networkfile:
    with open(pth.join(outputDir, f"networkNodesAtSafeDist_ffdi100{fireName}.geojson"), "w") as networkfile:
      networkfile.write(vectorToGeoJson(network))

    if RUNMAXFLOW:
      # Run flow-solver:
#      networkFlowSolver.run()
      SEMversion = 'SEM5'  # static-flow model based on maximum-flow method#, with each specified subflow from an injection-node to an exit-node being added to the arcs on that path
      print(f"SEMversion is {SEMversion}.")
##      inflowsByNodeID = {'matsimnode0': 700, 'matsimnode1': 500}
#      inflowsByNodeID = getInjectionNodesandInflows()  # injection-nodes must be inside fire perimeter; TODO: for now, this is defined above by hand
      inflowsandflowperiodsByNodeID = None  # not needed for SEM5
##      exitnodes = {'matsimnode5', 'matsimnode8'}
#      exitnodes = getExitNodes()
      subflowsToAssignedExitNodesByInjectionNodeID = simulDurationInHours = None
      SEMoutputGeoJSON = solversSEMversions.runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, subflowsToAssignedExitNodesByInjectionNodeID, simulDurationInHours)


      # Write data:
#      outName, outExt = pth.splitext(outputGeoJSONfilename)
#      print(f"outName is {outName}, outExt is {outExt}")
##      outVector = networkFlowSolver.getNetwork().convert(network_proj)
#      sys.exit()  # for debugging, exit before writing any files

#      with open(f"{outName}_{int(fire_max_time):06d}_{outExt}", "w") as outfile:
      with open(outputGeoJSONfilename, "w") as outfile:
##        outfile.write(vectorToGeoJson(outVector))
#        outfile.write(SEMoutputGeoJSON)
        json.dump(SEMoutputGeoJSON, outfile)


#  numFiresInside.write(pth.join(outputDir, f"numFiresInside_{int(fire_max_time):06d}.tif"))  # write raster-layer 'numFiresInside' to a file
  endTime = time()
  print(f"findEvacuationRisks.py ran for {endTime - startTime} seconds.")
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