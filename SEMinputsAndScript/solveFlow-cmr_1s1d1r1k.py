# Run static-evacuation model across several fires and towns.
## Run network-flow solver on simple scenario cmr_1s1d1r1k (Castlemaine region, one source, one destination, one road, one thousand agents).
# Adapted by Stephen Taylor from Vivian.py, written by Vivian Dabre (?) in 2020 (?)
import json

#from geostack.vector import Vector
from geostack.io import geoJsonToVector, vectorToGeoJson
from geostack.solvers import NetworkFlowSolver

#from geostack.core import ProjectionParameters

import sys
#from solversSEMversions import runSEM
import solversSEMversions
import time, gzip, csv

starttime = time.time()

# simulDurationInHours applies only to the SEM4 version of the model.
#simulDurationInHours = 0.016666667  # one minute's duration
#simulDurationInHours = 0.1  # six minutes' duration
#simulDurationInHours = 1.0
#simulDurationInHours = 2.0
#simulDurationInHours = 2.1  # in scenario cmr_3s2d, all wavefronts have exited the network by 2.0910 hours
#simulDurationInHours = 3.0
#simulDurationInHours = 5.0
#simulDurationInHours = 10.0
simulDurationInHours = 100.0  # to see the effects of running the simulation for an "infinite" length of time

#SEMversion = 'SEM1'  # original C++ version of SEM, with monotonic flow and an h-value at each node; the Newton-Raphson method finds solutions
#SEMversion = 'SEM2'  # adapted version of SEM1, with non-monotonic flow and an h-value at each node; a global optimiser finds solutions
#SEMversion = 'SEM3'  # version of SEM based on traffic-flow theory in chapters 7 and 8 of Treiber and Kesting (2013), with each vehicle using the shortest path to its nearest exit-node
#SEMversion = 'SEM4'  # dynamic-flow version of SEM3, with simulation running for a specified duration; not working yet for a set of specified subflows where each subflow travels from a given injection-node to a given exit-node
SEMversion = 'SEM5'  # static-flow model intended to run in shorter time than does SEM4, with each specified subflow from an injection-node to an exit-node being added to the arcs on that path
print(f"SEMversion is {SEMversion}.")

#SCENARIO = 'testdata'
#SCENARIO = 'onenodenetwork'  # a test-case, as the scenario with one node and no links ought to fail; this is because the number of injection-nodes is required to be both greater than zero and fewer than the total number of nodes, which makes it impossible to consistently define a set of injection-nodes in the one-node case
#SCENARIO = 'twonodenetwork'
#SCENARIO = 'threenodelinearnetwork'
##SCENARIO = 'threenodelinearnetworkFalseNegativeSEM2'
#SCENARIO = 'threenodelinearnetworkSEM2cannotSendSubflowsToAssignedExitNodes'
#SCENARIO = 'threenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes'  # bad example, as only solution is to split assigned subflow from 2 to 0 into two different paths
#SCENARIO = 'fivenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes'
#SCENARIO = 'threenodenetwork-LinkDirectionsAffectShortestPath'
#SCENARIO = 'threenodenetworkTwoExitNodes'
#SCENARIO = 'fournodelinearnetworkmissingendpoint'
#SCENARIO = 'fournodelinearnetwork'  # TODO: test this next
#SCENARIO = 'fournodenetwork-ThirdInjectionNodeHasNoPathToExit'
#SCENARIO = 'fournodenetworkTwoExitNodes'  # useful for developing SEM5
#SCENARIO = 'fivenodelinearnetwork'
#SCENARIO = 'fournodetree'
#SCENARIO = 'sixnodenontree'
#SCENARIO = 'sixnodeTwoWayLink'
#SCENARIO = 'sevennodenetworkwithcycle'
#SCENARIO = 'ninenodenontreeTwoInjectionNodesTwoExitNodes'
#SCENARIO = 'cmr_1s1d1r'
#SCENARIO = 'cmr_3s2d'  # if with total inflows of 3005, called 'cmr_3s2d3k'
#SCENARIO = 'cmr1full'  # scenario with full Castlemaine population?
SCENARIO = 'sixnodeTwoWayLink'

def extractAssignedSubflowsFromCSV( populationArchetypesCSVfile ):
#  subflowsToAssignedExitNodesByInjectionNodeID = {'258070680': {'977389105': 500, '236395531': 500}, '259608655': {'977389105': 500, '236395531': 500}, '60303151': {'977389105': 500, '236395531': 505}}
#  exitnodes = {'977389105', '236395531'}  # the exit-nodes as identified in README_3s2d3k.md; these nodes correspond respectively to D1 and D2
  with gzip.open(populationArchetypesCSVfile, mode='rt', newline='') as f:
    csvObject = csv.reader(f, delimiter = ',', quotechar='"')
    headers = next(csvObject, None)  # read header-row into a list of column-names
    print(f"headers is {headers}")
#    print(f"len(headers) is {len(headers)}")
    for i, row in enumerate(csvObject):
#      print("row is", row)
#      print("len(row) is", len(row))
      rowdata = {}
#      for j in range(len(row)):
#        rowdata[headers[j]] = row[j]
      for j, datum in enumerate(row):
#        print("j is", j)
        rowdata[headers[j]] = datum
      print("rowdata is", rowdata)
      geographicalCoordinate = rowdata['Geographical.Coordinate']
      evacLocationPreference = rowdata['EvacLocationPreference']
      print(f"geographicalCoordinate is '{geographicalCoordinate}', evacLocationPreference is '{evacLocationPreference}'")
      if i == 1:
        sys.exit()
  return (subflowsToAssignedExitNodesByInjectionNodeID, exitnodes)

# TODO: allow user to enter the name of the input-file via a command-line option (use ArgumentParser or similar):
if SCENARIO == 'testdata':
#  JSONnetworkfilename = 'test_data_1.zip'  # ST: has over 106,000 points, so probably too large
  JSONnetworkfilename = 'test_data_1.geojson'  # ST: unzipped from test_data_1.zip (alternative to loading in the .zip file directly)

elif SCENARIO == 'onenodenetwork':
  JSONnetworkfilename = 'onenodenetwork.geojson'  # two-node network with right-hand node missing, in the input-file 'twonodenetwork-missingEndpoint.geojson'
  inflowsByNodeID = {}
#  exitnodes = set()  # the empty set
  exitnodes = {0}

elif SCENARIO == 'twonodenetwork':
#  JSONnetworkfilename = 'twonodenetwork.geojson'
#  JSONnetworkfilename = 'twonodenetworkDuplicateEndpoint.geojson'
#  JSONnetworkfilename = 'twonodenetworkTwoDuplicateEndpoints.geojson'
#  JSONnetworkfilename = 'twonodenetworkMissingEndpoint.geojson'  # two-node network with right-hand node missing in the input-file
  JSONnetworkfilename = 'twonodenetworkNoMissingEndpoint.geojson'  # two-node network with both nodes given explicitly in the input-file
#  inflowNodeID = 0  # a single inflow at node 0
  inflowsByNodeID = {'matsimnode0': 400}  # this flow exceeds the single link's capacity
#  inflowsByNodeID = {'matsimnode0': 401}  # this flow exceeds the single link's capacity
#  inflowsByNodeID = {'matsimnode0': 501}  # this flow exceeds the single link's capacity
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode1': 400}}
  exitnodes = {'matsimnode1'}

elif SCENARIO == 'threenodelinearnetwork':
  JSONnetworkfilename = 'threenodelinearnetworkMissingEndpoint.geojson'  # three-node network with rightmost node missing
#  inflowsByNodeID = {0: 500}
#  inflowsByNodeID = {0: 1000}
  inflowsByNodeID = {'matsimnode0': 1000}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode2': 1000}}
  exitnodes = {'matsimnode2'}

#elif SCENARIO == 'threenodelinearnetworkFalseNegativeSEM2':
elif SCENARIO == 'threenodelinearnetworkSEM2cannotSendSubflowsToAssignedExitNodes':
#  JSONnetworkfilename = 'threenodelinearnetworkFalseNegativeSEM2.geojson'  # three-node network (leftmost and rightmost nodes being exit-nodes) that formerly, in the case of flow being free to go to any exit-node, demonstrated SEM2 finding a false negative (a solution with no congestion where in reality traffic might all take the short left link) - it now demonstrates something similar to a false negative, namely that in the case of the injection-node having one subflow, assigned to exit-node 0, SEM2 cannot send all flow to the assigned exit-node (even when the lefthand edge has link-capacity sufficiently high to accommodate that flow); here, all traffic has been assigned to take the short left link but SEM2 allows only a solution in which some traffic takes the right link towards an unassigned exit-node
  JSONnetworkfilename = 'threenodelinearnetworkSEM2cannotSendSubflowsToAssignedExitNodes.geojson'  # three-node network (leftmost and rightmost nodes being exit-nodes) that demonstrates that SEM2 cannot send every subflow to its assigned exit-node (even when link-capacities are sufficiently high to accommodate all flows): in the case of the injection-node having one subflow, assigned to exit-node 0, SEM2 cannot send all flow to the assigned exit-node even when the lefthand edge has link-capacity sufficiently high to accommodate that flow; SEM2 could solve this network correctly by allowing arbitrary h-values at exit-nodes, but scenario 'threenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes.geojson' below demonstrates that that improvement wouldn't allow all four-node networks to be solved correctly
  inflowsByNodeID = {'matsimnode1': 1000}
#  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode1': {'matsimnode0': 600, 'matsimnode2': 400}}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode1': {'matsimnode0': 1000}}
  exitnodes = {'matsimnode0', 'matsimnode2'}  # left- and righthand nodes, respectively

elif SCENARIO == 'threenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes':
  JSONnetworkfilename = 'threenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes.geojson'  # three-node network that demonstrates SEM2 cannot send every subflow to its assigned exit-node (even when link-capacities are sufficiently high to accommodate all flows), even if we allow arbitrary h-values at exit-nodes. (We don't implement arbitrary h-values at exit-nodes, but use this example to deduce that that wouldn't give SEM2 the ability to correctly solve this network and its set of assigned subflows.)
  inflowsByNodeID = {'matsimnode1': 1000, 'matsimnode2': 1000}
#  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode1': {'matsimnode0': 600, 'matsimnode2': 400}}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode1': {'matsimnode0': 1000}, 'matsimnode2': {'matsimnode0': 1000}}
  exitnodes = {'matsimnode0'}  # leftmost node

elif SCENARIO == 'fivenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes':
  JSONnetworkfilename = 'fivenodenetworkSEM2cannotSendSubflowsToAssignedExitNodes.geojson'  # five-node network that demonstrates SEM2 cannot send every subflow to its assigned exit-node (even when link-capacities are sufficiently high to accommodate all flows), even if we allow arbitrary h-values at exit-nodes. (We don't implement arbitrary h-values at exit-nodes, but use this example to deduce that that wouldn't give SEM2 the ability to correctly solve this network and its set of assigned subflows.)
  inflowsByNodeID = {'matsimnode1': 100, 'matsimnode2': 1000}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode1': {'matsimnode0': 100}, 'matsimnode2': {'matsimnode0': 500, 'matsimnode4': 500}}
  exitnodes = {'matsimnode0', 'matsimnode4'}

elif SCENARIO == 'threenodenetwork-LinkDirectionsAffectShortestPath':
  JSONnetworkfilename = 'threenodenetwork-LinkDirectionsAffectShortestPath.geojson'  # three-node network with the topmost link shorter in terms of traversal-time but going in the wrong direction (from the exit-node at right back to the injection-node at top left)
  inflowsByNodeID = {'matsimnode0': 1000}  # top-left node
#  inflowsByNodeID = {'matsimnode0': 1000, 'matsimnode1': 0}  # top-left and middle-bottom nodes
  exitnodes = {'matsimnode2'}  # topright node

elif SCENARIO == 'threenodenetworkTwoExitNodes':
  JSONnetworkfilename = 'threenodenetworkTwoExitNodes.geojson'  # three-node network with one injection-flow at the lefthand node directed to either of two exit-nodes on the right
#  inflowsByNodeID = {'matsimnode0': 1000}  # lefthand node
#  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode1': 700, 'matsimnode2': 300}}
  inflowsByNodeID = {'matsimnode0': 1600}  # lefthand node
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode1': 900, 'matsimnode2': 700}}
  exitnodes = {'matsimnode1', 'matsimnode2'}  # righthand nodes, top and bottom respectively

elif SCENARIO == 'fournodenetworkTwoExitNodes':
  JSONnetworkfilename = 'fournodenetworkTwoExitNodes.geojson'  # four-node network with one injection-flow at the lefthand node directed to either of two exit-nodes on the right
#  inflowsByNodeID = {'matsimnode0': 1000}  # lefthand node
#  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode2': 700, 'matsimnode3': 300}}
  inflowsByNodeID = {'matsimnode0': 1000}  # lefthand node
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode2': 800, 'matsimnode3': 200}}
#  inflowsByNodeID = {'matsimnode0': 1400}  # lefthand node
#  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode2': 900, 'matsimnode3': 500}}
  exitnodes = {'matsimnode2', 'matsimnode3'}  # righthand nodes, top and bottom respectively

elif SCENARIO == 'fournodelinearnetworkmissingendpoint':
  JSONnetworkfilename = 'fournodelinearnetworkMissingEndpoint.geojson'  # four-node linear network with rightmost node missing
  inflowsByNodeID = {0: 700}
#  inflowsByNodeID = {0: 700, 1: 200}
  exitnodes = {3}

elif SCENARIO == 'fournodelinearnetwork':
  JSONnetworkfilename = 'fournodelinearnetwork.geojson'  # four-node linear network with rounded-off values such as 50,000 and 20,000 metres for link-lengths
#  inflowsByNodeID = {'matsimnode0': 700, 'matsimnode1': 500}
  inflowsByNodeID = {'matsimnode0': 400, 'matsimnode1': 200}
  exitnodes = {'matsimnode3'}

elif SCENARIO == 'fournodenetwork-ThirdInjectionNodeHasNoPathToExit':
  JSONnetworkfilename = 'fournodenetwork-ThirdInjectionNodeHasNoPathToExit.geojson'  # four-node network in which rightmost injection-node can reach no exit-node because only a one-way link reaches it
#  inflowsByNodeID = {'matsimnode0': 1000, 'matsimnode3': 500}
  inflowsByNodeID = {'matsimnode0': 1000, 'matsimnode1': 0, 'matsimnode3': 500}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode2': 1000}, 'matsimnode3': {'matsimnode2': 500}}
  exitnodes = {'matsimnode2'}  # middle node on upper linear part of network

elif SCENARIO == 'fivenodelinearnetwork':
  JSONnetworkfilename = 'fivenodelinearnetwork.geojson'  # five-node linear network with right-hand node missing
  inflowsByNodeID = {'matsimnode0': 700}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode4': 700}}
  exitnodes = {'matsimnode4'}

elif SCENARIO == 'fournodetree':
#  JSONnetworkfilename = 'fournodetreeFourthNodedeleted.geojson'
  JSONnetworkfilename = 'fournodetree.geojson'
##  inflowsByNodeID = {1: 200, 3: 100}  # desired labelling of nodes
#  inflowsByNodeID = {0: 200, 1: 100}  # to conform with default node-labelling that begins at 0; TODO: improve this so that labels can be any integers not necessarily including zero?
  inflowsByNodeID = {'matsimnode0': 300, 'matsimnode1': 400, 'matsimnode2': 100}
#  inflowsByNodeID = {0: 200, 1: 401}
  exitnodes = {'matsimnode3'}

elif SCENARIO == 'sixnodenontree':
  JSONnetworkfilename = 'sixnodenontree.geojson'
  inflowsByNodeID = {'matsimnode0': 601, 'matsimnode2': 1}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode5': 601}, 'matsimnode2': {'matsimnode5': 1}}
#  inflowsByNodeID = {'matsimnode0': 500, 'matsimnode2': 30}
#  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode5': 500}, 'matsimnode2': {'matsimnode5': 30}}
  exitnodes = {'matsimnode5'}

elif SCENARIO == 'sixnodeTwoWayLink':
  JSONnetworkfilename = 'sixnodeTwoWayLink.geojson'
  inflowsByNodeID = {'matsimnode0': 200, 'matsimnode1': 450}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode5': 200}, 'matsimnode1': {'matsimnode4': 450}}
  exitnodes = {'matsimnode4', 'matsimnode5'}

elif SCENARIO == 'sevennodenetworkwithcycle':
  JSONnetworkfilename = 'sevennodenetworkwithcycle.geojson'
  inflowsByNodeID = {'matsimnode0': 100, 'matsimnode1': 300, 'matsimnode5': 600}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode6': 100}, 'matsimnode1': {'matsimnode3': 300}, 'matsimnode5': {'matsimnode4': 200, 'matsimnode6': 400}}
  exitnodes = {'matsimnode3', 'matsimnode4', 'matsimnode6'}

elif SCENARIO == 'ninenodenontreeTwoInjectionNodesTwoExitNodes':
  JSONnetworkfilename = 'ninenodenontreeTwoInjectionNodesTwoExitNodes.geojson'
  inflowsByNodeID = {'matsimnode0': 800, 'matsimnode6': 400}
  subflowsToAssignedExitNodesByInjectionNodeID = {'matsimnode0': {'matsimnode5': 600, 'matsimnode8': 200}, 'matsimnode6': {'matsimnode5': 100, 'matsimnode8': 300}}
  exitnodes = {'matsimnode5', 'matsimnode8'}

elif SCENARIO == 'cmr_1s1d1r':
##  JSONnetworkfilename = 'cmr_1s1d1r_network_linksLinkIDsadded.geojson'  # symbolic link to ~/staticevacuation/condaenvGeostack/evac-eval/modelinputsCmr_1s1d1r1k/sem/cmr_1s1d1r_network.links.geojson' created by ees2sem.py. In this version of the input-file, each link's properties contain that link's ABM/MATSIM link-ID
#  JSONnetworkfilename = 'cmr_1s1d1r_network.linksFinalNodeDeleted.geojson'
  JSONnetworkfilename = 'cmr_1s1d1r_network.linksFinalNodeDeleted-linksInTopologicalOrder.geojson'
#  inflowNodeID = 3  # a single inflow at node 3
#  inflowsByNodeID = {3: 1000}
#  inflowsByNodeID = {5:300, 6:700}
#  inflowsByNodeID = {3: 1000, 5:300, 6:700}
  inflowsByNodeID = {3: 1000, 4:1500, 5:300, 6:700}
#  inflowsByNodeID = {3: 10000}
  exitnodes = set()  # the empty set  # FIXME: which is the exit-node? (Is the one that's missing as an end-point from one or more LineStrings)

elif SCENARIO == 'cmr_3s2d':
#  JSONnetworkfilename = 'cmr_3s2d-linksInTopologicalOrder.geojson'
  JSONnetworkfilename = 'cmr_3s2d3k_network.links.geojson'
  inflowsByNodeID = {'258070680': 1000, '259608655':1000, '60303151':1005}  # the exit-nodes as identified in README_3s2d3k.md; these nodes correspond respectively to S1, S2, and S3
  subflowsToAssignedExitNodesByInjectionNodeID = {'258070680': {'977389105': 500, '236395531': 500}, '259608655': {'977389105': 500, '236395531': 500}, '60303151': {'977389105': 500, '236395531': 505}}
  exitnodes = {'977389105', '236395531'}  # the exit-nodes as identified in README_3s2d3k.md; these nodes correspond respectively to D1 and D2

elif SCENARIO == 'cmr1full':
  JSONnetworkfilename = 'cmr1full_network.links.geojson'
#  populationArchetypesCSVfile = '/home/tay373/staticevacuation/condaenvGeostack/evac-eval/modelinputsCmr1full/sem/population-archetypes.csv.gz'
#  starttime_extractAssignedSubflowsFromCSV  = time.time()
#  (subflowsToAssignedExitNodesByInjectionNodeID, exitnodes) = extractAssignedSubflowsFromCSV( populationArchetypesCSVfile )
##  endtime_extractAssignedSubflowsFromCSV = time.time()
##  timeToExtractAssignedSubflowsFromCSV = endtime_extractAssignedSubflowsFromCSV - starttime_extractAssignedSubflowsFromCSV
#  timeToExtractAssignedSubflowsFromCSV = time.time() - starttime_extractAssignedSubflowsFromCSV
#  print(f"Extraction of assigned subflows and exit-nodes from populationArchetypesCSVfile {populationArchetypesCSVfile} took {timeToExtractAssignedSubflowsFromCSV:.5f} seconds.")
  assignedSubflowsGeoJSONfile = 'assignedSubflows.geojson'
  (subflowsToAssignedExitNodesByInjectionNodeID, exitnodes) = readAssignedSubflowsAndExitNodesFromGeoJSON( assignedSubflowsGeoJSONfile )

print(f"JSONnetworkfilename is {JSONnetworkfilename}")

#if SEMversion == 'SEM4':
inflowsandflowperiodsByNodeID = {}
for i in inflowsByNodeID:
#  inflowsandflowperiodsByNodeID[i] = {'inflow': inflowsByNodeID[i], 'inflowstarttime': 0, 'duration': 1.0}
  inflowsandflowperiodsByNodeID[i] = {'inflow': inflowsByNodeID[i], 'inflowstarttime': 0, 'duration': 2.0}


#if inflowAtSingleNode % 1000 == 0:
if len(inflowsByNodeID) == 1:
  if list( inflowsByNodeID.values() )[0] % 1000 == 0:
#    inflowstring = str( int(inflowAtSingleNode / 1000) ) + 'k'  # string representing the inflow, to be used in the name of the output-file
    inflowstring = str( int(list( inflowsByNodeID.values() )[0] / 1000) ) + 'k'  # string representing the inflow, to be used in the name of the output-file
  else:
#    inflowstring = str(inflowAtSingleNode)
    inflowstring = str( list( inflowsByNodeID.values() )[0] )
else:
  inflowstring = ('-').join( [str(v) for v in inflowsByNodeID.values()] )
print(f"inflowstring is {inflowstring}")


# Load network from JSON-file:
if SEMversion == 'SEM1':
  network = geoJsonToVector(JSONnetworkfilename)  # in this version of the input-file, each link's properties include its ABM/MATSIM link-ID, capacity, diameter (probably not needed), free speed, length, number of lanes, and 'oneway' (Boolean?)

# A flow of 1.0 into a node seemingly corresponds to one vehicle per hour in the ABM:
###  network.setProperty(3, "flow", 10000)  # set inflow at node 3 - corresponds to 10,000 vehicles per hour in the ABM
##  network.setProperty(3, "flow", 1000)  # set inflow at node 3  # TODO: this node-ID should not be hand-coded: obtain it somehow from the input GeoJSON-file (it's the northernmost node, so has the least negative latitude)
#  network.setProperty(inflowNodeID, "flow", inflowAtSingleNode)  # set inflow at a single node
  for nodeID, inflowValue in inflowsByNodeID:
    network.setProperty(nodeID, "flow", inflowValue)  # set inflow at a single node

elif SEMversion in ('SEM2', 'SEM3', 'SEM4', 'SEM5'):
#  inflowsByNodeID = {inflowNodeID: inflowAtSingleNode}  # in general, will contain all non-zero (positive) inflows and the corresponding node-IDs
  pass


# Create solver:
if SEMversion == 'SEM1':
  # N.B. "constant" is used only when flow is monotonic in density, and only by the Hazen-Williams and Manning Open Channel models; it takes a floating-point value that is the link's pipe-roughness coefficient C (see 'Swift: a GPU-based coupled hydrodynamic/hydraulic framework for urban flood-prediction', 2015).
  # N.B. "defaultLinkType" takes an integer value that Geostack maps via the C++ enum NetworkSegmentType::Type as
  #   0: Undefined,
  #   1: HazenWilliams,
  #   2: ManningOpenChannel,
  #   3: Logarithmic,
  #   4: SqrtExp.
  networkConfig = {
    "constant": 100.0,
    "defaultLinkType": 1
  }
  networkFlowSolver = NetworkFlowSolver()
  networkFlowSolver.init(network, json.dumps(networkConfig))
  
  networkFlowSolver.run()  # run solver
  
  network = networkFlowSolver.getNetwork()  # get flow-solver's internal network


elif SEMversion in ('SEM2', 'SEM3', 'SEM4', 'SEM5'):
  outputGeoJSON = solversSEMversions.runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, subflowsToAssignedExitNodesByInjectionNodeID, simulDurationInHours)


#for idLS in network.getLineStringIndexes():  # ST: provide an ID for each link, that assigned by the internal network - but (COMPLETED) need each link's MATSIM ID as used within the ABM
#  network.setProperty(idLS, 'id', idLS)  # add the solver's internal network's ID for each line-string as a property of that feature

# Re-project network in EPSG:28355, as the original network for the ABM is in that projection and otherwise the SEM results will be in lon-lat (EPSG:4326):
# But re-projecting to EPSG:28355 probably won't solve Leorey's (round-off?) errors, because the input-file cmr_1s1d1r_network_links-linkIDsadded.geojson is in lon-lat and the Mount Alexander input applications/networkFlow/data/mtAlexanderShire/mount_alexander_shire_network_2018.links.geojson is also in lon-lat (although not its bounding-box - a Geostack bug?). These input-files' being in lon-lat results from ees2sem.py writing out a network's GeoJSON in EPSG:4326, which could be changed via Geostack Python-bindings to EPSG:28355, but James Hilton wrote on 16th March 2021 that "The geojson standard now specifies EPSG:4326 as standard and Geostack conforms to this."  Dhirendra Singh wrote on the same day that "if the expected projection for Geojson is 4326, we should avoid workflows where we force a different projection out (even if it is doable in Geostack). As James said, such un-standard files will eventually cause issues in standard platforms."
# Hence, the following two lines are commented out for now.
#proj_EPSG28355 = ProjectionParameters.from_proj4("+proj=utm +zone=55 +south +ellps=GRS80")
#network = network.convert(proj_EPSG28355)

# TODO: allow user to enter the name of the output-file via a command-line option.
##outputfilename = 'resultsSEM-cmr_1s1d1r1k.geojson'
#outputfilename = f'results{SEMversion}-cmr_1s1d1r{inflowstring}.geojson'
assert JSONnetworkfilename[-8:] == '.geojson'
inputname = JSONnetworkfilename[:-8]
outputfilename = f'results{SEMversion}-{inputname}-{inflowstring}-oneline.geojson'  # all on one line, so not human-readable
print(f"networkFlowSolver has run; writing to {outputfilename} ...")

print("Running sys.exit() so that no output-file is written."); sys.exit()  # useful in debugging, for exiting before a file is written

# Write results:
#with open('resultsSEM-cmr_1s1d1r1k.geojson', 'w') as f:
with open(outputfilename, 'w') as f:
  if SEMversion == 'SEM1':
    f.write(vectorToGeoJson(network))
#  f.write(vectorToGeoJson(network, enforceProjection=False))  # don't force projection of output to EPSG:4326 - this is probably unwise, given James Hilton's advice of 16th March 2021 that "The geojson standard now specifies EPSG:4326 as standard and Geostack conforms to this. If you want you can reproject and force to write out to MGA zone 55, but this will likely not work with other GIS tools (the projection isn’t stored with the geojson so they can’t tell which projection it is)."
  elif SEMversion in ('SEM2', 'SEM3', 'SEM4', 'SEM5'):
#    f.write(outputGeoJSON)
    json.dump(outputGeoJSON, f)

endtime = time.time()
runtime = endtime - starttime
print(f"Entire run took {runtime:.5f} seconds.")
