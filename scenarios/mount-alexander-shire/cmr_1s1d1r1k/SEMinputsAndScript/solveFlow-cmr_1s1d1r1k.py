# Run network-flow solver on simple scenario cmr_1s1d1r1k (Castlemaine region, one source, one destination, one road, one thousand agents).
# Adapted by Stephen Taylor from Vivian.py, written by Vivian Dabre (?) in 2020 (?)
import json

#from geostack.vector import Vector
from geostack.io import geoJsonToVector, vectorToGeoJson
from geostack.solvers import NetworkFlowSolver

#from geostack.core import ProjectionParameters

import sys
#from solversSEM2and3and4 import runSEM2or3or4
import solversSEM2and3and4

#SEMversion = 'SEM1'  # original C++ version of SEM, with monotonic flow and an h-value at each node; the Newton-Raphson method finds solutions
#SEMversion = 'SEM2'  # second version of SEM, with non-monotonic flow and an h-value at each node; a global optimiser finds solutions
#SEMversion = 'SEM3'  # third version of SEM, with each vehicle using its shortest path to an exit-node
SEMversion = 'SEM4'  # dynamic-flow version of SEM3, with simulation running for a specified duration

#SCENARIO = 'testdata'
SCENARIO = 'cmr_1s1d1r'
#SCENARIO = 'onenodenetwork'
#SCENARIO = 'twonodenetwork'
#SCENARIO = 'threenodelinearnetwork'
#SCENARIO = 'fournodelinearnetwork'
#SCENARIO = 'fivenodelinearnetwork'
#SCENARIO = 'fournodetree'
#SCENARIO = 'sixnodenontree'

# TODO: allow user to enter the name of the input-file via a command-line option (use ArgumentParser or similar):
if SCENARIO == 'testdata':
#  inputfilename = 'test_data_1.zip'  # ST: has over 106,000 points, so probably too large
  inputfilename = 'test_data_1.geojson'  # ST: unzipped from test_data_1.zip (alternative to loading in the .zip file directly)

elif SCENARIO == 'cmr_1s1d1r':
##  inputfilename = 'cmr_1s1d1r_network_linksLinkIDsadded.geojson'  # symbolic link to ~/staticevacuation/condaenvGeostack/evac-eval/modelinputsCmr_1s1d1r1k/sem/cmr_1s1d1r_network.links.geojson' created by ees2sem.py. In this version of the input-file, each link's properties contain that link's ABM/MATSIM link-ID
#  inputfilename = 'cmr_1s1d1r_network.linksFinalNodeDeleted.geojson'
  inputfilename = 'cmr_1s1d1r_network.linksFinalNodeDeleted-linksInTopologicalOrder.geojson'
#  inflowNodeID = 3  # a single inflow at node 3
#  inflowAtSingleNode = 1000
  inflowsByNodeID = {3: 1000}
#  inflowsByNodeID = {3: 10000}
  exitnodes = set()  # the empty set  # FIXME: which is the exit-node? (Is the one that's missing as an end-point from one or more LineStrings)

elif SCENARIO == 'onenodenetwork':  # a test-case scenario, as the scenario with one node and no links ought to fail; this is because the number of injection-nodes is required to be both greater than zero and fewer than the total number of nodes, which makes it impossible to consistently define a set of injection-nodes in the one-node case
  inputfilename = 'onenodenetwork.geojson'  # two-node network with right-hand node missing, in the input-file 'twonodenetwork-missingEndpoint.geojson'
  inflowsByNodeID = {}
#  exitnodes = set()  # the empty set
  exitnodes = {0}

elif SCENARIO == 'twonodenetwork':
#  inputfilename = 'twonodenetwork.geojson'
#  inputfilename = 'twonodenetworkDuplicateEndpoint.geojson'
#  inputfilename = 'twonodenetworkTwoDuplicateEndpoints.geojson'
  inputfilename = 'twonodenetworkMissingEndpoint.geojson'  # two-node network with right-hand node missing, in the input-file 'twonodenetwork-missingEndpoint.geojson'
#  inflowNodeID = 0  # a single inflow at node 0
##  inflowAtSingleNode = 400  # to use with input-file 'twonodenetwork-missingEndpoint.geojson'
##  inflowAtSingleNode = 500
#  inflowAtSingleNode = 501
##  inflowAtSingleNode = 600  # this flow exceeds the single link's capacity
#  inflowsByNodeID = {0: 400}
  inflowsByNodeID = {0: 501}
  exitnodes = {1}

elif SCENARIO == 'threenodelinearnetwork':
  inputfilename = 'threenodelinearnetworkMissingEndpoint.geojson'  # three-node network with right-hand node missing
#  inflowsByNodeID = {0: 500}
  inflowsByNodeID = {0: 1000}
  exitnodes = {2}

elif SCENARIO == 'fournodelinearnetwork':
  inputfilename = 'fournodelinearnetworkMissingEndpoint.geojson'  # four-node linear network with right-hand node missing
  inflowsByNodeID = {0: 700}
#  inflowsByNodeID = {0: 700, 1: 200}
  exitnodes = {3}

elif SCENARIO == 'fivenodelinearnetwork':
  inputfilename = 'fivenodelinearnetworkMissingEndpoint.geojson'  # five-node linear network with right-hand node missing
  inflowsByNodeID = {0: 700}
  exitnodes = {4}

elif SCENARIO == 'fournodetree':
  inputfilename = 'fournodetreeFourthNodedeleted.geojson'
##  inflowsByNodeID = {1: 200, 3: 100}  # desired labelling of nodes
#  inflowsByNodeID = {0: 200, 1: 100}  # to conform with default node-labelling that begins at 0; TODO: improve this so that labels can be any integers not necessarily including zero?
  inflowsByNodeID = {0: 300, 1: 400, 2: 100}
#  inflowsByNodeID = {0: 200, 1: 401}
  exitnodes = {3}

elif SCENARIO == 'sixnodenontree':
  inputfilename = 'sixnodenontreeNode5deleted.geojson'
##  inflowNodeID = 0  # a single inflow at node 0
##  inflowAtSingleNode = 500
#  inflowsByNodeID = {0: 500}
  inflowsByNodeID = {0: 601, 2: 1}
  exitnodes = {5}

print(f"inputfilename is {inputfilename}")

#if SEMversion == 'SEM4':
inflowsandflowperiodsByNodeID = {}
for i in inflowsByNodeID:
  inflowsandflowperiodsByNodeID[i] = {'inflow': inflowsByNodeID[i], 'starttime': 0, 'duration': 1.0}


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
  network = geoJsonToVector(inputfilename)  # in this version of the input-file, each link's properties include its ABM/MATSIM link-ID, capacity, diameter (probably not needed), free speed, length, number of lanes, and 'oneway' (Boolean?)

# A flow of 1.0 into a node seemingly corresponds to one vehicle per hour in the ABM:
###  network.setProperty(3, "flow", 10000)  # set inflow at node 3 - corresponds to 10,000 vehicles per hour in the ABM
##  network.setProperty(3, "flow", 1000)  # set inflow at node 3  # TODO: this node-ID should not be hand-coded: obtain it somehow from the input GeoJSON-file (it's the northernmost node, so has the least negative latitude)
#  network.setProperty(inflowNodeID, "flow", inflowAtSingleNode)  # set inflow at a single node
  for nodeID, inflowValue in inflowsByNodeID:
    network.setProperty(nodeID, "flow", inflowValue)  # set inflow at a single node

elif SEMversion in ('SEM2', 'SEM3', 'SEM4'):
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


elif SEMversion in ('SEM2', 'SEM3', 'SEM4'):
#  outputGeoJSON = solversSEM2and3and4.runSEM2or3or4( SEMversion, inputfilename, inflowsByNodeID)  # set inflow at node 0, for two-node network with right-hand node missing from input-file 'twonodenetwork-missingEndpoint.geojson'
  outputGeoJSON = solversSEM2and3and4.runSEM2or3or4( SEMversion, inputfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes)  # set inflow at node 0, for two-node network with right-hand node missing from input-file 'twonodenetwork-missingEndpoint.geojson'


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
assert inputfilename[-8:] == '.geojson'
inputname = inputfilename[:-8]
outputfilename = f'results{SEMversion}-{inputname}-{inflowstring}.geojson'
print(f"networkFlowSolver has run; writing to {outputfilename} ...")

#print("sys.exit()"); sys.exit()  # useful in debugging for exiting before a file is written

# Write results:
#with open('resultsSEM-cmr_1s1d1r1k.geojson', 'w') as f:
with open(outputfilename, 'w') as f:
  if SEMversion == 'SEM1':
    f.write(vectorToGeoJson(network))
#  f.write(vectorToGeoJson(network, enforceProjection=False))  # don't force projection of output to EPSG:4326 - this is probably unwise, given James Hilton's advice of 16th March 2021 that "The geojson standard now specifies EPSG:4326 as standard and Geostack conforms to this. If you want you can reproject and force to write out to MGA zone 55, but this will likely not work with other GIS tools (the projection isn’t stored with the geojson so they can’t tell which projection it is)."
  elif SEMversion in ('SEM2', 'SEM3', 'SEM4'):
#    f.write(outputGeoJSON)
    json.dump(outputGeoJSON, f)
