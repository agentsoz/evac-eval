# Run network-flow solver on simple scenario cmr_1s1d1r1k (Castlemaine region, one source, one destination, one road, one thousand agents).
# Adapted by Stephen Taylor from Vivian.py, written by Vivian Dabre (?) in 2020 (?)
import json

#from geostack.vector import Vector
from geostack.io import geoJsonToVector, vectorToGeoJson
from geostack.solvers import NetworkFlowSolver

#from geostack.core import ProjectionParameters

import sys

# Load network from JSON-file:  # TODO: allow user to enter the name of the input-file via a command-line option (use ArgumentParser or similar)
network = geoJsonToVector('cmr_1s1d1r_network_links-linkIDsadded.geojson')  # in this version of the input-file, each link's properties contain that link's ABM/MATSIM link-ID

#network = geoJsonToVector('test_data_1.zip')  # ST: has over 106,000 points, so probably too large
#network = geoJsonToVector('test_data_1.geojson')  # ST: unzipped from test_data_1.zip (alternative to loading in the .zip file directly)

# Does a flow of 1 into a node correspond to one thousand (1k) vehicles in the ABM?
##network.setProperty(0, "flow", 1)  # set inflow at node 0
#network.setProperty(3, "flow", 1)  # set inflow at node 3  # TODO: this node-ID should not be hand-coded: obtain it somehow from the input GeoJSON-file
network.setProperty(3, "flow", 1000)  # set inflow at node 3  # TODO: this node-ID should not be hand-coded: obtain it somehow from the input GeoJSON-file

#network.setProperty(0, "flow", 10)  # set inflow at node 0 - does this correspond to 10,000 vehicles in the ABM?

# Create solver:
networkConfig = {
    "constant": 100.0,
    "defaultLinkType": 1
}
networkFlowSolver = NetworkFlowSolver()
networkFlowSolver.init(network, json.dumps(networkConfig))

networkFlowSolver.run()  # run solver

network = networkFlowSolver.getNetwork()  # get flow-solver's internal network

#for idLS in network.getLineStringIndexes():  # ST: provide an ID for each link, that assigned by the internal network - but (COMPLETED) need each link's MATSIM ID as used within the ABM
#    network.setProperty(idLS, 'id', idLS)  # add the solver's internal network's ID for each line-string as a property of that feature

# Re-project network in EPSG:28355, as the original network for the ABM is in that projection and otherwise the SEM results will be in lon-lat (EPSG:4326):
# But re-projecting to EPSG:28355 probably won't solve Leorey's (round-off?) errors, because the input-file cmr_1s1d1r_network_links-linkIDsadded.geojson is in lon-lat and the Mount Alexander input applications/networkFlow/data/mtAlexanderShire/mount_alexander_shire_network_2018.links.geojson is also in lon-lat (although not its bounding-box - a Geostack bug?). These input-files' being in lon-lat results from ees2sem.py writing out a network's GeoJSON in EPSG:4326, which could be changed via Geostack Python-bindings to EPSG:28355, but James Hilton wrote on 16th March 2021 that "The geojson standard now specifies EPSG:4326 as standard and Geostack conforms to this."  Dhirendra Singh wrote on the same day that "if the expected projection for Geojson is 4326, we should avoid workflows where we force a different projection out (even if it is doable in Geostack). As James said, such un-standard files will eventually cause issues in standard platforms."
# Hence, the following two lines are commented out for now.
#proj_EPSG28355 = ProjectionParameters.from_proj4("+proj=utm +zone=55 +south +ellps=GRS80")
#network = network.convert(proj_EPSG28355)

#print("sys.exit()"); sys.exit()  # useful in debugging for exiting before a file is written

outputfilename = 'resultsSEM-cmr_1s1d1r1k.geojson'  # TODO: allow user to enter the name of the output-file via a command-line option
#outputfilename = 'resultsSEM-cmr_1s1d1r10k.geojson'
#print("networkFlowSolver has run; writing to resultsSEM-cmr_1s1d1r1k.geojson ...")
print(f"networkFlowSolver has run; writing to {outputfilename} ...")

# Write results:
#with open('resultsSEM-cmr_1s1d1r1k.geojson', 'w') as f:
with open(outputfilename, 'w') as f:
    f.write(vectorToGeoJson(network))
#    f.write(vectorToGeoJson(network, enforceProjection=False))  # don't force projection of output to EPSG:4326 - this is probably unwise, given James Hilton's advice of 16th March 2021 that "The geojson standard now specifies EPSG:4326 as standard and Geostack conforms to this. If you want you can reproject and force to write out to MGA zone 55, but this will likely not work with other GIS tools (the projection isn’t stored with the geojson so they can’t tell which projection it is)."
