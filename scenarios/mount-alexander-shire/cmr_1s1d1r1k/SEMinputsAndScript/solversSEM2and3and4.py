#!/usr/bin/python3

# TODO: parameter-validation for all functions.

import scipy.optimize
#print("scipy.__version__ is", scipy.__version__)  # need a later version such as 1.4.1 in order to use multivariate scipy.optimize.newton(), i.e. call it on a vector-valued function and a vector-valued initial estimate
import sys
#sys.exit()
import numpy as np
import math
import scipy.linalg
import random
import networkx as nx
import matplotlib.pyplot as plt
import time
import scipy.stats
import json

#meanVehicleLength_Km = .005  # these two shouldn't be global variables, so their definitions have been moved inside function 'setuptrafficflowproblem()'
#meanMinGapinStoppedTraffic_Km = .004
#EPS = 1e-3  # this shouldn't be a global variable, so its definition has been moved inside each of the functions 'adjustSpeedandDensityToLink()' and 'flowEnteringNode()'


def phiOnLinkrs_inVehiclesPerKmPerLane( lK_rs, meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km ):  # Calculate value of \phi_{rs} for link {r,s} given K_{rs} and two constants; 'meanMinGapinStoppedTraffic_Km' is the gap between the front of a car and the rear of the car immediately ahead of it
  speedCameraPrecision_kmperhr = 2
  stepsize = 1.0
#  if speedCameraPrecision_kmperhr >= lK_rs:  # these two lines of code might be needed if lK_rs is ever <= speedCameraPrecision_kmperhr, but choosing lK_rs to be small is so far useful only for forcing a link to have low capacity, which from the numerical point of view seems to be better achieved by setting the link's number of lanes to a low value such as .01. (Forcing a link to have low capacity is useful for creating an effective dummy link adjacent to an injection-node, each of which must have edge-degree >= 2.)
#    speedCameraPrecision_kmperhr = lK_rs - stepsize  # might be incorrect
  prevApprox = 0.0
  currentApprox = 1e99
  while abs(currentApprox - prevApprox) > 5e-4:
    prevApprox = currentApprox
    lowspeed = np.arange(stepsize, speedCameraPrecision_kmperhr + stepsize, stepsize)  # values are between, essentially, 0 and speedCameraPrecision_kmperhr
    phi = (meanMinGapinStoppedTraffic_Km + meanVehicleLength_Km)**(-1) / np.log( lK_rs / lowspeed )
    currentApprox = sum(phi) / len(phi)
    print("currentApprox is %.5f" % currentApprox) 
    stepsize *= .1
  return currentApprox


def geoJSONtoAdjacencyMatrix( JSONnetworkfilename, exitnodes ):
  with open(JSONnetworkfilename, 'r') as f:
    jsonobj = json.load(f)
#  print("jsonobj is", jsonobj)
  assert jsonobj["type"] == 'FeatureCollection'
  numDuplicatePoints = 0  # count duplicate Point-features
  numDuplicateLinestrings = 0  # count duplicate LineString-features
  pointsWithCoords = {}
  linestringsWithCoords = {}
  currentPointID = 0
  currentLinestringID = 0
  for feature in jsonobj["features"]:
    assert feature["type"] == 'Feature'
    geometry = feature["geometry"]
    if geometry['type'] == "Point":
      key = tuple( geometry['coordinates'] )
      if key not in pointsWithCoords:  # add each non-duplicate Point
        pointsWithCoords[ key ] = currentPointID
        currentPointID += 1
      else:
        print(f"A Point with coordinates {key} already exists in pointsWithCoords")
        numDuplicatePoints += 1
    elif geometry['type'] == "LineString":
      assert len( geometry['coordinates'] ) == 2
      key = []  # construct a dictionary key, a two-tuple of two-tuples, from the LineString's pair of Point-coordinates
      for pointcoords in geometry['coordinates']:
        key.append( tuple(pointcoords) )  # append this two-tuple
      key = tuple(key)  # two-tuple of two-tuples
#      print(f"key is {key}")
      if key not in linestringsWithCoords:
        linestringsWithCoords[ key ] = {'id': currentLinestringID}
        currentLinestringID += 1
#        for propertyname in feature['properties']:
#          linestringsWithCoords[key][propertyname] = feature['properties'][propertyname]
        linestringsWithCoords[key]['linkcapacity'] = float( feature['properties']['capacity'] )
        linestringsWithCoords[key]['diameter'] = int( feature['properties']['diameter'] )
        linestringsWithCoords[key]['freespeed'] = float( feature['properties']['freespeed'] )
        linestringsWithCoords[key]['length'] = float( feature['properties']['length'] )
        linestringsWithCoords[key]['matsim_linkID'] = feature['properties']['matsim_linkID']
        linestringsWithCoords[key]['oneway'] = bool(int( feature['properties']['oneway'] ))  # convert value of feature['properties']['oneway'] to a Boolean: "1" becomes 'True' (and also "2", "3", "-1", etc.) and "0" becomes 'False', but "0.0" and any other string not expressing an integer will raise a ValueError. TODO: decide whether and how the 'oneway' property should affect which solutions are allowed
        linestringsWithCoords[key]['permlanes'] = int( feature['properties']['permlanes'] )
      else:
        print(f"A LineString with coordinates {key} already exists in linestringsWithCoords")
        numDuplicateLinestrings += 1
  print("pointsWithCoords is", pointsWithCoords)
  if numDuplicatePoints > 0:
    pluralIndicator = "s" if numDuplicatePoints > 1 else ""
    print(f"WARNING: {numDuplicatePoints} duplicate point{pluralIndicator} ignored.")
#    raise Exception(f"WARNING: {numDuplicatePoints} duplicate point{pluralIndicator} ignored.")
  print("linestringsWithCoords is", linestringsWithCoords)
  if numDuplicateLinestrings > 0:
    pluralIndicator = "s" if numDuplicateLinestrings > 1 else ""
    print(f"WARNING: {numDuplicateLinestrings} duplicate LineStrings ignored.")
    raise Exception(f"WARNING: {numDuplicateLinestrings} duplicate Linestring{pluralIndicator} ignored.")  # TODO: should probably comment out this line, but will wait for more testing to see whether the exception is ever raised
#  exitnodes = []  # as of 16th April 2021, the previously-defined exit-nodes are an input-parameter to this function ('geoJSONtoAdjacencyMatrix')
  for l in linestringsWithCoords:  # add any missing Points at the ends of LineStrings
    for xy in l:
      if xy not in pointsWithCoords:
        print(f"Endpoint {xy} of LineString {l} is not in pointsWithCoords: adding it")
        pointsWithCoords[ xy ] = currentPointID
#        exitnodes.append( currentPointID )  # include the missing Point as an exit-node: all non-missing-point exit-nodes are already specified as a parameter of this function ('geoJSONtoAdjacencyMatrix')  # TODO: confirm that this is what is generally desired
        exitnodes.add( currentPointID )  # include the missing Point as an exit-node: all non-missing-point exit-nodes are already specified as a parameter of this function ('geoJSONtoAdjacencyMatrix')  # TODO: confirm that this is what is generally desired
        print(f"WARNING: including originally-missing endpoint {xy} of LineString {l} as an exit-node: exitnodes is now {exitnodes}")
        currentPointID += 1
      else:
#        print(f"Point {xy} is in pointsWithCoords")
        pass
  print("After addition of any missing points, pointsWithCoords is", pointsWithCoords)
  N = len(pointsWithCoords)
  adjacency = np.zeros( shape=(N,N) )  # construct adjacency matrix
  pointcoordsbyID = {}
  for (xy, pointID) in pointsWithCoords.items():
    pointcoordsbyID[pointID] = xy
#  linkcoordsbyID = {}  # probably don't need, as 'linestringsWithCoords' contains all information about links
#  if flowfunc == 'xExpnegx':  # TODO: must eventually generalise to allow any (non-monotonic in density) flow-function
  freespeed = {}  # free speed per link, in kilometres per hour
  permlanes = {}  # number of lanes per link
  linklength = {}  # length of each link, in kilometres
  linkcapacity = {}  # capacity of each link, in vehicles per hour
  matsimlinkID = {}  # MATSim link-ID of each link
  oneway = {}  # one-way status of each link, either 'True' or 'False'
  for l in linestringsWithCoords:
    a = pointsWithCoords[ l[0] ]
    b = pointsWithCoords[ l[1] ]
    adjacency[a, b] = 1
    adjacency[b, a] = 1
    freespeed[ frozenset({a,b}) ] = linestringsWithCoords[l]['freespeed'] * 3600/1000  # convert from metres per second to kilometres per hour
    permlanes[ frozenset({a,b}) ] = linestringsWithCoords[l]['permlanes']
    linklength[ frozenset({a,b}) ] = linestringsWithCoords[l]['length'] / 1000  # convert from metres to kilometres
    linkcapacity[ frozenset({a,b}) ] = linestringsWithCoords[l]['linkcapacity']
    matsimlinkID[ frozenset({a,b}) ] = linestringsWithCoords[l]['matsim_linkID']
    oneway[ frozenset({a,b}) ] = linestringsWithCoords[l]['oneway']
    linestringsWithCoords[l]['endpointNodeIDs'] = [a, b]  # add to each link in 'linestringsWithCoords' a list of its Points' node-IDs
  print("linestringsWithCoords is", linestringsWithCoords)
  return (adjacency, exitnodes, pointcoordsbyID, linestringsWithCoords, freespeed, permlanes, linklength, linkcapacity, matsimlinkID, oneway)


#def setuptrafficflowproblem( JSONnetworkfilename, b0 ):
def setuptrafficflowproblem( JSONnetworkfilename, exitnodes ):
  terminatorPressure = 0.0  # pressure h at every exit-node, which can be either a leaf-node or not
  eps = 1e-12
#  eps_flowbalance = 1e-3  # James Hilton's value

#  flowfunc = 'HazenWilliams'
#  flowfunc = 'concaveQuadratic'
#  flowfunc = 'xExpnegx'
  flowfunc = 'triangular'
  print("flowfunc is %s" % flowfunc)
#  if flowfunc == 'HazenWilliams':  # these values aren't correct for all links
##    lK = 1.93394e-7
##    lk = 1.0/1.852
#    lK = 2.0e-7
#    lK = 100.0
#    lk = .5
#  elif flowfunc == 'concaveQuadratic':
##    lK = math.sqrt(2.0e-7)
#    lK = 10.0
#    lk = None
#    hvalmax = 4.0

#  if lk != None:
#    print("lK %.8f, lk %.5f" % (lK, lk))
#  else:
#    print("lK %.8f, lk %s" % (lK, lk))

#  starttime = time.time()

  flowijfunc = {'HazenWilliams': HazenWilliamsflowfunc, 'concaveQuadratic': concaveQuadraticflowfunc, 'xExpnegx': xExpnegxflowfunc, 'triangular': triangularflowfunc}
  derivijfunc = {'HazenWilliams': HazenWilliamsderivfunc, 'concaveQuadratic': concaveQuadraticderivfunc, 'xExpnegx': xExpnegxderivfunc, 'triangular': triangularderivfunc}

#  linklength = {}  # length of each link in kilometres; entries will be added if needed

  meanVehicleLength_Km = .005
  meanMinGapinStoppedTraffic_Km = .004

  (adjacency, exitnodes, pointcoordsbyID, linestringsWithCoords, lK, numlanes, linklength, linkcapacity, matsimlinkID, oneway) = geoJSONtoAdjacencyMatrix( JSONnetworkfilename, exitnodes )
  print("adjacency is", adjacency)
#  leafnodes = np.sum(adjacency, axis=1) == 1  # np.sum(adjacency, axis=1) is row-wise sums of 'adjacency'
#  print("leafnodes is", leafnodes)
#  leafnodeindices = [i for (i, x) in enumerate(leafnodes) if x]
#  print("leafnodeindices is", leafnodeindices)
  print("exitnodes is", exitnodes)
  print("pointcoordsbyID is", pointcoordsbyID)
#  b0 = np.array([400])  # the inflow at each (non-leaf) node is proportional to the population there, and is measured in vehicles per hour
  print("lK is", lK)
  print("numlanes is", numlanes)
  print("linklength is", linklength)
  print("linkcapacity is", linkcapacity)
  print("matsimlinkID is", matsimlinkID)
  print("oneway is", oneway)
  if flowfunc == 'xExpnegx':
#    CALCULATEPHI = True  # use the original method that applies when links' capacities are not known in advance
    CALCULATEPHI = False  # for when the links' capacities are known in advance, use the capacity to calculate that link's value for phi
    phi = {}
    tau = {}
    for edge in numlanes.keys():
#      linkcapacity = {}  # TODO: is this line needed?
      if CALCULATEPHI:
        phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
        tau[edge] = phi[edge] * numlanes[edge]
        linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
      else:
        phi[edge] = linkcapacity[edge] / (lK[edge] * numlanes[edge] * math.exp(-1))
        tau[edge] = phi[edge] * numlanes[edge]
      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}

  elif flowfunc == 'concaveQuadratic':
##    hvalmax = {frozenset({0,1}): 3.0, frozenset({1,2}): 2.0, frozenset({2,3}): 2.0, frozenset({2,4}): 1.0}
#    hvalmax = {}
##    lK = {frozenset({0,1}): 4.0, frozenset({1,2}): 8.0, frozenset({2,3}): 8.0, frozenset({2,4}): 4.0}
    tau = {}
#    linkcapacity = {}
    maximumdensityperlane = {}
    for edge in lK.keys():
#      linkcapacity[edge] = lK[edge] * hvalmax[edge] / 4.0
#      hvalmax[edge] = 4.0 * linkcapacity[edge] / lK[edge]
      maximumdensityperlane[edge] = 4.0 * linkcapacity[edge] / (lK[edge] * numlanes[edge])
#      tau[edge] = hvalmax[edge] / 2
      tau[edge] = numlanes[edge] * maximumdensityperlane[edge] / 2
#    argsflowfunc = {'hvalmax': hvalmax, 'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}
    argsflowfunc = {'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength, 'numlanes': numlanes, 'maximumdensityperlane': maximumdensityperlane}

  elif flowfunc == 'triangular':
#    rhoC = {}  # density per lane at capacity
    tau = {}  # in terms of the notation on p. 94 of Treiber and Kesting, tau is numlanes * rho_C
#    maximumdensityperlane = {}
    maximumdensityperlane = 120.0  # in vehicles per kilometre per lane; a constant over all edges - this value taken from Table 8.1, p. 93, 'Traffic Flow Dynamics' (2013) by Treiber and Kesting
    print(f"maximumdensityperlane is {maximumdensityperlane:.4f}")
#    T = 1.4 / 3600  # in hours; a constant over all edges - this value taken from p. 93 of Treiber and Kesting
#    print(f"T is {T:.9f}")
    T = {}
    for edge in lK.keys():
      # N.B. 'linkcapacity' is capacity flow over all lanes, NOT per lane. Similarly, 'tau' is the density over all lanes at capacity flow.
      T[edge] = numlanes[edge] / linkcapacity[edge] - 1/(lK[edge]*maximumdensityperlane)  # in hours; from eq. (8.15), p. 93 of Treiber and Kesting
#      rhoC[edge] = 1 / (lK[edge] * T + 1/maximumdensityperlane)  # density per lane at capacity; from p. 94 of Treiber and Kesting
#      tau[edge] = numlanes[edge] / (lK[edge] * T[edge] + 1/maximumdensityperlane)  # density over entire link-width at capacity; from eq. (8.16), p. 94 of Treiber and Kesting
      tau[edge] = linkcapacity[edge] / lK[edge]  # density over entire link-width at capacity; from eq. (8.13), p. 93 of Treiber and Kesting, substituting Q=C/I and rho_{free} = tau/I
      print(f"For edge {edge}, slope of free-flow part of triangular fundamental diagram is {linkcapacity[edge]/tau[edge]}; slope of congestion part of triangular fundamental diagram is {-linkcapacity[edge]/(maximumdensityperlane*numlanes[edge]-tau[edge])}.")
#    argsflowfunc = {'hvalmax': hvalmax, 'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}
    argsflowfunc = {'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength, 'numlanes': numlanes, 'maximumdensityperlane': maximumdensityperlane, 'T': T}



#  networktype = 'linear'


#  adjacency = np.array([[0, 1], [1, 0]])  # adjacency matrix of linear road-network with two nodes
#  print("adjacency is", adjacency)
#  exitnodes = [1]
##  b0 = np.array([1200])  # the inflow at each (non-leaf) node is proportional to the population there, and is measured in vehicles per hour
#  b0 = np.array([400])  # the inflow at each (non-leaf) node is proportional to the population there, and is measured in vehicles per hour
#  v0 = np.array([10])  # assumed speed of injected traffic at the injection-node, in kilometres per hour
#  if flowfunc == 'xExpnegx':
##    lK = {frozenset({0,1}): 120}  # free speed per link, in kilometres per hour: this contains one highway
#    lK = {frozenset({0,1}): 50}  # free speed per link, in kilometres per hour: this contains one highway
##    numlanes = {frozenset({0,1}): 4}  # number of lanes per link; setting 'numlanes' to a small value such as .01 makes that edge a 'dummy' edge, which is useful because each injection-node must have edge-degree >= 2
#    numlanes = {frozenset({0,1}): 1}  # number of lanes per link; setting 'numlanes' to a small value such as .01 makes that edge a 'dummy' edge, which is useful because each injection-node must have edge-degree >= 2
#    linklength = {frozenset({0,1}): 50}  # length of each link in kilometres for the false-negative example
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}

#  adjacency = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])  # linear network with three nodes
#  exitnodes = [2]
##  exitnodes = [1]
###  b0 = np.array([60/.055])  # population-inflow (demand) at the non-leaf node, in vehicles per hour
##  b0 = np.array([3000])  # population-inflow (demand) at the non-leaf node, in vehicles per hour - yields a false positive
##  b0 = np.array([6000])  # population-inflow (demand) at the non-leaf node, in vehicles per hour - yields a false negative
##  b0 = np.array([500, 0])  # population-inflow (demand) at the non-leaf node, in vehicles per hour
#  b0 = np.array([0, 500])  # population-inflow (demand) at the non-leaf node, in vehicles per hour
##  v0 = np.array([50])  # assumed speed of injected traffic at the injection-node, in kilometres per hour - TODO: is this value realistic? Does the ABM have values for speeds of injected traffic?
#  v0 = np.array([10, 10])  # assumed speed of injected traffic at the injection-nodes, in kilometres per hour
#  if flowfunc == 'xExpnegx':
#    lK = {frozenset({0,1}): 120, frozenset({1,2}): 50}  # free speed per link, in kilometres per hour, for the false-positive example: this contains one highway and one minor road
#    numlanes = {frozenset({0,1}): 4, frozenset({1,2}): 1}  # number of lanes per link for the false-positive example; setting 'numlanes' to a small value such as .01 makes that edge a 'dummy' edge, which is useful because each injection-node must have edge-degree >= 2  # link-capacities are about 3969.4 and 505.6, respectively  #: no solution with zero total congestion is found when minimising 'fnormsquared'  #and still no solution with zero total congestion is found when minimising 'fnormplustotalcongestion'
##    lK = {frozenset({0,1}): 120, frozenset({1,2}): 120}  # free speed per link, in kilometres per hour, for the false-negative example: this contains two highway links that are identical except for the lefthand link being much shorter than the righthand one (but the links' lengths aren't taken into account by the current traffic-model flow)
##    numlanes = {frozenset({0,1}): 4, frozenset({1,2}): 4}  # number of lanes per link for the false-negative example
##    linklength = {frozenset({0,1}): 10, frozenset({1,2}): 20}  # length of each link in kilometres for the false-negative example
#    linklength = {frozenset({0,1}): 50, frozenset({1,2}): 20}  # length of each link in kilometres for the false-negative example
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
###      phi[edge] = omega1_carsperkmperlane / math.log(lK[edge]/1)  # first attempt at calculating phi, which is probably incorrect
##      if lK[edge] == 50:
##        phi[edge] = 27.488  # TODO: implement in Python the Octave code that was used to calculate these two values, so that we'll have access to a value of phi for any arbitrary value of lK
##      elif lK[edge] == 100:
##        phi[edge] = 23.362
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}


#  adjacency = np.array([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]])  # linear network with four nodes (outdated info: on which false positives are found, i.e. when minimising 'fnormplustotalcongestion' no soln is found with zero total congestion)
##  exitnodes = [2]
##  exitnodes = [3]
#  exitnodes = [0, 3]  # exits are the leaf-nodes, although this is no longer required
##  b0 = np.array([15, 0])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
##  b0 = np.array([60/.055, 0])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
##  b0 = np.array([1200, 200])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
#  b0 = np.array([300, 200])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
##  b0 = np.array([500, 0, 0])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
##  b0 = np.array([0, 500, 0])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
##  v0 = np.array([10, 10, 10])  # assumed speed of injected traffic at the injection-nodes, in kilometres per hour
#  v0 = np.array([50, 50, 50])  # assumed speed of injected traffic at the injection-nodes, in kilometres per hour
#  if flowfunc == 'xExpnegx':
#    lK = {frozenset({0,1}): 120, frozenset({1,2}): 50, frozenset({2,3}): 50}  # free speed per link, in kilometres per hour
##    numlanes = {frozenset({0,1}): 4, frozenset({1,2}): 1, frozenset({2,3}): 1}  # link-capacities are about 40.9, 5.1, and 5.1, respectively  #: still no solution with zero total congestion is found when minimising 'fnormplustotalcongestion')
#    numlanes = {frozenset({0,1}): 2, frozenset({1,2}): 1, frozenset({2,3}): 1}  # link-capacities are about 20.4, 5.1, and 5.1, respectively
##    linklength = {frozenset({0,1}): 10, frozenset({1,2}): 20, frozenset({2,3}): 20}  # length of each link in kilometres for the false-positive example
#    linklength = {frozenset({0,1}): 100, frozenset({1,2}): 20, frozenset({2,3}): 20}  # length of each link in kilometres for the false-positive example
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}


#  adjacency = np.array([[0, 1, 0, 0, 0], [1, 0, 1, 0, 0], [0, 1, 0, 1, 0], [0, 0, 1, 0, 1], [0, 0, 0, 1, 0]])  # linear network with five nodes
##  b0 = np.array([5.685, 7.58, 9.00125])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
#  b0 = np.array([-5.685, 7.58, 9.00125])  # try with a negative inflow, ie. outflow, at one non-leaf node (measured in vehicles per hour)



#  networktype = 'highway'
#
#  adjacency = np.array([[0, 1, 0, 0, 0], [1, 0, 1, 0, 0], [0, 1, 0, 1, 1], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0]])  # 'highway' network with five nodes, four nodes on a high-capacity highway and a fifth node adjacent to the second of the highway's two non-terminal nodes
#  exitnodes = [3, 4]
###  b0 = np.array([3.5, 0.0])  # 4 solutions found from 11^2 starting-points with 0 <= hi <= 10
##  b0 = np.array([60/.055, 0])
##  b0 = np.array([5000, 100, 200])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
#  b0 = np.array([1000, 100, 200])  # population-inflow (demand) at the non-leaf nodes, in vehicles per hour
#  v0 = np.array([50, 50, 50])  # assumed speed of injected traffic at the injection-nodes, in kilometres per hour
#  if flowfunc == 'concaveQuadratic':
#    hvalmax = {frozenset({0,1}): 3.0, frozenset({1,2}): 2.0, frozenset({2,3}): 2.0, frozenset({2,4}): 1.0}
#    lK = {frozenset({0,1}): 4.0, frozenset({1,2}): 8.0, frozenset({2,3}): 8.0, frozenset({2,4}): 4.0}
#    tau = {}
#    for edge in lK.keys():
#      tau[edge] = hvalmax[edge] / 2
#    linkcapacity = {}
#    for edge in lK.keys():
#      linkcapacity[edge] = lK[edge] * hvalmax[edge] / 4.0
#    argsflowfunc = {'hvalmax': hvalmax, 'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity}
#  elif flowfunc == 'HazenWilliams':
#    lK = 100.0
#    lk = .5
#    argsflowfunc = {'lK': lK, 'lk': lk}
#  elif flowfunc == 'xExpnegx':
##    lK = {frozenset({0,1}): 3*math.exp(1), frozenset({1,2}): 4*math.exp(1), frozenset({2,3}): 4*math.exp(1), frozenset({2,4}): 1*math.exp(1)}  # values lK chosen to make link-capacities 3, 4, 4, and 1, respectively
#    lK = {frozenset({0,1}): 120, frozenset({1,2}): 140, frozenset({2,3}): 140, frozenset({2,4}): 50}  # free speed per link, in kilometres per hour
##    numlanes = {frozenset({0,1}): 3*math.exp(1), frozenset({1,2}): 4*math.exp(1), frozenset({2,3}): 4*math.exp(1), frozenset({2,4}): 1*math.exp(1)}  # values numlanes chosen to make link-capacities 3, 4, 4, and 1, respectively
#    numlanes = {frozenset({0,1}): 2, frozenset({1,2}): 3, frozenset({2,3}): 3, frozenset({2,4}): 1}  # link-capacities are about 20.4, 5.1, and 5.1, respectively
#    linklength = {frozenset({0,1}): 10, frozenset({1,2}): 10, frozenset({2,3}): 10, frozenset({2,4}): 20}  # length of each link in kilometres
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}



#  networktype = 'star'

#  adjacency = np.array([[0, 1, 0, 0], [1, 0, 1, 1], [0, 1, 0, 0], [0, 1, 0, 0]])  # y-shaped network with four nodes
#  b0 = np.array([15.0])  # population-inflow (demand) at the non-leaf node, in vehicles per hour
#  if flowfunc == 'concaveQuadratic':
#    hvalmax = {frozenset({0,1}): 2.0, frozenset({1,2}): 1.0, frozenset({1,3}): 6.0}

#  adjacency = np.array([[0, 1, 0, 0, 0], [1, 0, 1, 0, 0], [0, 1, 0, 1, 1], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0]])  # y-shaped network with five nodes
##  b0 = np.array([2.0, 8.0])  # population-inflows (demands) at the non-leaf nodes, in vehicles per hour
#  b0 = np.array([2.0, 26.673])  # population-inflows (demands) at the non-leaf nodes, in vehicles per hour
#  if flowfunc == 'xExpnegx':
#    lK = {frozenset({0,1}): 120, frozenset({1,2}): 50, frozenset({2,3}): 50, frozenset({2,4}): 50}  # free speed per link, in kilometres per hour
#    numlanes = {frozenset({0,1}): 4, frozenset({1,2}): 1, frozenset({2,3}): 1, frozenset({2,4}): 1}  # link-capacities are about ?.? and ?.?, respectively #: no solution with zero total congestion is found when minimising 'fnormsquared'  #and still no solution with zero total congestion is found when minimising 'fnormplustotalcongestion'
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity}

#  adjacency = np.array([[0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0], [0, 0, 0, 1, 1, 0], [0, 0, 1, 0, 0, 0], [0, 1, 1, 0, 0, 1], [0, 0, 0, 0, 1, 0]])  # m-shaped network with six nodes
##  b0 = np.array([2000.0, 2600.0, 0.0])  # population-inflows (demands) at the non-leaf nodes, in vehicles per hour
#  b0 = np.array([200.0, 260.0, 0.0])  # population-inflows (demands) at the non-leaf nodes, in vehicles per hour
#  if flowfunc == 'xExpnegx':
#    lK = {frozenset({0,1}): 50, frozenset({1,4}): 120, frozenset({2,3}): 50, frozenset({2,4}): 120, frozenset({4,5}): 50}  # free speed per link, in kilometres per hour
##    lK = {frozenset({0,1}): 3, frozenset({1,4}): 120, frozenset({2,3}): 3, frozenset({2,4}): 120, frozenset({4,5}): 50}  # free speed per link, in kilometres per hour
##    numlanes = {frozenset({0,1}): 1, frozenset({1,4}): 4, frozenset({2,3}): 1, frozenset({2,4}): 4, frozenset({4,5}): 1}  # number of lanes per link
#    numlanes = {frozenset({0,1}): 0.01, frozenset({1,4}): 4, frozenset({2,3}): 0.01, frozenset({2,4}): 4, frozenset({4,5}): 1}  # number of lanes per link; setting 'numlanes' to a small value such as .01 makes that edge a 'dummy' edge, which is useful because each injection-node must have edge-degree >= 2
#    linklength = {frozenset({0,1}): 1, frozenset({1,4}): 10, frozenset({2,3}): 1, frozenset({2,4}): 10, frozenset({4,5}): 20}  # length of each link in kilometres
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}



#  networktype = 'nontree'  # for testing shortest-path approach
#
#  adjacency = np.array([[0, 1, 1, 0, 0, 0], [1, 0, 0, 0, 1, 0], [1, 0, 0, 1, 0, 0], [0, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 1, 0]])  # non-tree network with six nodes
#  exitnodes = [5]
##  b0 = np.array([200, 0, 0, 0, 0])  # population-inflows (demands) at the non-leaf nodes, in vehicles per hour - an inflow b0=200 seems unrealistically low.
##  b0 = np.array([506, 0, 0, 0, 0])  # b0=506.0 is the highest inflow for which the global optimiser has found a flow-balanced solution.
##  b0 = np.array([507, 0, 0, 0, 0])  # b0=507.0 is the lowest inflow for which the global optimiser has failed to find a flow-balanced solution.
##  b0 = np.array([407, 0, 0, 100, 0])
#  b0 = np.array([1900, 0, 0, 100, 0])
##  b0 = np.array([1200, 0, 0, 0, 0])  # an inflow b0=1200 seems realistic
#  v0 = np.array([50, 50, 50, 50, 50])  # assumed speed of injected traffic at the injection-nodes, in kilometres per hour
#  if flowfunc == 'xExpnegx':
#    lK = {frozenset({0,1}): 120, frozenset({1,4}): 120, frozenset({0,2}): 20, frozenset({2,3}): 50, frozenset({3,4}): 50, frozenset({4,5}): 50}  # free speed per link, in kilometres per hour
#    numlanes = {frozenset({0,1}): 2, frozenset({1,4}): 2, frozenset({0,2}): 1, frozenset({2,3}): 1, frozenset({3,4}): 1, frozenset({4,5}): 1}  # number of lanes per link
#    linklength = {frozenset({0,1}): 6, frozenset({1,4}): 6, frozenset({0,2}): 4, frozenset({2,3}): 4, frozenset({3,4}): 4, frozenset({4,5}): 5}  # length of each link in kilometres
#    phi = {}
#    tau = {}
#    linkcapacity = {}
#    for edge in numlanes.keys():
#      phi[edge] = phiOnLinkrs_inVehiclesPerKmPerLane( lK[edge], meanMinGapinStoppedTraffic_Km, meanVehicleLength_Km )
#      tau[edge] = phi[edge] * numlanes[edge]
#      linkcapacity[edge] = lK[edge] * phi[edge] * numlanes[edge] * math.exp(-1)
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
#    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}


# (End of network definitions)



#  if flowfunc == 'concaveQuadratic':
#    print("hvalmax is %s" % hvalmax)
#    print("lK is %s" % lK)
  for key in argsflowfunc.keys():
    print("argsflowfunc['%s'] is %s" % (key, argsflowfunc[key]))
#  print("adjacency is", adjacency)
  N = adjacency.shape[0]  # number of nodes in road-network
  print("N is %d" % N)
##  h = np.zeros(N)  #h = np.zeros(shape=(N,1))
##  print("np.sum(adjacency, axis=1) == 1 is", np.sum(adjacency, axis=1) == 1)
##  print("np.where( np.sum(adjacency, axis=1) == 1 ) is", np.where( np.sum(adjacency, axis=1) == 1 ))
#  leafnodesBool = np.sum(adjacency, axis=1) == 1  # np.sum(adjacency, axis=1) is row-wise sums of 'adjacency'
#  print("leafnodesBool is", leafnodesBool)

##  exitnodes = np.where( np.sum(adjacency, axis=1) == 1 ).tolist()  # np.sum(adjacency, axis=1) gives row-wise sum of 'adjacency'
#  exitnodes = [i for i in range(N) if leafnodesBool[i]]
  print("exitnodes is", exitnodes)
#  injectionnodes = [i for i in range(N) if not leafnodesBool[i]]
  injectionnodes = [i for i in range(N) if i not in exitnodes]
  print("injectionnodes is", injectionnodes)
  indexamonginjectionorexitnodes = np.zeros(N, dtype=int)
  for k, i in enumerate(list(injectionnodes)):
#    print("k %d, i %d" % (k, i))
    indexamonginjectionorexitnodes[i] = k
  for k, i in enumerate(list(exitnodes)):
#    print("k %d, i %d" % (k, i))
    indexamonginjectionorexitnodes[i] = k
  print("indexamonginjectionorexitnodes is %s" % indexamonginjectionorexitnodes)
  neighboursofj = {}
  for j in range(N):
    neighboursofj[j] = [i for i in range(N) if adjacency[i,j] == 1]
#    print("Node %d has neighbours %s" % (j, neighboursofj[j]))
  print("neighboursofj is %s" % neighboursofj)
#  h[exitnodes] = terminatorPressure
  hExitNodes = terminatorPressure*np.ones(len(exitnodes))
#  print("h is", h)
  print("hExitNodes is", hExitNodes)
  assert len(injectionnodes) > 0
  assert len(injectionnodes) < N
##  iteration = 0
#  numcallsoffuncf = 0
####  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, networktype, adjacency, b0, N, leafnodesBool, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0)
###  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, networktype, adjacency, b0, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0)
##  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0, pointcoordsbyID)
#  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, argsflowfunc, pointcoordsbyID)
  return (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, argsflowfunc, pointcoordsbyID, linestringsWithCoords)


def sign(x):
#  print("sign(): x is %s, or %.10f" % (repr(x), x))
#  return (x > 0) - (x < 0)  # gives 'TypeError: numpy boolean subtract, the `-` operator, is not supported, use the bitwise_xor, the `^` operator, or the logical_xor function instead.'
  if x == 0.0:
    return 0
  else:
    return math.copysign(1, x)

#def HazenWilliamsflowfunc(i, j, hi, hj, argsflowfunc):
def HazenWilliamsflowfunc(i, j, hiMinushj, argsflowfunc):
#  hval = abs(hi-hj)
  hval = abs(hiMinushj)
#  hsign = sign(hi - hj)
  hsign = sign(hiMinushj)
  lK = argsflowfunc['lK']
  lk = argsflowfunc['lk']
  return hsign*(hval*lK)**lk
def HazenWilliamsderivfunc(hval, i, j, argsflowfunc):
  lK = argsflowfunc['lK']
  lk = argsflowfunc['lk']
  return lk * lK**lk * (hval+eps)**(lk-1.0)
#def logflowfunc(i, j, hi, hj, argsflowfunc):
def logflowfunc(i, j, hiMinushj, argsflowfunc):
#  hval = abs(hi-hj)
  hval = abs(hiMinushj)
#  hsign = sign(hi - hj)
  hsign = sign(hiMinushj)
  lK = argsflowfunc['lK']
  lk = argsflowfunc['lk']
  return hsign*lK*math.log(lk*hval+1.0)

#def concaveQuadraticflowfunc(i, j, hi, hj, argsflowfunc):
def concaveQuadraticflowfunc(i, j, hiMinushj, argsflowfunc):
#  hval = abs(hi-hj)
  hval = abs(hiMinushj)
#  hvalmax = argsflowfunc['hvalmax']
##  if isinstance(hvalmax, float):
##    hvalmax_current = hvalmax
##  elif isinstance(hvalmax, dict):
##    hvalmax_current = hvalmax[frozenset({i,j})]
#  hvalmax_current = hvalmax[frozenset({i,j})]
  numlanes = argsflowfunc['numlanes']
  numlanes_current = numlanes[frozenset({i,j})]
  maximumdensityperlane = argsflowfunc['maximumdensityperlane']
  maximumdensityperlane_current = maximumdensityperlane[frozenset({i,j})]
#  if hval > hvalmax_current:
  if hval > numlanes_current*maximumdensityperlane_current:
    return 0.0  # ensure that sign of flow q isn't opposite to (i.e., the negative of) sign of hi-hj: want non-negative flow for positive hi-hj, and non-positive flow for negative hi-hj
#  hsign = sign(hi - hj)
  hsign = sign(hiMinushj)
  lK = argsflowfunc['lK']
#  if isinstance(lK, float):
#    lK_current = lK
#  elif isinstance(lK, dict):
#    lK_current = lK[frozenset({i,j})]
  lK_current = lK[frozenset({i,j})]
#  return hsign*hval*lK_current*(1-hval/hvalmax_current)  # maximum |q| is lK*hvalmax_current/4
  return hsign*hval*lK_current*(1-hval/(numlanes_current*maximumdensityperlane_current))  # maximum |q| is lK*numlanes_current*maximumdensityperlane_current/4

def concaveQuadraticderivfunc(hval, i, j, argsflowfunc):
#  hvalmax = argsflowfunc['hvalmax']
##  if isinstance(hvalmax, float):
##    hvalmax_current = hvalmax
##  elif isinstance(hvalmax, dict):
##    hvalmax_current = hvalmax[frozenset({i,j})]
#  hvalmax_current = hvalmax[frozenset({i,j})]
  lK = argsflowfunc['lK']
#  if isinstance(lK, float):
#      lK_current = lK
#    elif isinstance(lK, dict):
#    lK_current = lK[frozenset({i,j})]
  lK_current = lK[frozenset({i,j})]
#  return lK_current * (1 - 2*hval/hvalmax_current)
  numlanes = argsflowfunc['numlanes']
  numlanes_current = numlanes[frozenset({i,j})]
  maximumdensityperlane = argsflowfunc['maximumdensityperlane']
  maximumdensityperlane_current = maximumdensityperlane[frozenset({i,j})]
  return lK_current * (1 - 2*hval/(numlanes_current*maximumdensityperlane_current))

#def xExpnegxflowfunc(i, j, hi, hj, argsflowfunc):
def xExpnegxflowfunc(i, j, hiMinushj, argsflowfunc):
#  print("i %d, j %d, hi %.10f, hj %.10f" % (i, j, hi, hj))
#  hval = abs(hi-hj)
  hval = abs(hiMinushj)
#  hsign = sign(hi - hj)
  hsign = sign(hiMinushj)
  lK = argsflowfunc['lK']
#  print("xExpnegxflowfunc(): lK is %s" % lK)
  if isinstance(lK, float):
    lK_current = lK
  elif isinstance(lK, dict):
    lK_current = lK[frozenset({i,j})]
#  print("xExpnegxflowfunc(): lK_current is %s" % lK_current)
  numlanes = argsflowfunc['numlanes']
#  print("xExpnegxflowfunc(): numlanes is %s" % numlanes)
  if isinstance(numlanes, float):
    numlanes_current = numlanes
  elif isinstance(numlanes, dict):
    numlanes_current = numlanes[frozenset({i,j})]
  phi = argsflowfunc['phi']
  if isinstance(phi, float):
    phi_current = phi
  elif isinstance(phi, dict):
    phi_current = phi[frozenset({i,j})]
  return hsign*hval*lK_current*math.exp(-hval/(phi_current*numlanes_current))  # maximum |q| is lK*numlanes*exp(-1)

def xExpnegxderivfunc(hval, i, j, argsflowfunc):
  lK = argsflowfunc['lK']
  if isinstance(lK, float):
    lK_current = lK
  elif isinstance(lK, dict):
    lK_current = lK[frozenset({i,j})]
  numlanes = argsflowfunc['numlanes']
  if isinstance(numlanes, float):
    numlanes_current = numlanes
  elif isinstance(numlanes, dict):
    numlanes_current = numlanes[frozenset({i,j})]
  phi = argsflowfunc['phi']
  if isinstance(phi, float):
    phi_current = phi
  elif isinstance(phi, dict):
    phi_current = phi[frozenset({i,j})]
  return lK_current * (1 - hval/(phi_current*numlanes_current)) * math.exp(-hval / (phi_current*numlanes_current))

def triangularflowfunc(i, j, hiMinushj, argsflowfunc):
  hval = abs(hiMinushj)
  numlanes = argsflowfunc['numlanes']
  numlanes_current = numlanes[frozenset({i,j})]
  maximumdensityperlane = argsflowfunc['maximumdensityperlane']  # a constant over all edges
#  maximumdensityperlane_current = maximumdensityperlane[frozenset({i,j})]
#  if hval > hvalmax_current:
  if hval > numlanes_current*maximumdensityperlane:
    return 0.0  # ensure that sign of flow q isn't opposite to (i.e., the negative of) sign of hi-hj: want non-negative flow for positive hi-hj, and non-positive flow for negative hi-hj
  hsign = sign(hiMinushj)
  lK = argsflowfunc['lK']
  lK_current = lK[frozenset({i,j})]
#  T = argsflowfunc['T']  # a constant over all edges
  T = argsflowfunc['T']
  T_current = T[frozenset({i,j})]
  tau = argsflowfunc['tau']
  tau_current = tau[frozenset({i,j})]
##  return hsign*hval*lK_current*(1-hval/hvalmax_current)  # maximum |q| is lK*hvalmax_current/4
#  return hsign*hval*lK_current*(1-hval/(numlanes_current*maximumdensityperlane_current))  # maximum |q| is lK*numlanes_current*maximumdensityperlane_current/4
#  if hval <= numlanes_current * rhoC_current:
  if hval <= tau_current:
#    return hsign*lK_current*hval*numlanes_current  # maximum |q| is numlanes_current*lK/(lK*T+1/maximumdensityperlane)  # from p. 93 of Treiber and Kesting
    return hsign*lK_current*hval  # maximum |q| is numlanes_current*lK/(lK*T+1/maximumdensityperlane)  # from p. 93 of Treiber and Kesting
#  else:  # hval > numlanes_current * rhoC_current
  else:  # hval > tau_current
##    return numlanes_current/T*(1-hsign*hval/maximumdensityperlane)
#    return (1/T)*(1-hsign*hval/(numlanes_current*maximumdensityperlane))
    return (1/T_current)*(1-hsign*hval/(numlanes_current*maximumdensityperlane))

def triangularderivfunc(hval, i, j, argsflowfunc):
  lK = argsflowfunc['lK']
  lK_current = lK[frozenset({i,j})]
#  return lK_current * (1 - 2*hval/hvalmax_current)
  numlanes = argsflowfunc['numlanes']
  numlanes_current = numlanes[frozenset({i,j})]
  maximumdensityperlane = argsflowfunc['maximumdensityperlane']
#  T = argsflowfunc['T']  # a constant over all edges
  T = argsflowfunc['T']
  T_current = T[frozenset({i,j})]
#  if hval <= numlanes_current * rhoC_current:
  if hval <= tau_current:
#    return lK_current * numlanes_current
    return lK_current
  elif tau_current < hval <= numlanes_current*maximumdensityperlane:
##    return -numlanes_current/(T*maximumdensityperlane)
#    return -1/(T*numlanes_current*maximumdensityperlane)
    return -1/(T_current*numlanes_current*maximumdensityperlane)
  # TODO: return what if hval > numlanes_current*maximumdensityperlane ?



def debug(debuggingMessage, DEBUG = False):
#  DEBUG = False  # flag whether to print debugging-messages
  if DEBUG:
    print(debuggingMessage)

def f(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  debug("f(): h is %s, type(h) is %s" % (h, type(h)))
  fvect = []
  for j in injectionnodes:  # construct Sum_i q_{ij} = b_j for injection-nodes j
    totalinflowj = 0.0
    hj = h[indexamonginjectionorexitnodes[j]]
    for i in neighboursofj[j]:
      if i in exitnodes:
        hi = hExitNodes[indexamonginjectionorexitnodes[i]]
      else:
#        print("h is %s" % h)
        hi = h[indexamonginjectionorexitnodes[i]]
#      debug("i %d, j %d, hi %.10f, hj %.10f" % (i, j, hi, hj))
      debug("i %d, j %d" % (i, j))
      debug("hi %.10f, hj %.10f" % (hi, hj))
#      qij = flowijfunc[flowfunc](i, j, hi, hj, argsflowfunc)
      qij = flowijfunc[flowfunc](i, j, hi-hj, argsflowfunc)
      debug("qij is %.5f" % qij)
      totalinflowj += qij
    totalinflowj += b0[indexamonginjectionorexitnodes[j]]
    fvect.append( totalinflowj )
#    print("fvect[%d] is" % j, fvect[j])
    debug("totalinflowj at (injection-)node %d is %.10f" % (j, totalinflowj))
  debug("fvect is %s" % fvect)
  return np.array(fvect)

def fnormsquared(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  normsquared = scipy.linalg.norm( f(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )**2
#  debug("normsquared is %.10f" % normsquared)
  return normsquared

def fnorm(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  norm = scipy.linalg.norm( f(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )
#  debug("norm is %.10f" % norm)
  return norm

#  congestiononlink = linkcongestion(r, s, qrs, abs(hr-hs), argsflowfunc)
#  congestiononlink = linkcongestion(p, j, flowenteringparent, densityinlink, argsflowfunc)
def linkcongestion(r, s, flow, density, argsflowfunc):
  # The current measure of congestion, which is zero if the link's density is <= density at capacity and is the link's capacity less its actual flow otherwise, is problematic; for example, a link into which a congestion-zone has flowed upsttream might have very low flow only because its head-node has very low demand (i.e. injection-flow), and this will be measured to be high congestion. A better measure might follow from assessing congestion not on the links but at the injection-nodes, taking it to be the demand (injection-flow) less the total inflow into links of the network (noting that total inflow will be less than demand where one or more of the adjacent links has had its flow reduced by an upstream-propagating congestion-zone).
  print(f"linkcongestion(): calculating congestion on link {frozenset({r,s})} from flow {flow} and density {density}...")
  humpdensity = argsflowfunc['tau'][frozenset({r,s})]
#  print("tau[%s] is %.4f" % ((r,s), argsflowfunc['tau'][(r,s)]))
  print("tau[%s] is %.4f" % ((r,s), argsflowfunc['tau'][frozenset({r,s})]))
##  if abs(hr-hs) >= argsflowfunc['tau'][edge]:
#  if density >= argsflowfunc['tau'][(r,s)]:
  if density >= humpdensity:
###    flow = flowijfunc[flowfunc](r, s, hr-hs, argsflowfunc)
##    flow = flowijfunc[flowfunc](r, s, density, argsflowfunc)
#    linkcapacity = argsflowfunc['linkcapacity'][(r,s)]
    linkcapacity = argsflowfunc['linkcapacity'][frozenset({r,s})]
    congestiononlink = linkcapacity - abs(flow)
  else:
    congestiononlink = 0.0
#  print("congestiononlink is %.4f" % congestiononlink)
  return congestiononlink

def congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
#  print("hInjectionNodes is %s" % hInjectionNodes)
  congestionij = {}
  totalcongestion = 0.0
#  print("indexamonginjectionorexitnodes is %s" % indexamonginjectionorexitnodes)
#  if flowfunc == 'concaveQuadratic':
#    edges = argsflowfunc['lK'].keys()
#  elif flowfunc == 'xExpnegx':
#    edges = argsflowfunc['numlanes'].keys()
  if flowfunc in ('concaveQuadratic', 'xExpnegx', 'triangular'):
    edges = argsflowfunc['lK'].keys()
#  print("edges is", edges)
  for edge in edges:
#    print("congestioninnetwork(): edge is %s" % edge)
    (r, s) = edge
    if r in exitnodes:
      hr = hExitNodes[indexamonginjectionorexitnodes[r]]
    else:
      hr = hInjectionNodes[indexamonginjectionorexitnodes[r]]
    if s in exitnodes:
      hs = hExitNodes[indexamonginjectionorexitnodes[s]]
    else:
      hs = hInjectionNodes[indexamonginjectionorexitnodes[s]]
    debug("r is %d, s is %d, hr is %.3f, hs is %.3f" % (r, s, hr, hs))
    if flowfunc == 'concaveQuadratic':
      debug("hvalmax[%s] is %.4f" % (edge, argsflowfunc['hvalmax'][edge]))
    linkcapacity = argsflowfunc['linkcapacity'][edge]
#    print("linkcapacity is %s" % linkcapacity)
    qrs = flowijfunc[flowfunc](r, s, hr-hs, argsflowfunc)
#    print("tau[%s] is %.4f" % (edge, argsflowfunc['tau'][edge]))
#    if abs(hr-hs) >= argsflowfunc['tau'][edge]:
#      qrs = flowijfunc[flowfunc](r, s, hr-hs, argsflowfunc)
#      congestiononlink = linkcapacity - abs(qrs)
##      print("qrs is %.4f, congestiononlink is %.4f" % (qrs, congestiononlink))
#    else:
#      congestiononlink = 0.0
    congestiononlink = linkcongestion(r, s, qrs, abs(hr-hs), argsflowfunc)
#    print("congestiononlink is %.4f" % congestiononlink)
    congestionij[edge] = congestiononlink
    totalcongestion += congestiononlink
#  elif flowfunc == 'xExpnegx':
##    argsflowfunc = {'numlanes': numlanes}
#    for edge in argsflowfunc['numlanes'].keys():
##      numlanes_current = argsflowfunc['numlanes'][frozenset({i,j})]
#      numlanes_current = argsflowfunc['numlanes'][edge]
#      linkcapacity =  argsflowfunc['numlanes'][edge] * math.exp(-1)
#      print("Flow-capacity of link %s is %.3f" % (edge, linkcapacity))
#  return totalcongestion
  return (congestionij, totalcongestion)

def totalcongestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  (congestionij, totalcongestion) = congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
  return totalcongestion

def fnormlesstotalcongestion(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  totalcongestionfactor = .2  # TODO: what should this value be? - consider using Lagrange multipliers
  return fnorm(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) - totalcongestionfactor*totalcongestioninnetwork(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)

def fnormplustotalcongestion(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  totalcongestionfactor = .2  # TODO: what should this value be? - consider using Lagrange multipliers
  return fnorm(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) + totalcongestionfactor*totalcongestioninnetwork(h, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)

eps_GlobalMinconvergence = 7e-8  # empirically obtained value that seems (so far) to yield solutions comparable to those given by Newton-Raphson with eps_NRconvergence = 1e-3 (James Hilton's value)
def printlocalminimum(h, fAtMinimum, accepted):
#  print("Local minimum %.4e at h %s: accepted is %s" % (fAtMinimum, h, accepted))
  # Stop the basinhopping routine if the function-value at this minimum is smaller than a threshold:
  if fAtMinimum < eps_GlobalMinconvergence and accepted: return True

def runsolver(NEWTONRAPHSON, N, injectionnodes, exitnodes, indexamonginjectionorexitnodes, terminatorPressure, hExitNodes, neighboursofj, flowfunction, derivfunction, hJmin, hJmax, hinitInjectionNodes, argsflowfunc, ENERGYFUNCTION, flowijfunc, flowfunc, b0):
  if not NEWTONRAPHSON:  # use a minimum-finding method to find the minimum of ||f(h)||^2; the chosen method seeks either a local or a global minimum, and either uses the Jacobian or is a Jacobian-less method such as the secant method or a hybrid method, perhaps approximating the Jacobian:
#    hInjectionNodes = scipy.optimize.newton(f, np.zeros(len(injectionnodes)), tol=1e-9, maxiter=100)  # uses the secant method if "derivative" 'fprime' of f isn't given, otherwise uses the Newton-Raphson method; however, it's not clear from the documentation of 'newton' that 'fprime' is a Jacobian or that this method can handle multivariate problems. (Untried: if the "second-order derivative" 'fprime2' of f is also given, then uses Halley’s method.) Using 'tol' and 'maxiter', as for example in 'tol=1e-9, maxiter=100', allows a more precise solution to be found, at the cost of running more than the default number (fifty) of iterations.
#    hInjectionNodes = scipy.optimize.fsolve(f, np.zeros(len(injectionnodes)))  # when a function 'fprime' to compute the Jacobian of f isn't specified, the Jacobian is calculated by a forward-difference approximation ('fsolve' is a wrapper around MINPACK’s 'hybrd' algorithm, which finds a zero of a system of N non-linear functions in N variables by a modification of the Powell hybrid method). This 'fsolve' method requires fewer iterations than 'newton' when neither is supplied with an 'fprime' function for the Jacobian/"derivative".

##    hInjectionNodes = scipy.optimize.fmin_bfgs(fnormsquared, np.zeros(len(injectionnodes)))  # getting trapped in local minimum?
#    res = scipy.optimize.minimize(fnormsquared, np.zeros(len(injectionnodes)), tol=1.0)

##    hJ0 = random.choices(range(hJmin, hJmax+1), k=3)  # initial guess at solution
    hJ0 = [0] * len(injectionnodes)  # default (zero-vector) initial guess at solution

    hJ0 = hinitInjectionNodes
    print("Initial guess is h = %s" % repr(hJ0))

    # Specify bounds on hJ (hInjectionNodes), as required by L-BFGS-B:
#    hJmin = int(-4)
#    hJmax = int(4)

#    hJmin = int(-5)  # values used to obtain problematic solution in which h_1 is negative despite h_2 and b_2 being positive
#    hJmax = int(4)
#    hJmin = int(-18)  # values used to obtain a second problematic solution in which h_1 is large and negative; implies h_1 could be arbitrarily large and negative and still yield a solution
#    hJmax = int(5)

#    hJminvals = [hJmin] * 2
    hJminvals = [hJmin] * len(injectionnodes)
    hJmaxvals = [hJmax] * len(injectionnodes)
    bounds = [(low, high) for low, high in zip(hJminvals, hJmaxvals)]
    # Method L-BFGS-B is for bound-constrained minimisation; works when (check this) the problem is smooth and bounded (is it?):
#    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, tol=.001)
#    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, tol=.0001)
#    minimizer_kwargs = dict(method="L-BFGS-B", tol=.0001)  # drop the bounds-constraint: seems to make the global optimiser run in less time.
#    minimizer_kwargs = dict(args=(injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0), method="L-BFGS-B", tol=.0001)
    if ENERGYFUNCTION == 'conservationOfMassResidual':
      tol=.0001
    elif ENERGYFUNCTION == 'conservationOfMassResidualPlusTotalCongestion':
      tol=.01
    elif ENERGYFUNCTION == 'conservationOfMassResidualLessTotalCongestion':
#      tol=.01  # TODO: establish a good value
      pass
    minimizer_kwargs = dict(args=(injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc), method="L-BFGS-B", tol=tol)  # 'gtol' option of L-BFGS-B is set to value of 'tol'
#    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, tol=.00005)
#    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, options=dict(ftol=.001))

##    minimizer_kwargs = dict(method="Powell")  # Powell method is for unconstrained minimisation; takes about seven times as many iterations as L-BFGS-B to converge

##    minimizer_kwargs = dict(method="Nelder-Mead", tol=None)  # Nelder-Mead method uses simplex algorithm somehow.
#    minimizer_kwargs = dict(method="Nelder-Mead", tol=.00001)  # Nelder-Mead method uses simplex algorithm somehow; 'tol' is chosen to be sufficiently small for the solution found to equal the theoretical solution to four decimal places.

#    print("minimizer_kwargs is %s" % minimizer_kwargs)

    if ENERGYFUNCTION == 'conservationOfMassResidual':
#      res = scipy.optimize.basinhopping(fnormsquared, hJ0, minimizer_kwargs=minimizer_kwargs, niter=200)  # run for 200 basin-hopping iterations rather than the default of one hundred
#      res = scipy.optimize.basinhopping(fnormsquared, hJ0, minimizer_kwargs=minimizer_kwargs, niter=400)  # run for larger number of basin-hopping iterations rather than the default of one hundred
#      res = scipy.optimize.basinhopping(fnormsquared, hJ0, minimizer_kwargs=minimizer_kwargs, callback=printlocalminimum, niter=1)
      res = scipy.optimize.basinhopping(fnormsquared, hJ0, minimizer_kwargs=minimizer_kwargs, callback=printlocalminimum)  # minimise deviance from mass-conservation
    elif ENERGYFUNCTION == 'conservationOfMassResidualPlusTotalCongestion':
      res = scipy.optimize.basinhopping(fnormplustotalcongestion, hJ0, minimizer_kwargs=minimizer_kwargs, callback=printlocalminimum)  # minimise deviance from mass-conservation and maximise total congestion in network
    elif ENERGYFUNCTION == 'conservationOfMassResidualLessTotalCongestion':
      res = scipy.optimize.basinhopping(fnormlesstotalcongestion, hJ0, minimizer_kwargs=minimizer_kwargs, callback=printlocalminimum)  # minimise deviance from mass-conservation and maximise total congestion in network

    hInjectionNodes = res.x  # the solution-array returned
    print("res is %s" % res)


  else:  # use re-implementation of Newton-Raphson method in James Hilton's code 'gs_network_flow.cpp':
    # Extreme values of MAXITERS and eps_NRconvergence were tried in an attempt to obtain, starting from h0=(0, 1.7, 3.5, 0), a good approximation to the theoretical second solution (0, 2, 4, 0) for the four-node network with b_2 = 10; it took MAXITERS = and eps_NRconvergence = 1e-x to get to the theoretical solution.
#    MAXITERS = 100  # the values on this line and the next find h=(0, 2.0559, 4.1117, 0) in .003 seconds, 3 iterations
#    eps_NRconvergence = 1e-3
#    MAXITERS = 100  # the values on this line and the next find h=(0, 2.0505, 4.1011, 0) in .014 seconds, 19 iterations
#    eps_NRconvergence = 1e-6
#    MAXITERS = 100  # the values on this line and the next find h=(0, 2.0339, 4.0679, 0) in .074 seconds, 100 iterations
#    eps_NRconvergence = 1e-7
#    MAXITERS = 800  # the values on this line and the next find h=(0, 2.0088, 4.0176, 0) in .56 seconds, 785 iterations
#    eps_NRconvergence = 1e-9
#    MAXITERS = 32000  # the values on this line and the next find h=(0, 2.0003, 4.0006, 0) in 21 seconds, 28754 iterations
#    eps_NRconvergence = 1e-15
#    MAXITERS = 64000  # the values on this line and the next find h=(0, 2.0002, 4.0003, 0) in 36 seconds, 51213 iterations
#    eps_NRconvergence = 1e-16
#    MAXITERS = 128000  # the values on this line and the next find h=(0, 2.0001, 4.0002, 0) in 65 seconds, 91146 iterations
#    eps_NRconvergence = 1e-17
#    MAXITERS = 256000  # the values on this line and the next find h=(0, 2.0000, 4.0001, 0) in 117 seconds, 162155 iterations
#    eps_NRconvergence = 1e-18

#    MAXITERS = 100  # James Hilton's number (large enough to yield convergence with the original monotonic flow-function)
    MAXITERS = 400  # larger number, to force Newton-Raphson to reach the same solutions (in the sense of their having balanced flow) that the global minimiser finds, when the flow-function is the (non-monotonic) concave quadratic
#    eps_NRconvergence = 1e-3  # James Hilton's value (small enough to yield convergence with the original monotonic flow-function)
    eps_NRconvergence = 7e-11  # smaller value, to force Newton-Raphson to reach the same solutions (in the sense of their having balanced flow) that the global minimiser finds, when the flow-function is the (non-monotonic) concave quadratic

    h = np.zeros(N)
    h[exitnodes] = terminatorPressure
#    hInjectionNodes = terminatorPressure*np.ones( len(injectionnodes) )
#    h[injectionnodes] = 0.0  # shouldn't need as h is set to np.zeros(N) above
###    h[1] = 2.5  # specify a non-zero initial guess for h
###    h[2] = 0.0
#    h[1], h[2], h[3] = (0, 1.0, 4.0)
#    h[1], h[2] = (1.0, 4.0)
    for ind, jn in enumerate(injectionnodes):
      h[jn] = hinitInjectionNodes[ind]
    print("Initial guess is h = %s" % h)
    dh_dot = dh0_dot = 1.0
#    b = np.zeros(N)
    b = terminatorPressure*np.ones(N)
    d = np.zeros(N)
    iteration = 0
    while True:
      print("iteration %d:" % iteration)
      A = np.zeros( shape=(N,N) )
#      b = np.array(b0)  # make deep copy of b0, so modifying b won't affect b0
      for k in range(N):  # node N-1 won't have any neighbours not already considered
        if k in exitnodes:
#          b[k] = 0.0
          b[k] = terminatorPressure
        else:
          b[k] = b0[indexamonginjectionorexitnodes[k]]
        d[k] = 0.0  # reset diagonal
      debug("h is %s" % h)
#      if lk != None:
#        print("lK %.8f, lk %.5f" % (lK, lk))
#      else:
#        print("lK %.8f, lk %s" % (lK, lk))
      for i in range(N):
        if i in exitnodes:
          hi = hExitNodes[indexamonginjectionorexitnodes[i]]
        else:
#          hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
          hi = h[i]
        for j in [k for k in range(i+1,N) if k in neighboursofj[i]]:
          if j in exitnodes:
            hj = hExitNodes[indexamonginjectionorexitnodes[j]]
          else:
#            hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
            hj = h[j]
          hval = abs(hi-hj)
          print("i %d, j %d, hi %.3f, hj %.3f, hval %.3f" % (i, j, hi, hj, hval))
          qval = dval = 0.0
##          qval = flowijfunc[flowfunc](i, j, hi, hj)
##          dval = derivijfunc[flowfunc](hval)
##          qval = flowfunction(i, j, hi, hj, hvalmax, lK)
##          qval = flowfunction(i, j, hi, hj)
#          qval = flowfunction(i, j, hi, hj, argsflowfunc)
          qval = flowfunction(i, j, hi-hj, argsflowfunc)
#          dval = derivfunction(hval)
#          dval = derivfunction(hval, i, j)
          dval = derivfunction(hval, i, j, argsflowfunc)
          print("qval is %.5f, dval is %.5f" % (qval, dval))
          if i in injectionnodes and j in injectionnodes:  # we set A[i,j] and A[j,i] only for an edge {i,j} with non-leaf endpoint-nodes. On the other hand, for an edge {i,j} adjacent to a exit-node we don't set A[i,j] or A[j,i], the (i,j)th and (j,i)th elements of J_F(h); this is because, if w.l.o.g. i is the exit-node, then F_i = h_i and so A[i,j] = \partial F_i / \partial h_j = 0, and in the system of equations -A*dh = F(h), A[j,i] will be multiplied by dh(i), which is zero (as h(i) is fixed at 'terminatorPressure'), meaning that A[j,i] might as well be zero.
            A[i, j] = -dval
            A[j, i] = -dval
            debug("A[%d,%d] and A[%d,%d] set to -dval of %.5f" % (i, j, j, i, -dval))
          if i in injectionnodes:
            d[i] += dval
            debug("d[%d] += dval of %.5f" % (i, dval))
            b[i] -= qval
            debug("b[%d] -= qval of %.5f" % (i, qval))
          if j in injectionnodes:
            d[j] += dval
            debug("d[%d] += dval of %.5f" % (j, dval))
            b[j] += qval
            debug("b[%d] += qval of %.5f" % (j, qval))
#          if i in injectionnodes:
#            d[i] += dval
#            print("d[%d] += dval of %.5f" % (i, dval))
#            b[i] -= qval
#            print("b[%d] -= qval of %.5f" % (i, qval))
#            A[i, j] = -dval
#            print("A[%d,%d] set to -dval of %.5f" % (i, j, -dval))
#          if j in injectionnodes:
#            d[j] += dval
#            print("d[%d] += dval of %.5f" % (j, dval))
#            b[j] += qval
#            print("b[%d] += qval of %.5f" % (j, qval))
#            A[j, i] = -dval
#            print("A[%d,%d] set to -dval of %.5f" % (j, i, -dval))
      for i in range(N):
        if d[i] == 0.0:
          A[i, i] = 1.0
        else:
          A[i, i] = d[i]
      print("A is %s" % A)
      print("b is %s" % b)
      # Now b = F(h), and A = -J_F: we solve -J_F*dh = F(h).
      try:
        dh = scipy.linalg.solve(A, b)
      except np.linalg.LinAlgError as e:
        print("=========== ALERT: THE SOLUTION FAILS at iteration %d: np.linalg.LinAlgError exception raised. ========================" % iteration)
        return (e, None, None, False, iteration)  # e, hInjectionNodes, qij, flowbalanced
      print("dh is %s" % dh)
      h += dh  # if h falls outside the bounds if any, e.g. |h_i| <= 4 for flowfunc = 'concaveQuadratic', should h be adjusted to fall within the bounds? If so, should h be placed at the boundary or should it be "reflected" back within the bounds? (Probably the former, I imagine, as I don't know a way to justify "reflection," which might anyway cause the adjusted h to fall outside the bounds in the opposite direction.) But I believe that h will fall outside the bounds only if there is no root within the bounds for it to find, so adjusting h would merely delay the information that no root has been found.
      try:
#        dh_dot = scipy.linalg.norm( dh )  # incorrect?
        dh_dot = scipy.linalg.norm( dh )**2
      except ValueError as e:  # e.g., dh contains 'inf' or 'NaN'
        print("=========== ALERT: THE SOLUTION FAILS at iteration %d: ValueError exception raised, perhaps by 'inf' or 'NaN'. ========================" % iteration)
        return (e, None, None, False, iteration)  # e, hInjectionNodes, qij, flowbalanced
      except OverflowError as e:  # e.g., (34, 'Numerical result out of range')
        print("=========== ALERT: THE SOLUTION FAILS at iteration %d: OverflowError exception raised, perhaps because numerical result is out of range. ========================" % iteration)
        return (e, None, None, False, iteration)  # e, hInjectionNodes, qij, flowbalanced
      print("dh_dot is %.5f" % dh_dot)
      print("h is %s" % h)
      if iteration == 0:
        print("iteration is %d, so setting dh0_dot to dh_dot=%.5f" % (iteration, dh_dot))
        dh0_dot = dh_dot
      iteration += 1
#      if dh_dot/dh0_dot <= 1e-3 or iteration >= MAXITERS:
      if dh_dot/dh0_dot <= eps_NRconvergence or iteration >= MAXITERS:
        break
    print("Solution for all nodes is h =", h)
    hInjectionNodes = h[injectionnodes]
#    print("fnormsquared(hInjectionNodes) is %.6e" % fnormsquared(hInjectionNodes))
##    print("fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0) is %.6e" % fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0))
#    print("fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) is %.6e" % fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
    print("fnorm(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) is %.6e" % fnorm(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
    normsratio = dh_dot/dh0_dot
#    if dh_dot/dh0_dot > eps_NRconvergence and iteration >= MAXITERS:
    if normsratio > eps_NRconvergence and iteration >= MAXITERS:
      print("=========== ALERT: Newton-Raphson method has not converged: has reached maximum number (%d) of iterations without ||dh||^2 (%.12f) falling to %.4e * ||dh0||^2 (||dh0||^2 = %.5f); ||dh||^2/||dh0||^2 = %.6e" % (MAXITERS, dh_dot, eps_NRconvergence, dh0_dot, normsratio))
#    if dh_dot/dh0_dot <= eps_NRconvergence:
    if normsratio <= eps_NRconvergence:
      print("Newton-Raphson method has converged, in that ||dh||^2 (%.12f) has fallen to no more than %.4e * ||dh0||^2 (||dh0||^2 = %.5f); ||dh||^2/||dh0||^2 = %.6e" % (dh_dot, eps_NRconvergence, dh0_dot, normsratio))


  print("Solution for injection-nodes is hInjectionNodes =", hInjectionNodes)
#  print("Inflows at nodes for this solution is %s" % f(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0))
  print("Inflows at nodes for this solution is %s" % f(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
  if NEWTONRAPHSON:
    print("Flow network took %d Newton-Raphson iterations." % iteration)
#  print("Function f was called %d times." % numcallsoffuncf)
##  print("The solution-finding method took %.5f seconds." % (endtime-starttime))

##   Calculate Qmax (maximum flow) and totaloutflow
#   Calculate Qabsmax (maximum absolute flow) and totaloutflow
  Qabsmax = 0.0
  qij = {}  # store flows on all edges
#  congestionij = {}  # store congestion-values on all edges
  for i in range(N-1):
    if i in exitnodes:
      hi = hExitNodes[indexamonginjectionorexitnodes[i]]
    else:
      hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
    for j in [k for k in range(i+1,N) if k in neighboursofj[i]]:
      if j in exitnodes:
        hj = hExitNodes[indexamonginjectionorexitnodes[j]]
      else:
        hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
#      hval = abs(hi-hj)
      print("i %d, j %d, hi %.3f, hj %.3f, hval %.3f" % (i, j, hi, hj, abs(hi-hj)))
#      qij[i, j] = flowijfunc[flowfunc](i, j, hi, hj, argsflowfunc)
      qij[i, j] = flowijfunc[flowfunc](i, j, hi-hj, argsflowfunc)
      print("i %d, j %d, hi %.5f, hj %.5f, qij[%d,%d] is %.5f" % (i, j, hi, hj, i, j, qij[i,j]))
      Qabsmax = max(Qabsmax, abs(qij[i,j]))
  print("Qabsmax is", Qabsmax)

  outflowsByNodeID = {}  # outflow at an exit-node is defined to be the negative of its total inflow
  totaloutflow = 0.0
  for t in exitnodes:
    outflowsByNodeID[t] = 0.0
#    totaloutflowj = 0.0
#    ht = h[indexamonginjectionorexitnodes[t]]
    ht = hExitNodes[indexamonginjectionorexitnodes[t]]
#    print("Node %d has neighbours %s" % (t, neighboursofj))
#    assert len(neighboursofj[t]) == 1  # no longer valid, as can now choose a non-leaf node to be a exit- (i.e. exit-) node
    for i in neighboursofj[t]:
#      assert i in injectionnodes  # might be invalid now that any node can be an exit-node: can have adjacent exit-nodes
##      if i in exitnodes:
###        hi = hExitNodes[indexamonginjectionorexitnodes[i]]
##       pass
##      else:
##        hi = h[indexamonginjectionorexitnodes[i]]
#      hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
      if i in injectionnodes:
        hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
      else:
        hi = hExitNodes[indexamonginjectionorexitnodes[i]]
      hval = abs(hi-ht)
      if (i,t) in qij:
        qit = qij[i, t]
      else:
        qit = -qij[t, i]
      print("i %d, t %d, hi %.3f, ht %.3f, hval %.3f, qit %.5f" % (i, t, hi, ht, hval, qit))
#      print("qit is %.5f" % qit)
      outflowsByNodeID[t] -= qit  # outflow is the negative of inflow
      totaloutflow -= qit  # this is really the negative of total outflow
  print("b0 (inflow at injection-nodes) is", b0)
#  if flowfunc == 'concaveQuadratic':
#    print("hvalmax is %s" % hvalmax)
  totalinflow = sum(b0)
  print("totalinflow is", totalinflow, " totaloutflow is", totaloutflow)
  print("totaloutflow + totalinflow is", totaloutflow + totalinflow)
  print("(totaloutflow + totalinflow)/Qabsmax is", (totaloutflow + totalinflow)/Qabsmax)
  EPS_FLOWBALANCE = 1e-3  # James Hilton's value
#  if abs(totaloutflow+totalinflow)/Qmax > 1.0e-3:
  if abs(totaloutflow+totalinflow)/Qabsmax > EPS_FLOWBALANCE:
    if NEWTONRAPHSON:
      print("=========== ALERT: THE SOLUTION FAILS at iteration %d: flow does not balance over network. ========================" % iteration)
    else:
      print("=========== ALERT: THE SOLUTION FAILS: flow does not balance over network. ========================")
#    flowbalanced = 'flow does not balance over network'
    flowbalanced = False
  else:
    if NEWTONRAPHSON:
      print("SUCCESS: METHOD HAS YIELDED A GOOD SOLUTION at iteration %d - flow balances over the network. =======" % iteration)
    else:
      print("SUCCESS: METHOD HAS YIELDED A GOOD SOLUTION - flow balances over the network. =======")
#    flowbalanced = 'flow balances over network'
    flowbalanced = True
  exception = None  # no numpy.linalg.LinAlgError exception was raised while trying to solve the linear system
  if not NEWTONRAPHSON:
    iteration = None
#  print("hInjectionNodes is", hInjectionNodes)
#  print("hExitNodes is", hExitNodes)
  (congestionij, totalnetworkcongestion) = congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
  return (exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, iteration, outflowsByNodeID)




def generateAllIndices( n, numpointsperdimension ):  # return a list of n-tuples, each tuple containing n integers between zero and (numpointsperdimension-1) inclusive, with each possible such tuple included in the list exactly once
  if n == 1:
#    return [(val,) for val in range(maxindex+1)]
    return [(val,) for val in range(numpointsperdimension)]
  elif n > 1:
#    return [ (val,) + t for val in range(maxindex+1) for t in generateAllIndices(n-1, maxindex) ]
    return [ (val,) + t for val in range(numpointsperdimension) for t in generateAllIndices(n-1, numpointsperdimension) ]



def FloydWarshall( edges, vertices, edgeweights ):
  """
  A dynamic-programming algorithm: find shortest path between every pair of vertices, given the distance between each pair of adjacent vertices as stored in 'edgeweights'. In the returned arrays, pathdistance[i,j] is the shortest-path distance from i to j, and nextvertex[i,j] is the next vertex after i in the shortest path from i to j.
  """
  N = len(vertices)
#  pathdistance = 9e99 * np.ones( shape=(N,N) )  # set each path-distance to "infinity"
  pathdistance = np.full( shape=(N,N), fill_value=9e99  )  # set each path-distance to "infinity"
  nextvertex = np.full( shape=(N,N), fill_value=None )
  for e in edges:
    eAsList = list(e)
#    print("eAsList is", eAsList)
    (a, b) = eAsList
#    pathdistance[a, b] = edgeweights[e]
    pathdistance[a, b] = edgeweights[a, b]
    nextvertex[a, b] = b
  for v in vertices:
    pathdistance[v, v] = 0
    nextvertex[v, v] = v
#  print("FloydWarshall(): pathdistance is", pathdistance)
#  print("FloydWarshall(): nextvertex is", nextvertex)
  for k in vertices:
    for i in vertices:
      for j in vertices:
        if pathdistance[i, j] > pathdistance[i, k] + pathdistance[k, j]:
          pathdistance[i, j] = pathdistance[i, k] + pathdistance[k, j]
          nextvertex[i, j] = nextvertex[i, k]
#          print("k %d: pathdistance[%d, %d] is %s" % (k, i, j, pathdistance[i,j]))
#          print("nextvertex[%d, %d] is %s" % (i, j, nextvertex[i,j]))
  return (pathdistance, nextvertex)


def constructShortestPath(i, j, nextvertex, edgeweight):
  if nextvertex[i, j] == None:
    return []
  path = [i]
#  while i != j:
#    i = nextvertex[i, j]
#    path.append( i )
  pathlength = 0.0
  a = i
  while a != j:
    b = nextvertex[a, j]
    pathlength += edgeweight[a, b]
    path.append( b )
    a = b
  return (path, pathlength)


def adjustSpeedandDensityToLink( flowfunction, p, j, densityenteringparent, argsflowfunc, flowenteringparent, speedenteringparent ):
  EPS = 1e-3
#  EPS = 1e-0  # might be all the precision that's required
  densityinlink = densityenteringparent
  previousspeedinlink = 1e9
  speedinlink = speedenteringparent
  numiterations = 0
#  if flowfunc in ('concaveQuadratic', 'xExpnegx', 'triangular'):
#    humpdensity = argsflowfunc['tau'][frozenset({p,j})]
#  elif flowfunc == 'triangular':
#    humpdensity = numlanes*argsflowfunc['rhoC'][frozenset({p,j})]
  humpdensity = argsflowfunc['tau'][frozenset({p,j})]
#  if densityinlink >= humpdensity:
  if densityinlink > humpdensity:
    print(f"WARNING: density in link {{{p},{j}}} of {densityinlink:.6f} is over the hump {humpdensity:.6f} by a gap of {(densityinlink-humpdensity):.6e}.")
  # Iterate until speedinlink converges:
  CALCULATESPEEDBEFOREDENSITY = True
#  CALCULATESPEEDBEFOREDENSITY = False
  print("CALCULATESPEEDBEFOREDENSITY is", CALCULATESPEEDBEFOREDENSITY)
  while abs(speedinlink - previousspeedinlink) > EPS:
    previousspeedinlink = speedinlink
    # First, as speed is a function of density, the speed adjusts in accordance with the density and the conditions (free speed, number of lanes) on the link
    theoreticalflowinlink = flowfunction(p, j, densityinlink, argsflowfunc)  # flow that, for the given density, would result from the conditions (free speed, number of lanes) on this link
    print(f"iteration {numiterations}: given density of {densityinlink:.3f}, theoreticalflowinlink is {theoreticalflowinlink:.3f}")
    if CALCULATESPEEDBEFOREDENSITY:
#      if abs(densityinlink) >= EPS:
      if abs(densityinlink) >= EPS or abs(theoreticalflowinlink) >= EPS:
        speedinlink = theoreticalflowinlink / densityinlink  # v = q/rho: the speed adjusts in accordance with the density and the conditions (free speed, number of lanes) on the link
#        print(f"adjustSpeedandDensityToLink(): speedinlink adjusted to {speedinlink:.3E}")
        print(f"speedinlink adjusted to {theoreticalflowinlink:.3f}/{densityinlink:.3f} = {speedinlink:.3f}")
      else:  # both flow and density in link are close to zero: speedinlink retains its existing value
        pass
#      assert speedinlink >= 0.0
      assert speedinlink >= -EPS, f"speedinlink is {speedinlink:.3f}, numiterations is {numiterations}"
      # Second, the density adjusts in accordance with the original entering flow and the new speed on the link (density = flow/speed)
      if abs(speedinlink) >= 1e-18:
        densityinlink = flowenteringparent / speedinlink  # the density adjusts in accordance with the new adjusted speed and the constraint that flow is conserved
        print(f"densityinlink adjusted to {flowenteringparent:.3f}/{speedinlink:.3f} = {densityinlink:.3f}")
      else:
        densityinlink = 9e99  # TODO: is it valid to set this to a high value, or should it be set to infinity/inf/NAN?
#        densityinlink = 1e18
#        densityinlink = 9.99999999999e17
        print(f"densityinlink adjusted to {densityinlink:.3f}")
      print(f"densityinlink is {densityinlink:.3f}")
    else:  # CALCULATESPEEDBEFOREDENSITY == False
      if abs(densityinlink) >= EPS or abs(theoreticalflowinlink) >= EPS:
        densityinlink = theoreticalflowinlink / speedinlink  # rho = q/v: the density adjusts in accordance with the speed and the conditions (free speed, number of lanes) on the link
        print(f"densityinlink adjusted to {theoreticalflowinlink:.3f}/{speedinlink:.3f} = {densityinlink:.3f}")
      else:  # both flow and density in link are close to zero: densityinlink retains its existing value
        pass
      assert densityinlink >= -EPS, f"densityinlink is {densityinlink:.3f}, numiterations is {numiterations}"
      # Second, the speed adjusts in accordance with the original entering flow and the new density on the link (speed = flow/density)
      if abs(densityinlink) >= 1e-18:
        speedinlink = flowenteringparent / densityinlink  # the speed adjusts in accordance with the new adjusted density and the constraint that flow is conserved
        print(f"speedinlink adjusted to {flowenteringparent:.3f}/{densityinlink:.3f} = {speedinlink:.3f}")
      else:
        speedinlink = 9e99  # TODO: is it valid to set this to a high value, or should it be set to infinity/inf/NAN?
#        speedinlink = 1e18
        print(f"speedinlink adjusted to {speedinlink:.3f}")
      print(f"speedinlink is {speedinlink:.3f}")
    numiterations += 1
  print(f"abs(theoreticalflowinlink - flowenteringparent) is {abs(theoreticalflowinlink - flowenteringparent):.3E}")
  EPS_FLOWINLINK = 2e-2
  assert abs(theoreticalflowinlink - flowenteringparent) <= EPS_FLOWINLINK
  print(f"adjustSpeedandDensityToLink(): numiterations is {numiterations}")
  print(f"adjustSpeedandDensityToLink(): {densityinlink:.6f}*{speedinlink:.6f} = {densityinlink*speedinlink:.6f}")
  return (speedinlink, densityinlink, numiterations)


def densityonfreeflowbranchoftriangularfundamentaldiagram( freeflow, edge, argsflowfunc ):
  freespeed = argsflowfunc['lK'][edge]
  return freeflow / freespeed  # density over entire link-width; from eq. (8.13), p. 94 of Treiber and Kesting, substituting Q=q_tot/I and rho_{free} = rho_tot/I


def densityoncongestedbranchoftriangularfundamentaldiagram( congestedflow, edge, argsflowfunc ):
  T = argsflowfunc['T'][edge]
  numlanes = argsflowfunc['numlanes'][edge]
  maximumdensityperlane = argsflowfunc['maximumdensityperlane']
  return numlanes * (1 - T*congestedflow/numlanes) * maximumdensityperlane


def propagationSpeedOfShockFront( upstreamflow, upstreamdensity, downstreamflow, downstreamdensity ):
  EPS = 1e-3
  assert abs(downstreamdensity - upstreamdensity) > EPS
  return (downstreamflow - upstreamflow) / (downstreamdensity - upstreamdensity)


def propagateSEM3congestionUpstreamToParents( i, totalflowleavingnodeiDownstreamOfCongestion, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes ):
  print(f"totalflowleavingnodeiDownstreamOfCongestion at node {i} is {totalflowleavingnodeiDownstreamOfCongestion}")
  initialtotalflowEnteringNodeiFromItsParents = sum([qij[p, i] for p in parentsOnShortestPathsFromInjectionNodes[i]])
  print(f"Initial total flow entering node {i} from its parent-nodes is {initialtotalflowEnteringNodeiFromItsParents:.3f}")
  totalflowEnteringNodei = 0.0
  if i in injectionnodes:
    injectionflowati = b0[indexamonginjectionorexitnodes[i]]
    if injectionflowati > 0.0:
      print(f"Adding to totalflowEnteringNodei the injectionflow at node {i} of {injectionflowati:.5f}")
    totalflowEnteringNodei += injectionflowati
  for p in parentsOnShortestPathsFromInjectionNodes[i]:
##    qij[p, i] *= totalflowleavingnodeiDownstreamOfCongestion / initialtotalflowEnteringNodei  # divide congested flow between the links entering node i: the flow in each decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
#    qij[p, i] *= (totalflowleavingnodeiDownstreamOfCongestion - injectionflowati) / initialtotalflowEnteringNodeiFromItsParents  # divide congested flow between the links entering node i: the flow in each decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
#    congestedflowpi = qij[p, i]
    initialflowpi = qij[p, i]
    print(f"initial flow[{p}, {i}] is {initialflowpi}")
    if initialflowpi == 0.0:  # no need to consider a parent-link that has no incoming flow, as congestion can't affect it  # FIXME: should check for closeness to zero of flow's absolute value
      print(f"initialflowpi of {initialflowpi} is zero, which cannot be reduced by congestion: there is no real congestion in link ({p}, {i}).")
      continue
#    congestedflowpi = qij[p, i] * (totalflowleavingnodeiDownstreamOfCongestion - injectionflowati) / initialtotalflowEnteringNodeiFromItsParents  # divide congested flow between the links entering node i: the flow in each decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
    congestedflowpi = qij[p, i] * totalflowleavingnodeiDownstreamOfCongestion / (initialtotalflowEnteringNodeiFromItsParents + injectionflowati)  # divide congested flow between the links entering node i: the flow in each link, and also by implication the injection-flow, decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
    print(f"congested flow[{p}, {i}] is {congestedflowpi}")
    print(f"As congestion-zone propagates upstream, flow in link [{p}, {i}] decreases from {initialflowpi} to {congestedflowpi}")
    qij[p, i] = congestedflowpi
    totalflowEnteringNodei += congestedflowpi
    edge = frozenset({p, i})
#    densityij[frozenset({p, i})] = densityoncongestedbranchoftriangularfundamentaldiagram( qij[p, i], edge, argsflowfunc )
    congesteddensitypi = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowpi, edge, argsflowfunc )
    densityij[edge] = congesteddensitypi
    print(f"densityij IS NOW {densityij}")
    upstreamfreeflowpi = initialflowpi
    upstreamfreedensitypi = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamfreeflowpi, edge, argsflowfunc )
    print(f"Calculating propagation-speed: upstreamfreeflowpi is {upstreamfreeflowpi}, upstreamfreedensitypi is {upstreamfreedensitypi}, congestedflowpi is {congestedflowpi}, congesteddensitypi is {congesteddensitypi}")
    shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamfreeflowpi, upstreamfreedensitypi, congestedflowpi, congesteddensitypi )
    if shockFrontPropagationSpeed < 0.0:
      upordownstreamText = 'travels upstream'
    elif shockFrontPropagationSpeed > 0.0:
      upordownstreamText = 'travels downstream'
    else:
      upordownstreamText = 'is stationary'
    print(f"shockFrontPropagationSpeed is {shockFrontPropagationSpeed}: the front {upordownstreamText}.")
#    assert shockFrontPropagationSpeed < 0.0
    if shockFrontPropagationSpeed == 0.0:  # can be zero if, for example, the link's free flow upstream of the shock-front is zero, as then the "congested" flow is also zero
      print(f"upstreamfreeflowpi {upstreamfreeflowpi} = congestedflowpi {congestedflowpi}: there is no real congestion.")
#    speedij[frozenset({p, i})] = qij[p, i] / densityij[edge]
    speedij[edge] = congestedflowpi / congesteddensitypi
    print(f"Density in link {edge} is now {densityij[edge]}, speed is now {speedij[edge]}.")
    propagateSEM3congestionUpstreamToParents( p, congestedflowpi, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )

  if totalflowleavingnodeiDownstreamOfCongestion < totalflowEnteringNodei:
#    print(f"ALERT: totalflowEnteringNodei of {totalflowEnteringNodei} not equal to totalflowleavingnodeiDownstreamOfCongestion of {totalflowleavingnodeiDownstreamOfCongestion}.")
    print(f"ALERT: at node {i}, totalflowleavingnodeiDownstreamOfCongestion {totalflowleavingnodeiDownstreamOfCongestion} < totalflowEnteringNodei {totalflowEnteringNodei}: node {i} must be an injection-node, and there is a shock where traffic-flow enters the network there.")
    assert i in injectionnodes
#  assert totalflowEnteringNodei == totalflowleavingnodeiDownstreamOfCongestion
  assert totalflowleavingnodeiDownstreamOfCongestion <= totalflowEnteringNodei


#def propagateSEM4congestionUpstreamToParents( upstreamPropagatingShockFront, a, p, congestedFlowThatHasPropagatedUpstreamToNodea, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes ):
def propagateSEM4congestionUpstreamToParents( upstreamPropagatingShockFront, a, timeThatShockFrontReachedNodea, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, flowstateijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes, propagatingwavefronts, nexttimestepduration ):
  EPS = 1e-3
  congestedFlowThatHasPropagatedUpstreamToNodea = upstreamPropagatingShockFront['flow']
  print(f"congestedFlowThatHasPropagatedUpstreamToNodea to node {a} is {congestedFlowThatHasPropagatedUpstreamToNodea}")
  initialtotalflowEnteringNodeaFromItsParents = sum([qij[p, a] for p in parentsOnShortestPathsFromInjectionNodes[a]])
  print(f"Initial total flow entering node {a} from its parent-nodes is {initialtotalflowEnteringNodeaFromItsParents:.3f}")
  newtotalflowEnteringNodea = 0.0  # will include the injection-flow at a, if any
  if a in injectionnodes:
    injectionflowata = b0[indexamonginjectionorexitnodes[a]]
    if injectionflowata > 0.0:
      print(f"Adding to newtotalflowEnteringNodea the injectionflow at node {a} of {injectionflowata:.5f}")
    newtotalflowEnteringNodea += injectionflowata
  for p in parentsOnShortestPathsFromInjectionNodes[a]:
    initialflowpa = qij[p, a]
    print(f"initial flow[{p}, {a}] is {initialflowpa}")
    if initialflowpa == 0.0:  # no need to consider a parent-link that has no incoming flow, as congestion can't affect it  # FIXME: should check for closeness to zero of flow's absolute value
      print(f"initialflowpa of {initialflowpa} is zero, which cannot be reduced by congestion: there is no real congestion in link ({p}, {a}).")
      continue
    congestedflowpa = initialflowpa * congestedFlowThatHasPropagatedUpstreamToNodea / (initialtotalflowEnteringNodeaFromItsParents + injectionflowata)  # divide congested flow between the links entering node a: the flow in each link, and also by implication the injection-flow, decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
    print(f"congested flow[{p}, {a}] is {congestedflowpa}")
    print(f"As congestion-zone propagates upstream, flow in link [{p}, {a}] decreases from {initialflowpa} to {congestedflowpa}")
    qij[p, a] = congestedflowpa
    newtotalflowEnteringNodea += congestedflowpa
    edge = frozenset({p, a})
    congesteddensitypa = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowpa, edge, argsflowfunc )
    densityij[edge] = congesteddensitypa
    print(f"densityij is now {densityij}")
    upstreamflowpa = initialflowpa
    if flowstateijatnodej[p, a] == 'freeflowstate':
      upstreamdensitypa = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamflowpa, edge, argsflowfunc )
    elif flowstateijatnodej[p, a] == 'congestedflowstate':
      upstreamdensitypa = densityoncongestedbranchoftriangularfundamentaldiagram( upstreamflowpa, edge, argsflowfunc )
    print(f"Calculating propagation-speed: upstreamflowpa is {upstreamflowpa}, upstreamdensitypa is {upstreamdensitypa}, congestedflowpa is {congestedflowpa}, congesteddensitypa is {congesteddensitypa}")
    shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamflowpa, upstreamdensitypa, congestedflowpa, congesteddensitypa )
    if shockFrontPropagationSpeed < 0.0:
      upordownstreamText = 'travels upstream'
    elif shockFrontPropagationSpeed > 0.0:
      upordownstreamText = 'travels downstream'
    else:
      upordownstreamText = 'is stationary'
    print(f"shockFrontPropagationSpeed is {shockFrontPropagationSpeed}: the front {upordownstreamText}.")
    assert shockFrontPropagationSpeed < 0.0
    if shockFrontPropagationSpeed == 0.0:  # can be zero if, for example, the link's free flow upstream of the shock-front is zero, as then the "congested" flow is also zero
      print(f"upstreamflowpa {upstreamflowpa} = congestedflowpa {congestedflowpa}: there is no real congestion.")
#    speedij[frozenset({p, a})] = qij[p, a] / densityij[edge]
    speedij[edge] = congestedflowpa / congesteddensitypa
    flowstateijatnodej[p, a] = 'congestedflowstate'
    print(f"Density in link {edge} is now {densityij[edge]}, speed is now {speedij[edge]}, flowstateijatnodej[{p}, {a}] is {flowstateijatnodej[p, a]}.")
    lengthoflinkpa = argsflowfunc['linklength'][ frozenset({p,a}) ]
    timetillupstreamnodeonlink = -lengthoflinkpa / shockFrontPropagationSpeed
##    propagateSEM3congestionUpstreamToParents( p, congestedflowpa, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )
    propagatingwavefronts.append( {'time': timeThatShockFrontReachedNodea, 'currentlink': frozenset({p,a}), 'sourcenodeonlink': a, 'destinationnodeonlink': p, 'distancereacheddownlink': 0.0, 'speed': -shockFrontPropagationSpeed, 'flow': congestedflowpa, 'density': congesteddensitypa, 'timetilldestinationnodeonlink': lengthoflinkpa/(-shockFrontPropagationSpeed), 'nodesonpathtoanexit': 'upstream-propagating shock-front, so not applicable', 'flowaheadofwavefront': upstreamflowpa, 'densityaheadofwavefront': upstreamdensitypa} )

  upstreamPropagatingShockFront['timetilldestinationnodeonlink'] = 0.0
  upstreamPropagatingShockFront['distancereacheddownlink'] += upstreamPropagatingShockFront['speed'] * nexttimestepduration
  sourceNodeOfShockFront = upstreamPropagatingShockFront['sourcenodeonlink']
  lengthofcurrentlink = argsflowfunc['linklength'][ frozenset({sourceNodeOfShockFront, a}) ]
  assert abs(upstreamPropagatingShockFront['distancereacheddownlink'] - lengthofcurrentlink) <= EPS
  print(f"As it propagates upstream at node {a}, wave-front {upstreamPropagatingShockFront} is removed from propagatingwavefronts.")
  propagatingwavefronts.remove( upstreamPropagatingShockFront )
  print(f"propagatingwavefronts is {propagatingwavefronts}")

  if congestedFlowThatHasPropagatedUpstreamToNodea < newtotalflowEnteringNodea:
#    print(f"ALERT: newtotalflowEnteringNodea of {newtotalflowEnteringNodea} not equal to congestedFlowThatHasPropagatedUpstreamToNodea of {congestedFlowThatHasPropagatedUpstreamToNodea}.")
    print(f"ALERT: at node {a}, congestedFlowThatHasPropagatedUpstreamToNodea {congestedFlowThatHasPropagatedUpstreamToNodea} < newtotalflowEnteringNodea {newtotalflowEnteringNodea}: node {a} must be an injection-node, and there is a shock where traffic-flow enters the network there.")
    assert a in injectionnodes
#  assert newtotalflowEnteringNodea == congestedFlowThatHasPropagatedUpstreamToNodea
  assert congestedFlowThatHasPropagatedUpstreamToNodea <= newtotalflowEnteringNodea


#def flowandspeedEnteringNode( i, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnExit, flowijfunc, flowfunc ):
def flowEnteringNode( i, childofi, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes ):
  # Traffic injected at a node flows towards the next node on its shortest path towards an exit.
#  # At the beginning of the first link the flow and speed are known, which gives the density; for traffic flowing into the link this density gives, according to v = v(rho) for the link, the speed in the link. The flow is conserved between each injection-node and the end of the first link on its shortest path to an exit; if a link is divided for numerical reasons into sub-sections, then the flow is conserved between adjacent sections.
  # At an injection-node the inflow is known; this divides among any adjacent links that are used by shortest paths, and the portion of flow on each such link gives its density, on the free-flow side of its triangular fundamental diagram; the speed is the link's free speed V0 (see Figure 8.9, p. 92 of Treiber and Kesting).
  EPS = 1e-3
  print(f"Entering flowEnteringNode( {i}, ... ):")
  if i in exitnodes:
#    print(f"Node {i} is an exit-node, so has zero injection-flow and -density, and injection-speed is undefined.")
    print(f"Node {i} is an exit-node, so has zero injection-flow.")
#    injectionflow = injectiondensity = 0.0
    injectionflow = 0.0
##    injectionspeed = 'undefined'
##    print(f"At node {i} injectionflow is {injectionflow:.3f}, injectionspeed is {injectionspeed}, and injectiondensity is {injectiondensity:.3f}.")
#    print(f"At node {i} injectionflow is {injectionflow:.3f}.")
  else:
    injectionflow = b0[indexamonginjectionorexitnodes[i]]
    print(f"injectionflow is {injectionflow:.5f}")
#    injectionspeed = v0[indexamonginjectionorexitnodes[i]]
    speedinlinkofinjectedflow = argsflowfunc['lK'][frozenset({i,childofi})]
#    print(f"speedinlinkofinjectedflow is {speedinlinkofinjectedflow:.5f}")
    print(f"speed in link [{i}, {childofi}] of injected flow is {speedinlinkofinjectedflow}")
#    injectiondensity = injectionflow / speedinlinkofinjectedflow
##    print(f"injectiondensity is {injectiondensity:.5f}")
#    print(f"injectiondensity is {injectiondensity}")
#    print(f"At node {i} injectionflow is {injectionflow:.3f} and speedinlinkofinjectedflow is {speedinlinkofinjectedflow:.3f}, so injectiondensity is {injectiondensity:.3f}.")
    print(f"At node {i} injectionflow is {injectionflow:.3f} and speedinlinkofinjectedflow is {speedinlinkofinjectedflow:.3f}.")

  if parentsOnShortestPathsFromInjectionNodes[i] == set():
    print(f"Node {i} has no parent-nodes on shortest paths from injection-nodes.")  # N.B. node i might be an exit-node with no flow entering it
    if i in exitnodes:
      print(f"Node {i} is an exit-node, so should have zero injection-flow and injection-speed")
      assert injectionflow == 0.0
      assert speedinlinkofinjectedflow == 0.0
#      print(f"Node {i} is an exit-node; returning injection-flow {injectionflow:.3f} and speed in link of injected flow {speedinlinkofinjectedflow:.3f}")
      print(f"Node {i} is an exit-node; returning injection-flow {injectionflow:.3f}")
#      return injectionflow, speedinlinkofinjectedflow
      return injectionflow
#    print(f"Node {i} is an injection-node; returning injection-flow {injectionflow:.3f} and speed in link of injected flow {speedinlinkofinjectedflow:.3f}")
    print(f"Node {i} is an injection-node; returning injection-flow {injectionflow:.3f}")
#    return injectionflow, speedinlinkofinjectedflow
    return injectionflow

  # Node i has at least one parent-node on shortest paths from injection-nodes:
  totalflow = injectionflow  # total flow so far: flows from parent-nodes will be added to this
#  totaldensity = injectiondensity  # vehicles per longitudinal metre of road, no matter how many lanes it has laterally
#  speedsinlinks = []  # store speeds in all links flowing into i
  print(f"Node {i}'s parent-nodes on shortest paths from injection-nodes are {parentsOnShortestPathsFromInjectionNodes[i]}")
  for p in parentsOnShortestPathsFromInjectionNodes[i]:
##    flowenteringparent, speedenteringparent = flowandspeedEnteringNode( p, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnExit, flowijfunc, flowfunc )
#    flowenteringparent, speedenteringparent = flowandspeedEnteringNode( p, i, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc )
    flowenteringparent = flowEnteringNode( p, i, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes )
#    densityenteringparent = flowenteringparent / speedenteringparent
#    print(f"Node {p} has closest exit-node {closestExit[p]}, at distance {shortestDistanceToAnExit[p]:.4f} with shortest path {shortestPathToAnExit[p]}")
#    print(f"p is {p}, i is {i}")
    assert i == shortestPathToAnExit[p][1]  # i is the next node in the shortest path from p to an exit-node
    linkcapacity = argsflowfunc['linkcapacity'][frozenset({p,i})]
    if flowenteringparent > linkcapacity:
      # TODO: is it physically realistic to have traffic-flow in excess of the link-capacity use a different route to an exit, perhaps a second-shortest path? Note there might be several different paths all with the shortest length.
#      print(f"WARNING: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {i}) which has a capacity of {linkcapacity:.3f}; high congestion and density, and low speed, will result.")
      print(f"WARNING: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {i}) which has a capacity of {linkcapacity:.3f}: the link will have flow at this capacity, while a congestion-zone will propagate up the upstream link(s); the congestion-zone will have density to the right of the hump in its link's flow-function, and correspondingly low speed.")
      flowleavingparent = linkcapacity
      densityleavingparent = argsflowfunc['tau'][frozenset({p,i})]
      speedleavingparent = flowleavingparent / densityleavingparent
#      for gp in parentsOnShortestPathsFromInjectionNodes[p]:
      propagateSEM3congestionUpstreamToParents( p, flowleavingparent, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )
    else:  # entering flow is not above capacity
      flowleavingparent = flowenteringparent
      speedleavingparent = argsflowfunc['lK'][frozenset({p,i})]
      densityleavingparent = flowleavingparent / speedleavingparent

    print(f"The total outflow at node {p} of {flowleavingparent:.3f} veh/hr starts to flow towards node {i} with density {densityleavingparent:.3f} veh/km and speed v0={speedleavingparent:.3f} km/hr.")
    assert densityleavingparent * speedleavingparent == flowleavingparent
#    flowfunction = flowijfunc[flowfunc]
#    (speedinlink, densityinlink, numiterations) = adjustSpeedandDensityToLink( flowfunction, p, i, densityenteringparent, argsflowfunc, flowenteringparent, speedenteringparent )
#    print(f"After {numiterations} iterations, the speed adjusts to {speedinlink:.3f} km/hr and density to {densityinlink:.3f}.")
    (speedinlink, densityinlink) = (speedleavingparent, densityleavingparent)
    qij[p, i] = flowleavingparent
    densityij[frozenset({p, i})] = densityinlink
    print(f"AFTER DEFINING densityij[frozenset({{{p}, {i}}})], densityij is {densityij}")
    speedij[frozenset({p, i})] = speedinlink
#    if flowleavingparent > linkcapacity:
##      assert abs(speedinlink) < 1e-30  # true in tests so far
#      assert abs(speedinlink) < 1e-1  # all the precision that's required
    totalflow += flowleavingparent
#    totaldensity += densityinlink
#    speedsinlinks.append( speedinlink )
#    print(f"Adding {flowleavingparent:.3f} to totalflow and {densityinlink:.3f} to totaldensity that enter node {i}.")
    print(f"Adding {flowleavingparent:.3f} to totalflow that enters node {i}.")
#  speed = totalflow / totaldensity
#  assert abs(totaldensity) >= EPS  # TODO: is this always true? I put it here to catch any unexpected case in which it isn't
##  if abs(totaldensity) >= EPS:
##    speed = totalflow / totaldensity  # v = q/rho
##  elif abs(totalflow) < EPS:
##    speed = np.mean( speedsinlinks )  # both total flow and total density at node i are close to zero, so we set speed to the mean value of the speeds in inbound links (FIXME: a bit of a hack, although might well be fine in practice)
##  else:  # flow in link is far from zero, density in link is close to zero
##    speed = totalflow / totaldensity  # v = q/rho
#  print(f"At node {i}, totalflow is {totalflow:.3f} and totaldensity is {totaldensity:.3f}, so speed is {speed:.3f}.")
#  print(f"At node {i}, totalflow is {totalflow:.3f} and totaldensity is {totaldensity:.3f}.")
  print(f"At node {i}, totalflow is {totalflow:.3f}.")
#  return totalflow, speed
  return totalflow


def evolutionofwaveandshockfronts( inflowsandflowperiodsByNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, flowstateijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes, outflowsByNodeID ):
  EPS = 1e-6
  INITIALTIMESTAMP = 0.0
#  propagatingwavefronts = {}
#  for edge in argsflowfunc['lK']:
#    propagatingwavefronts[edge] = []
  propagatingwavefronts = []
  # Create at time t=0 the initial wave-fronts at injection-nodes:
  for a in inflowsandflowperiodsByNodeID:
    print(f"For node {a}, inflowsandflowperiodsByNodeID[a] is {inflowsandflowperiodsByNodeID[a]}")
    shortestpathtoanexit_froma = shortestPathToAnExit[a]
    b = shortestpathtoanexit_froma[1]  # b is the next node in the shortest path from a to an exit-node
#    initialflowab = qij[a, b]
    initialflowab = 0.0
    print(f"initial flow[{a}, {b}] is {initialflowab} veh/hr")
    assert initialflowab == 0.0  # FIXME: should check for closeness to zero of flow's absolute value
    print(f"In-flow at node {a} causes a downstream-propagating wave-front to form, with the in-flow - perhaps restricted by capacity of link ({a},{b}) - freely flowing upstream of the wave-front, and zero flow downstream of it.")
    downstreamflowab = initialflowab
    print(f"Downstream flow[{a}, {b}] is {downstreamflowab} veh/hr")
    linkcapacity = argsflowfunc['linkcapacity'][frozenset({a,b})]
    print(f"linkcapacity is {linkcapacity} veh/hr")
    if inflowsandflowperiodsByNodeID[a]['inflow'] <= linkcapacity:
      upstreamfreeflowab = inflowsandflowperiodsByNodeID[a]['inflow']
    else:
      upstreamfreeflowab = linkcapacity
      # If inflow exceeds capacity, report an upstream-propagating shock-wave created at the injection-node that affects only the hypothetical injectinon-"link" by which injection-traffic might be supposed to enter the network:
      print(f"ALERT: at injection-node node {a}, linkcapacity {linkcapacity} < injection-flow {inflowsandflowperiodsByNodeID[a]['inflow']}: there is a shock where traffic-flow enters the network there.")
    print(f"Upstream free flow[{a}, {b}] is {upstreamfreeflowab} veh/hr")
    edge = frozenset({a, b})
    upstreamfreedensityab = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamfreeflowab, edge, argsflowfunc )
#    qij[a, b] = downstreamflowab
    downstreamdensityab = densityonfreeflowbranchoftriangularfundamentaldiagram( downstreamflowab, edge, argsflowfunc )
    print(f"downstreamdensityab is {downstreamdensityab} veh/km")
#    if abs(downstreamflowab) <= EPS:
##      assert (a, b) not in qij
##      assert (b, a) not in qij
#      qij[a, b] = downstreamflowab
#      print(f"qij[{a}, {b}] is now {qij[a, b]}")
##    densityij[edge] = downstreamdensityab
##    print(f"densityij IS NOW {densityij} veh/km")
    print(f"As wave-front propagates downstream, flow in link [{a}, {b}] will increase from {initialflowab} to {upstreamfreeflowab} veh/hr and density increase from {downstreamdensityab} to {upstreamfreedensityab} veh/km.")
    print(f"Calculating propagation-speed: upstreamfreeflowab is {upstreamfreeflowab}, upstreamfreedensityab is {upstreamfreedensityab}, downstreamflowab is {downstreamflowab}, downstreamdensityab is {downstreamdensityab}")
    wavefrontPropagationSpeed = propagationSpeedOfShockFront( upstreamfreeflowab, upstreamfreedensityab, downstreamflowab, downstreamdensityab )
    if wavefrontPropagationSpeed < 0.0:
      upordownstreamText = 'will travel upstream'
    elif wavefrontPropagationSpeed > 0.0:
      upordownstreamText = 'will travel downstream'
    else:
      upordownstreamText = 'is stationary'
    print(f"wavefrontPropagationSpeed is {wavefrontPropagationSpeed} km/hr: the front {upordownstreamText}")
    assert wavefrontPropagationSpeed >= 0.0
    if wavefrontPropagationSpeed == 0.0:  # can be zero if, for example, the inflow from node a into link (a,b) (upstream of the "wave-front") is zero, because the downstream flow is of course also zero
      print(f"upstreamfreeflowab {upstreamfreeflowab} = downstreamflowab {downstreamflowab} (veh/hr): there is no real propagating wave-front.")
    assert wavefrontPropagationSpeed == argsflowfunc['lK'][ frozenset({a,b}) ]
#    speedij[edge] = upstreamfreeflowab / upstreamfreedensityab
    upstreamfreespeedinedgeab = upstreamfreeflowab / upstreamfreedensityab
    assert upstreamfreespeedinedgeab == argsflowfunc['lK'][frozenset({a,b})]
##    print(f"Density in link {edge} is now {densityij[edge]} veh/km, speed is now {speedij[edge]} km/hr.")
#    print(f"Density in link {edge} is now {upstreamfreedensityab} veh/km, speed is now {upstreamfreespeedinedgeab} km/hr.")
    print(f"Traffic-speed in link {edge} will increase from zero to {upstreamfreespeedinedgeab} km/hr.")
    lengthofcurrentlink = argsflowfunc['linklength'][ frozenset({a,b}) ]
    timetilldestinationnodeonlink = lengthofcurrentlink / wavefrontPropagationSpeed
#    if abs(wavefrontPropagationSpeed) > 0.0:  # record only the wave-fronts with non-zero speed
##    propagatingwavefronts[frozenset({a,b})].append( {'time': INITIALTIMESTAMP, 'speed': wavefrontPropagationSpeed, 'nodesonpathtoanexit': shortestpathtoanexit_froma} )
#    propagatingwavefronts.append( {'time': INITIALTIMESTAMP, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': a, 'destinationnodeonlink': b, 'distancereacheddownlink': 0.0, 'speed': wavefrontPropagationSpeed, 'flow': upstreamfreeflowab, 'density': upstreamfreedensityab, 'timetilldestinationnodeonlink': timetilldestinationnodeonlink, 'nodesonpathtoanexit': shortestpathtoanexit_froma} )
    propagatingwavefronts.append( {'time': INITIALTIMESTAMP, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': a, 'destinationnodeonlink': b, 'distancereacheddownlink': 0.0, 'speed': wavefrontPropagationSpeed, 'flow': upstreamfreeflowab, 'density': upstreamfreedensityab, 'timetilldestinationnodeonlink': timetilldestinationnodeonlink, 'nodesonpathtoanexit': shortestpathtoanexit_froma, 'flowaheadofwavefront': 0.0, 'densityaheadofwavefront': 0.0} )
  print(f"propagatingwavefronts is {propagatingwavefronts}")

  simulationclocktime = INITIALTIMESTAMP
  print(f"==============simulationclocktime is {simulationclocktime} hours ({simulationclocktime*3600} seconds).")
  cumulNumVehiclesToHaveDepartedLinkDuringSimulation = {}
  for (i, j) in qij:
    cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j] = 0
  while simulationclocktime < simulDuration and propagatingwavefronts:
    print(f"propagatingwavefronts is {propagatingwavefronts}")
#    nexttimestepduration = min([wv['timetilldestinationnodeonlink'] for wv in propagatingwavefronts])  # the next time-step's duration, during which all wave-fronts will propagate; at the end of the next time-step one (or more) of the wave-fronts will intersect with the next node on its shortest path to an exit
    minTimeForAWavefrontToReachNextNodeOnPath = min([wv['timetilldestinationnodeonlink'] for wv in propagatingwavefronts])  # initial candidate for the next time-step's duration, during which all wave-fronts will propagate; at the end of the minimum time just calculated, at least one of the wave-fronts will intersect with the next node on its shortest path to an exit
    wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes = [wv for wv in propagatingwavefronts if wv['timetilldestinationnodeonlink'] == minTimeForAWavefrontToReachNextNodeOnPath]  # every wave-front that will reach, at the end of the minimum time just calculated, the next node on its shortest path to an exit
    assert len(wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes) == 1  # TODO: modify this assertion if it turns out to be not always true, which it probably will eventually
    nextwavefronttopropagate = wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes[0]  # this wave-front, in the next time-step, will either reach the next node on its shortest path to an exit or will first intersect with a slower-moving wave-front propagating on the same link (in which case both wave-fronts are probably upstream-propagating shock-fronts; TODO: might not be generally true, but will assert it as a check)
    print(f"nextwavefronttopropagate, which in the next time-step will either reach the next node on its shortest path to an exit or will first intersect with a slower-moving wave-front propagating on the same link, is {nextwavefronttopropagate}")
    edgePropagatingOn = nextwavefronttopropagate['currentlink']
    b = nextwavefronttopropagate['destinationnodeonlink']
    otherWavefrontsOnSameLink = [wv for wv in propagatingwavefronts if wv['currentlink'] == edgePropagatingOn and wv != nextwavefronttopropagate]
    print(f"otherWavefrontsOnSameLink is {otherWavefrontsOnSameLink}")
#    if otherWavefrontsOnSameLink: sys.exit()
    if otherWavefrontsOnSameLink:
      speed_nwvtp = nextwavefronttopropagate['speed']
      distreacheddownlink_nwvtp = nextwavefronttopropagate['distancereacheddownlink']
#      catchuptime = {}  # for each other wave-front on the same link, the time it would be caught up to by 'nextwavefronttopropagate'
      minTimeForNextWavefrontToCatchAnotherOnSameLink = 9e99
      for owv in otherWavefrontsOnSameLink:  # calculate time till this wave-front would be caught up from behind by 'nextwavefronttopropagate'
        speed_owv = owv['speed']
        distreacheddownlink_owv = owv['distancereacheddownlink']
        assert speed_owv < speed_nwvtp  # the two wave-fronts' speeds should never be equal, unless they are both already at the same distance down the link, which shouldn't occur
        assert distreacheddownlink_owv > distreacheddownlink_nwvtp
#        assert speed_owv != speed_nwvtp  # the two wave-fronts' speeds should never be equal, unless they are both already at the same distance down the link, which shouldn't occur
        timeTillCaughtUp = (distreacheddownlink_owv - distreacheddownlink_nwvtp) / (speed_nwvtp - speed_owv)
        assert timeTillCaughtUp > 0.0
#        catchuptime[owv] = timeTillCaughtUp
        if timeTillCaughtUp < minTimeForNextWavefrontToCatchAnotherOnSameLink:
          minTimeForNextWavefrontToCatchAnotherOnSameLink = timeTillCaughtUp
          wavefrontCaughtUpFromBehind = owv
##      if catchuptime:
#      print("catchuptime is", catchuptime)
#      wavefrontCaughtUpFromBehind = min(catchuptime, catchuptime.get)  # obtain key in dictionary 'catchuptime' having the minimum value; the default for 'min' is to give the minimum key, but here the keys are ordered by the function catchuptime.get, i.e. ordered by the keys' dictionary-values (see https://docs.python.org/3/library/functions.html#min)
#      minTimeForNextWavefrontToCatchAnotherOnSameLink = catchuptime[wavefrontCaughtUpFromBehind]
      print(f"minTimeForNextWavefrontToCatchAnotherOnSameLink is {minTimeForNextWavefrontToCatchAnotherOnSameLink} hours")
      print(f"wavefrontCaughtUpFromBehind is {wavefrontCaughtUpFromBehind}")
      assert minTimeForNextWavefrontToCatchAnotherOnSameLink > 0.0
      assert nextwavefronttopropagate['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable'  # TODO: this assertion might prove to be not generally correct
      assert wavefrontCaughtUpFromBehind['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable'  # TODO: this assertion might prove to be not generally correct
      nexttimestepduration = minTimeForNextWavefrontToCatchAnotherOnSameLink
      pointReachedByNextPropagatingWaveFront = 'catch up to slower-propagating shock-wave on same link'
    else:  # there are no other wave-fronts propagating on the same link as 'nextwavefronttopropagate'
      wavefrontCaughtUpFromBehind = None
      minTimeForNextWavefrontToCatchAnotherOnSameLink = None
#      nextwavefronttoreachitsdestinationnode = nextwavefronttopropagate
      nexttimestepduration = min([wv['timetilldestinationnodeonlink'] for wv in propagatingwavefronts])  # the next time-step's duration, during which all wave-fronts will propagate; at the end of the next time-step one (or more) of the wave-fronts will intersect with the next node on its shortest path to an exit
      pointReachedByNextPropagatingWaveFront = f'reach destination-node {b} at end of link'

    print(f"nextwavefronttopropagate is {nextwavefronttopropagate}")
    if simulationclocktime + nexttimestepduration > simulDuration:
      nexttimestepduration = simulDuration - simulationclocktime
      pointReachedByNextPropagatingWaveFront = f'travel until the simulation clock-time reaches simulDuration of {simulDuration} hours'
    print(f"nexttimestepduration is {nexttimestepduration} hours")
    for (i, j) in qij:
      print(f"Current flow on link ({i}, {j}) between most downstream wave-front in that link and end of link at node {j} is qij[{i}, {j}]={qij[i, j]}.")
      numVehiclesDepartingLinkDuringThisTimestep = qij[i,j] * nexttimestepduration
      print(f"numVehiclesDepartingLinkDuringThisTimestep is {numVehiclesDepartingLinkDuringThisTimestep}")
      cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j] += numVehiclesDepartingLinkDuringThisTimestep
      print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation[{i}, {j}] is now {cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j]}")
#    wavefrontsthatreachtheirdestinationnodesAtNextClockTick = [wv for wv in propagatingwavefronts if wv['timetilldestinationnodeonlink'] == nexttimestepduration]  # the wave-front that will soonest reach its next node on its shortest path to an exit
#    print(f"wavefrontsthatreachtheirdestinationnodesAtNextClockTick is {wavefrontsthatreachtheirdestinationnodesAtNextClockTick}")
#    assert len(wavefrontsthatreachtheirdestinationnodesAtNextClockTick) == 1  # TODO: modify this assertion if it turns out to be not always true; in fact, len(wavefrontsthatreachtheirdestinationnodesAtNextClockTick) can be any integer >= 0
#    nextwavefronttopropagate = wavefrontsthatreachtheirdestinationnodesAtNextClockTick[0]
#    print(f"nextwavefronttopropagate is {nextwavefronttopropagate}")

#    otherpropagatingwavefronts = [wv for wv in propagatingwavefronts if wv['timetilldestinationnodeonlink'] > nexttimestepduration]

    assert nextwavefronttopropagate['speed'] > 0.0  # wave-front's speed in direction of its travel
    a = nextwavefronttopropagate['sourcenodeonlink']
#    b = nextwavefronttopropagate['destinationnodeonlink']
    edge = frozenset({a, b})
    currentlink = nextwavefronttopropagate['currentlink']
#    print(f"Wave-front propagates away from node {a} on link {currentlink} at speed {nextwavefronttopropagate['speed']} km/hr, taking {nexttimestepduration} hours to reach node {b}.")
    print(f"Wave-front propagates away from node {a} on link {currentlink} at speed {nextwavefronttopropagate['speed']} km/hr, taking {nexttimestepduration} hours to {pointReachedByNextPropagatingWaveFront}.")
    wavefrontflowab = nextwavefronttopropagate['flow']
    nextwavefronttopropagate['time'] = simulationclocktime + nexttimestepduration
    lengthofcurrentlink = argsflowfunc['linklength'][ currentlink ]
    print(f"lengthofcurrentlink is {lengthofcurrentlink} km")
    if pointReachedByNextPropagatingWaveFront[:22] == 'reach destination-node':
      # The chosen next wave-front to propagate reaches its current link's end-node at the end of this time-step, so propagate it further, either up- or downstream, into the shortest-path tree:
      if nextwavefronttopropagate['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable':
#        numnodesremainingtobetraversedonpath = 'upstream-propagating shock-front, so not applicable'
        numnodesremainingtobetraversedonpath = 0
      else: 
        indexInPathOfCurrentNode = nextwavefronttopropagate['nodesonpathtoanexit'].index(b)
        numnodesremainingtobetraversedonpath = len(nextwavefronttopropagate['nodesonpathtoanexit']) - 1 - indexInPathOfCurrentNode
      if numnodesremainingtobetraversedonpath > 0:  # update downstream-propagating wave-front as it passes through node a to the next link on its path
#        gotonextnodeonpath()
        c = nextwavefronttopropagate['nodesonpathtoanexit'][ indexInPathOfCurrentNode + 1]
#        print(f"Next onward node on shortest path is c = {c}")
        nextlinkonpath = frozenset({b, c})
#        print(f"nextlinkonpath is {nextlinkonpath}")
        nextlinkcapacity = argsflowfunc['linkcapacity'][nextlinkonpath]
        print(f"At node {b}, wave-front shall encounter capacity on next link {nextlinkonpath} of {nextlinkcapacity} veh/hr.")
        # Update wavefront's properties on the next downstream link:
        nextwavefronttopropagate['currentlink'] = nextlinkonpath
        nextwavefronttopropagate['sourcenodeonlink'] = b
        nextwavefronttopropagate['destinationnodeonlink'] = c
        nextwavefronttopropagate['distancereacheddownlink'] = 0.0
        nextwavefronttopropagate['speed'] = argsflowfunc['lK'][frozenset({b,c})]
        lengthofnextlink = argsflowfunc['linklength'][ nextlinkonpath ]
        print(f"lengthofnextlink is {lengthofnextlink} km")
        nextwavefronttopropagate['timetilldestinationnodeonlink'] = lengthofnextlink / nextwavefronttopropagate['speed']
        nextwavefronttopropagate['flowaheadofwavefront'] = 0.0
        nextwavefronttopropagate['densityaheadofwavefront'] = 0.0
        if nextlinkcapacity >= wavefrontflowab:
          print(f"This next link's capacity of {nextlinkcapacity} is no lower than the wave-front's flow of {nextwavefronttopropagate['flow']} veh/hr, so no shock-front or bottleneck results.")
          # The wave-front traverses node b and enters the next link, perhaps changing its flow, speed, and density:
#          freeflowinnextlink = nextlinkcapacity
          # Now that wave-front has reached the end of its current link, update qij, densityij, and speedij as they are at the link's departing-end:
          qij[a, b] = wavefrontflowab  # flow at the end of the current link
#          assert len( [wv for wv in propagatingwavefronts if wv['currentlink']==edge] ) == 1  # TODO: this assertion might prove to be not generally correct
          densityij[edge] = nextwavefronttopropagate['density']
          speedij[edge] = nextwavefronttopropagate['speed']
          flowstateijatnodej[a, b] = 'freeflowstate'
          print(f"qij[{a}, {b}] is now {qij[a, b]} veh/hr, densityij[{edge}] is {densityij[edge]}, speedij[{edge}] is {speedij[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
        else:
          print(f"This next link's capacity of {nextlinkcapacity} is lower than the wave-front's flow of {nextwavefronttopropagate['flow']} veh/hr, so a bottleneck results: an upstream-propagating shock-front is created, and the wave-front's flow will be reduced as it travels down the next link.")
          nextwavefronttopropagate['flow'] = nextlinkcapacity
#          congestedflowab = qij[a, b] * totalflowleavingnodeiDownstreamOfCongestion / (initialtotalflowEnteringNodeiFromItsParents + injectionflowati)  # TODO: divide congested flow between the links entering node b: the flow in each link, and also by implication the injection-flow, decreases proportionally
          congestedflowab = nextlinkcapacity
          print(f"congested flow[{a}, {b}] is {congestedflowab} veh/hr")
          print(f"As congestion-zone propagates upstream, flow in link [{a}, {b}] will decrease from {wavefrontflowab} to {congestedflowab} veh/hr")
          # Now that wave-front has reached the end of its current link, update qij, densityij, and speedij as they are at the link's departing-end:
          qij[a, b] = congestedflowab
#          print(f"edge is {edge}, propagatingwavefronts is {propagatingwavefronts}")
#          assert len( [wv for wv in propagatingwavefronts if wv['currentlink']==edge] ) == 1  # TODO: this assertion might prove to be not generally correct
#          totalflowEnteringNodei += congestedflowab
          edge = frozenset({a, b})
          congesteddensityab = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowab, edge, argsflowfunc )
#          densityij[edge] = congesteddensityab
#          print(f"densityij IS NOW {densityij} veh/km")
          densityij[edge] = congesteddensityab
          speedij[edge] = congestedflowab / congesteddensityab
          flowstateijatnodej[a, b] = 'congestedflowstate'
          print(f"qij[{a}, {b}] is now {qij[a, b]} veh/hr, densityij[{edge}] is {densityij[edge]}, speedij[{edge}] is {speedij[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
          upstreamfreeflowab = wavefrontflowab
          upstreamfreedensityab = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamfreeflowab, edge, argsflowfunc )
          print(f"Calculating propagation-speed of shock-front: upstreamfreeflowab is {upstreamfreeflowab}, upstreamfreedensityab is {upstreamfreedensityab}, congestedflowab is {congestedflowab}, congesteddensityab is {congesteddensityab}")
          shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamfreeflowab, upstreamfreedensityab, congestedflowab, congesteddensityab )
          if shockFrontPropagationSpeed < 0.0:
            upordownstreamText = 'travels upstream'
          elif shockFrontPropagationSpeed > 0.0:
            upordownstreamText = 'travels downstream'
          else:
            upordownstreamText = 'is stationary'
          print(f"shockFrontPropagationSpeed is {shockFrontPropagationSpeed} km/hr: the front {upordownstreamText}")
          assert shockFrontPropagationSpeed < 0.0
#          if shockFrontPropagationSpeed == 0.0:  # can be zero if, for example, the link's free flow upstream of the shock-front is zero, as then the "congested" flow is also zero
#            print(f"upstreamfreeflowab {upstreamfreeflowab} = congestedflowab {congestedflowab}: there is no real congestion.")
#          print(f"Density in link {edge} is now {densityij[edge]}, speed is now {speedij[edge]}.")
          timetillupstreamnodeonlink = -lengthofcurrentlink / shockFrontPropagationSpeed
          assert timetillupstreamnodeonlink > 0.0
          shockfrontdict = {'time': simulationclocktime + nexttimestepduration, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': b, 'destinationnodeonlink': a, 'distancereacheddownlink': 0.0, 'speed': -shockFrontPropagationSpeed, 'flow': congestedflowab, 'density': congesteddensityab, 'timetilldestinationnodeonlink': timetillupstreamnodeonlink, 'nodesonpathtoanexit': 'upstream-propagating shock-front, so not applicable', 'flowaheadofwavefront': upstreamfreeflowab, 'densityaheadofwavefront': upstreamfreedensityab}

          print(f"New upstream-propagating shock-front is {shockfrontdict}")
          propagatingwavefronts.append( shockfrontdict )
#          print(f"After appending new upstream-propagating shock-front, propagatingwavefronts is {propagatingwavefronts}")
          nextwavefronttopropagate['density'] = nextwavefronttopropagate['flow'] / nextwavefronttopropagate['speed']

      else:  # there are no more untraversed links remaining on this shortest path
        if nextwavefronttopropagate['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable':  # the wave-front is upstream-propagating
          parentnodesofb = parentsOnShortestPathsFromInjectionNodes[b]
          print(f"Node {b}'s parent-nodes on shortest paths from injection-nodes are {parentnodesofb}")
          if not parentnodesofb:  # b has no parent-nodes on the shortest paths
            nextwavefronttopropagate['timetilldestinationnodeonlink'] = 0.0
            nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['speed'] * nexttimestepduration
#            assert nextwavefronttopropagate['distancereacheddownlink'] == argsflowfunc['linklength'][ frozenset({a, b}) ]
            assert abs(nextwavefronttopropagate['distancereacheddownlink'] - lengthofcurrentlink) <= EPS
            print(f"As it exits the network at exit-node {b}, wave-front {nextwavefronttopropagate} is removed from propagatingwavefronts.")
            propagatingwavefronts.remove( nextwavefronttopropagate )
          else:  # update upstream-propagating shock-wave as it travels up node b's parent-links
#            assert len(parentsOnShortestPathsFromInjectionNodes[b]) == 1  # TODO: generalise to the case in which node b has more than one parent-link on the shortest paths from injection-nodes
#            p0 = parentsOnShortestPathsFromInjectionNodes[b][0]  # p0 is the single parent-node of b
            for p in parentnodesofb:
#              flowenteringparent = flowEnteringNode( p, i, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes )
              assert b == shortestPathToAnExit[p][1]  # b is the next node in the shortest path from p to an exit-node
#              linkcapacity = argsflowfunc['linkcapacity'][frozenset({p,i})]
#              if flowenteringparent > linkcapacity:
#                print(f"WARNING: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {i}) which has a capacity of {linkcapacity:.3f}: the link will have flow at this capacity, while a congestion-zone will propagate up the upstream link(s); the congestion-zone will have density to the right of the hump in its link's flow-function, and correspondingly low speed.")
              flowpb = qij[p, b]
#              densitypb = densityij[frozenset({p, b})]
#              speedpb = flowpb / densitypb
              timeThatShockFrontReachedNodea = simulationclocktime + nexttimestepduration
#              propagateSEM4congestionUpstreamToParents( nextwavefronttopropagate, b, p, flowpb, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )
              propagateSEM4congestionUpstreamToParents( nextwavefronttopropagate, b, timeThatShockFrontReachedNodea, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, flowstateijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes, propagatingwavefronts, nexttimestepduration )

        else:  # the wave-front is downstream-propagating
          assert b in exitnodes
          nextwavefronttopropagate['timetilldestinationnodeonlink'] = 0.0
          nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['speed'] * nexttimestepduration
          print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']}, lengthofcurrentlink is {lengthofcurrentlink}")
#          assert nextwavefronttopropagate['distancereacheddownlink'] == argsflowfunc['linklength'][ frozenset({a, b}) ]
          assert abs(nextwavefronttopropagate['distancereacheddownlink'] - lengthofcurrentlink) <= EPS
          qij[a, b] = nextwavefronttopropagate['flow']  # wave-front has reached the end of its current link, so update qij, densityij, and speedij
          densityij[edge] = nextwavefronttopropagate['density']
          speedij[edge] = nextwavefronttopropagate['speed']
          flowstateijatnodej[a, b] = 'freeflowstate'
          print(f"qij[{a}, {b}] is now {qij[a, b]} veh/hr, densityij[{edge}] is {densityij[edge]}, speedij[{edge}] is {speedij[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
          print(f"As it exits the network at exit-node {b}, wave-front {nextwavefronttopropagate} is removed from propagatingwavefronts.")
          propagatingwavefronts.remove( nextwavefronttopropagate )
          assert nextwavefronttopropagate['nodesonpathtoanexit'][-1] == b
          outflowsByNodeID[b] = -nextwavefronttopropagate['flow']

    elif pointReachedByNextPropagatingWaveFront == 'catch up to slower-propagating shock-wave on same link':
      nextwavefronttopropagate['timetilldestinationnodeonlink'] -= nexttimestepduration
      nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['speed'] * nexttimestepduration
      print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']}, lengthofcurrentlink is {lengthofcurrentlink}")
      upstreamflow = wavefrontCaughtUpFromBehind['flowaheadofwavefront']
      upstreamdensity = wavefrontCaughtUpFromBehind['densityaheadofwavefront']
      congestedflow = nextwavefronttopropagate['flow']
      congesteddensity = nextwavefronttopropagate['density']
      print(f"Calculating new propagation-speed of shock-front: upstreamflow is {upstreamflow}, upstreamdensity is {upstreamdensity}, congestedflow is {congestedflow}, congesteddensity is {congesteddensity}")
      shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamflow, upstreamdensity, congestedflow, congesteddensity )
      if shockFrontPropagationSpeed < 0.0:
        upordownstreamText = 'travels upstream'
      elif shockFrontPropagationSpeed > 0.0:
        upordownstreamText = 'travels downstream'
      else:
        upordownstreamText = 'is stationary'
      print(f"New shockFrontPropagationSpeed is {shockFrontPropagationSpeed} km/hr: the front {upordownstreamText}")
      assert shockFrontPropagationSpeed < 0.0
      nextwavefronttopropagate['flowaheadofwavefront'] = wavefrontCaughtUpFromBehind['flowaheadofwavefront']
      nextwavefronttopropagate['densityaheadofwavefront'] = wavefrontCaughtUpFromBehind['densityaheadofwavefront']
      print(f"As it's caught up from behind by upstream-propagating shock-front with speed {nextwavefronttopropagate['speed']}, wave-front {wavefrontCaughtUpFromBehind} is removed from propagatingwavefronts.")
      propagatingwavefronts.remove( wavefrontCaughtUpFromBehind )
      print(f"propagatingwavefronts is {propagatingwavefronts}")

    elif pointReachedByNextPropagatingWaveFront[:60] == 'travel until the simulation clock-time reaches simulDuration':
      nextwavefronttopropagate['timetilldestinationnodeonlink'] -= nexttimestepduration
      nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['speed'] * nexttimestepduration
      print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']} km")
      print(f"As its propagation ceases at the end of the simulation, wave-front {nextwavefronttopropagate} is removed from propagatingwavefronts.")
      propagatingwavefronts.remove( nextwavefronttopropagate )

    if nextwavefronttopropagate in propagatingwavefronts:  # if this wave-front hasn't just been removed from 'propagatingwavefronts'
#      print(f"Updated downstream-propagating wavefront is {nextwavefronttopropagate}")
      print(f"Updated propagating wavefront is {nextwavefronttopropagate}")
##    print(f"After updating downstream-propagating wavefront, propagatingwavefronts is {propagatingwavefronts}")
#    print(f"After updating downstream-propagating wavefront, len(propagatingwavefronts) is {len(propagatingwavefronts)}")
    print(f"After updating propagating wavefront, len(propagatingwavefronts) is {len(propagatingwavefronts)}")

#    otherpropagatingwavefronts = [wv for wv in propagatingwavefronts if wv['timetilldestinationnodeonlink'] > nexttimestepduration]
    otherpropagatingwavefronts = [wv for wv in propagatingwavefronts if wv['time'] == simulationclocktime]
    print(f"otherpropagatingwavefronts, i.e. with 'time' {simulationclocktime}, is {otherpropagatingwavefronts}")
#    assert len(otherpropagatingwavefronts) == len(propagatingwavefronts) - 1  # not generally true
    for wavefront in otherpropagatingwavefronts:
#      if wavefront['speed'] <= 0.0:
#        continue
      wavefront['time'] += nexttimestepduration
      wavefront['distancereacheddownlink'] += wavefront['speed'] * nexttimestepduration
      wavefront['timetilldestinationnodeonlink'] -= nexttimestepduration
      if wavefront['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable':
        direction = 'upstream'
      else:
        direction = 'downstream'
      print(f"Updated {direction}-propagating wavefront is {wavefront}")
  
    simulationclocktime += nexttimestepduration
    print(f"==============simulationclocktime is {simulationclocktime} hours ({simulationclocktime*3600} seconds).")
    print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation is {cumulNumVehiclesToHaveDepartedLinkDuringSimulation}")
#    sys.exit()

  if simulationclocktime < simulDuration:  # run final time-step, ending at the specified duration of the simulation
    assert propagatingwavefronts == []
    finaltimestepduration = simulDuration - simulationclocktime  # the final time-step's duration, by the start of which every propagating wave-front has departed the network: this implies that during the final time-step, the flow in (the entire length of) every directed link (i, j) adds to cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j]
    print(f"finaltimestepduration is {finaltimestepduration} hours")
    for (i, j) in qij:
      print(f"Current flow on entire length of link ({i}, {j}) is qij[{i}, {j}]={qij[i, j]}.")
      numVehiclesDepartingLinkDuringThisTimestep = qij[i,j] * finaltimestepduration
      print(f"numVehiclesDepartingLinkDuringThisTimestep on link ({i}, {j}) is {numVehiclesDepartingLinkDuringThisTimestep}")
      cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j] += numVehiclesDepartingLinkDuringThisTimestep
      print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation[{i}, {j}] is now {cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j]}")
    simulationclocktime += finaltimestepduration
    print(f"==============simulationclocktime is {simulationclocktime} hours ({simulationclocktime*3600} seconds).")

#  timeaveragedflowsinlinksoverSimulDuration = {}
#  return timeaveragedflowsinlinksoverSimulDuration
  return (cumulNumVehiclesToHaveDepartedLinkDuringSimulation, outflowsByNodeID)


#def runSEM2or3or4( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, simulDuration=14/3600):
def runSEM2or3or4( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, simulDuration=1.0):  # SEM4 specifies a duration (in hours) for the simulated propagation of wave- and shock-fronts
  starttime = time.time()
  print(f"duration of simulation is {simulDuration*3600} seconds.")
  print("inflowsByNodeID is", inflowsByNodeID)
##  (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0, pointcoordsbyID) = setuptrafficflowproblem( JSONnetworkfilename, b0 )
#  (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, argsflowfunc, pointcoordsbyID) = setuptrafficflowproblem( JSONnetworkfilename )
  (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, argsflowfunc, pointcoordsbyID, linestringsWithCoords) = setuptrafficflowproblem( JSONnetworkfilename, exitnodes )

  print("injectionnodes is", injectionnodes)
  b0 = np.zeros( len(injectionnodes) )
#  b0[3] = 1000.0  # TODO: should be determined via 'inflowsByNodeID'
#  for i in range(N):  # node N-1 won't have any neighbours not already considered
#  for i in injectionnodes:
  for i in inflowsByNodeID:
    b0[indexamonginjectionorexitnodes[i]] = inflowsByNodeID[i]
  print("b0 (inflow at injection-nodes) is", b0)
  assert len(b0) == len(injectionnodes)

  for i in [j for j in range(N) if j not in inflowsByNodeID]:
    inflowsByNodeID[i] = 0.0
  print("inflowsByNodeID is", inflowsByNodeID)

#  SPEEDOFINJECTEDTRAFFIC = 10  # assumed speed of injected traffic at every injection-node (speed of traffic as it enters an injection-node), in kilometres per hour
#  SPEEDOFINJECTEDTRAFFIC = 13  # assumed speed of injected traffic at every injection-node (speed of traffic as it enters an injection-node), in kilometres per hour
#  SPEEDOFINJECTEDTRAFFIC = 20  # assumed speed of injected traffic at every injection-node (speed of traffic as it enters an injection-node), in kilometres per hour
#  SPEEDOFINJECTEDTRAFFIC = 30  # assumed speed of injected traffic at every injection-node (speed of traffic as it enters an injection-node), in kilometres per hour

#  v0 = np.full( shape=len(b0), fill_value=SPEEDOFINJECTEDTRAFFIC)  # speed of traffic as it enters each injection-node, in kilometres per hour
#  print("v0 (speed of traffic as it enters an injection-node) is", v0)
  totalinflow = sum(b0)
  print("totalinflow is sum(b0) =", totalinflow)
  endtime = time.time()
  timetodefineproblemandnetworkandprocessgraph = endtime - starttime
  print("Definition of problem and network, and processing of graph, took %.5f seconds." % timetodefineproblemandnetworkandprocessgraph)
#  sys.exit()  # for debugging


  if SEMversion == 'SEM2':  # use the solver, rather than finding shortest-path routes
    starttime = time.time()
#    NEWTONRAPHSON = True  # use Newton-Raphson to find h satisfying f(h) = 0
    NEWTONRAPHSON = False  # use a (local- or global-) minimisation method to find the minimum of ||f(h)||^2
    ENERGYFUNCTION = 'conservationOfMassResidual'
#    ENERGYFUNCTION = 'conservationOfMassResidualPlusTotalCongestion'
#    ENERGYFUNCTION = 'conservationOfMassResidualLessTotalCongestion'

#     The solver, whether Newton-Raphson or global minimiser, can be run more than once in order to automatically find distinct solutions h.
    NUMINITPOINTSPERDIMENSION = 1  # run the solver, either Newton-Raphson or global minimiser, exactly once

#    NUMINITPOINTSPERDIMENSION = 6  # run the solver more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 11  # run the solver more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 21  # run the solver more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 101  # run the solver more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 1001  # run the solver more than once in order to automatically find distinct solutions h
##    hJmin = int(-4)
##    hJmax = int(0)
    hJmin = int(0)
#    hJmax = int(4)
#    hJmax = int(10)
#    hJmax = int(20)
#    hJmax = int(600)
    hJmax = int(10000)

#    if NUMGLOBALOPTIMISERRUNS > 1:
    if NUMINITPOINTSPERDIMENSION > 1:
#      if NUMGLOBALOPTIMISERRUNS < 2**len(injectionnodes):  # ensure that at each injection-node r, h_r takes on at least two different values
#        NUMGLOBALOPTIMISERRUNS_orig = NUMGLOBALOPTIMISERRUNS
#        NUMGLOBALOPTIMISERRUNS = 2**len(injectionnodes)
#        print("Increasing NUMGLOBALOPTIMISERRUNS from %d to %d" % (NUMGLOBALOPTIMISERRUNS_orig, NUMGLOBALOPTIMISERRUNS))
#      numpointsperdimension = int(round( NUMGLOBALOPTIMISERRUNS**(1/len(injectionnodes)) ))
#      print("numpointsperdimension is %d" % numpointsperdimension)
#      deltah = (hJmax-hJmin) / (NUMGLOBALOPTIMISERRUNS**(1/len(injectionnodes)) - 1)
      deltah = (hJmax-hJmin) / (NUMINITPOINTSPERDIMENSION - 1)
      print("deltah is %.5f" % deltah)
      hinitList = [tuple([deltah * i for i in il]) for il in generateAllIndices( len(injectionnodes), NUMINITPOINTSPERDIMENSION ) ]
      print("hinitList is %s" % hinitList)

#    else:  # NUMGLOBALOPTIMISERRUNS is 1
    else:  # NUMINITPOINTSPERDIMENSION is 1
      hinitInjectionNodes = len(injectionnodes) * [0.0]
      # Override the default zero-vector initial estimate:
#      hinitInjectionNodes = [90]
#      hinitInjectionNodes = [0, 6]
#      hinitInjectionNodes = [0, 10]
#      hinitInjectionNodes = [0, 12]
#      hinitInjectionNodes = [0.11791782, 0.04704803]
#      hinitInjectionNodes = [0.2, 0.2]
#      hinitInjectionNodes = [1.6, 0.25]
#      hinitInjectionNodes = [4, 2]
#      hinitInjectionNodes = [6, 2]
#      hinitInjectionNodes = [0.0, 1.0, 2.5]
#      hinitInjectionNodes = [0.0, 2.0, 2.0]
#      hinitInjectionNodes = [1.0, 3.0, 3.0]

#      hinitInjectionNodes = [10*random.random(), 10*random.random()]  # choose an initial estimate randomly (when trying to find several distinct solutions)
#      hinitInjectionNodes = [3*random.random(), 3*random.random()]  # choose an initial estimate randomly (when trying to find several distinct solutions)
#      hinitInjectionNodes = [5*random.random(), 5*random.random()]  # choose an initial estimate randomly (when trying to find several distinct solutions)

      hinitList = [ hinitInjectionNodes ]


    flowfunction = flowijfunc[flowfunc]
    derivfunction = derivijfunc[flowfunc]
    solutionsfound = []
    NUMGLOBALOPTIMISERRUNS = NUMINITPOINTSPERDIMENSION**len(injectionnodes)
    if NEWTONRAPHSON:
      maxNRiterationsUsedToFindGoodSoln = 0
    endtime = time.time()
    timetoinitialisestartingpoints = endtime - starttime

    sumofsolvertimes = 0.0
    for w in range(NUMGLOBALOPTIMISERRUNS):
      print("GLOBALOPTIMISER-RUN NUMBER %d:" % w)
      hinitInjectionNodes = hinitList[w]
      print("hinitInjectionNodes is %s" % repr(hinitInjectionNodes))
      starttime = time.time()
      (exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, numNRiterationsused, outflowsByNodeID) = runsolver(NEWTONRAPHSON, N, injectionnodes, exitnodes, indexamonginjectionorexitnodes, terminatorPressure, hExitNodes, neighboursofj, flowfunction, derivfunction, hJmin, hJmax, hinitInjectionNodes, argsflowfunc, ENERGYFUNCTION, flowijfunc, flowfunc, b0)
      endtime = time.time()
      print("congestionij is %s" % congestionij)  # these two lines check whether congestionij and totalnetworkcongestion are correctly returned by runsolver()
      print("Total congestion in network is %.3f" % totalnetworkcongestion)
      timetofindsoln = endtime - starttime
      sumofsolvertimes += timetofindsoln
      print("Solver-run %d: hInjectionNodes is %s" % (w, hInjectionNodes))
#      print("flowijfunc is %s" % flowijfunc)
      print("flowfunc %s, ENERGYFUNCTION %s, b0 %s" % (flowfunc, ENERGYFUNCTION, b0) )
#      print("argsflowfunc is %s" % argsflowfunc)
##      totalnetworkcongestion = totalcongestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )
#      (congestionij, totalnetworkcongestion) = congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
#      print("congestionij is %s" % congestionij)
#      print("Total congestion in network is %.3f" % totalnetworkcongestion)
      solutionsfound.append( (hinitInjectionNodes, exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, timetofindsoln) )
      if flowbalanced and NEWTONRAPHSON:
#        maxNRiterationsUsedToFindGoodSoln = max(numNRiterationsused, maxNRiterationsUsedToFindGoodSoln)
        if numNRiterationsused > maxNRiterationsUsedToFindGoodSoln:
          maxNRiterationsUsedToFindGoodSoln = numNRiterationsused
#          print("maxNRiterationsUsedToFindGoodSoln is now %d" % maxNRiterationsUsedToFindGoodSoln)
    print("All %d runs of the solution-finding method took %.5f seconds." % (NUMGLOBALOPTIMISERRUNS, sumofsolvertimes))
    starttime = time.time()
##    print("solutionsfound is %s" % solutionsfound)
#    flowbalancedsolutions = []
    validsolutions = []
    nonsolutionsfound = 0  # count the number of times the solver fails to find a valid, i.e. flow-balanced, solution
    totaltimetofindnonsolutions = 0.0  # keep track of the total time spent finding invalid, i.e. flow-unbalanced, solutions
    for w in range(NUMGLOBALOPTIMISERRUNS):
      (hinitInjectionNodes, exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, timetofindsoln) = solutionsfound[w]
#      print("Solver-run %d: hinitInjectionNodes is %s, exception is %s, flowbalanced is %s, hInjectionNodes is %s, qij is %s" % (w, hinitInjectionNodes, exception, flowbalanced, hInjectionNodes, qij))
#      print("qij is %s" % qij)
      if qij == None:
        nonsolutionsfound += 1
        totaltimetofindnonsolutions += timetofindsoln  # add to the total time spent finding this flow-balanced solution
        continue
      qij_nparray = np.array([ di[1] for di in sorted(qij.items()) ])
##      if flowbalanced and (hInjectionNodes >= hJmin).all() and (hInjectionNodes <= hJmax).all():  # if solution is flow-balanced and its h-values lie between the acceptable bounds, then add it to 'flowbalancedsolutions' so long as it's not too close to any previously-recorded flow-balanced solution
#      if flowbalanced:  # if solution is flow-balanced, then add it to 'flowbalancedsolutions' so long as it's not too close to any previously-recorded flow-balanced solution
#      if AUGMENTEDENERGYFUNCTION or flowbalanced and not AUGMENTEDENERGYFUNCTION:  # if solution is valid, then add it to 'validsolutions' so long as it's not too close to any previously-recorded valid solution
      if ENERGYFUNCTION in ('conservationOfMassResidualPlusTotalCongestion', 'conservationOfMassResidualLessTotalCongestion') or flowbalanced and ENERGYFUNCTION == 'conservationOfMassResidual':  # if solution is valid, then add it to 'validsolutions' so long as it's not too close to any previously-recorded valid solution
        solutionPreviouslyRecorded = False
#        for fbs in flowbalancedsolutions:
        for fbs in validsolutions:
#          if scipy.linalg.norm(hInjectionNodes - fbs['hInjectionNodes']) / scipy.linalg.norm(fbs['hInjectionNodes']) <= 4e-2:
          if scipy.linalg.norm(qij_nparray - fbs['qij']) / scipy.linalg.norm(fbs['qij']) <= 2e-3:
            solutionPreviouslyRecorded = True
            fbs['numtimesfound'] += 1  # this flow-balanced solution has been found one more time
            fbs['totaltimetofindsoln'] += timetofindsoln  # add to total time spent finding this flow-balanced solution
            break
        if not solutionPreviouslyRecorded:
#          flowbalancedsolutions.append( {'hinitInjectionNodes':hinitInjectionNodes, 'hInjectionNodes':hInjectionNodes, 'qij':qij_nparray, 'numtimesfound':1, 'totaltimetofindsoln':timetofindsoln} )  # have found this solution once so far, and have spent a total so far of 'timetofindsoln' finding it
          validsolutions.append( {'hinitInjectionNodes':hinitInjectionNodes, 'hInjectionNodes':hInjectionNodes, 'qij':qij_nparray, 'congestionij':congestionij, 'totalnetworkcongestion':totalnetworkcongestion, 'numtimesfound':1, 'totaltimetofindsoln':timetofindsoln} )  # have found this solution once so far, and have spent a total so far of 'timetofindsoln' finding it
      else:
        nonsolutionsfound += 1
        totaltimetofindnonsolutions += timetofindsoln  # add to the total time spent finding this flow-balanced solution
    endtime = time.time()
    timetopostprocess = endtime - starttime
    print("Initialisation of starting-points took %.5f seconds." % timetoinitialisestartingpoints)
    print("Post-processing of solutions took %.5f seconds.\n" % timetopostprocess)
##    flowbalancedsolutions = [(hinitInjectionNodes, hInjectionNodes, qij) for (hinitInjectionNodes, exception, hInjectionNodes, qij, flowbalanced) in solutionsfound if flowbalanced and (hInjectionNodes >= hJmin).all() and (hInjectionNodes <= hJmax).all()]
#    print("Distinct flowbalancedsolutions are\n - %s" % '\n - '.join(["%d. %s" % (i, str(fbs)) for (i, fbs) in enumerate(flowbalancedsolutions)]))
    print("Distinct validsolutions are\n - %s" % '\n - '.join(["%d. %s" % (i, str(fbs)) for (i, fbs) in enumerate(validsolutions)]))
#    totaltimetofindvalidsolutions = sum([fbs['totaltimetofindsoln'] for fbs in flowbalancedsolutions])
    totaltimetofindvalidsolutions = sum([fbs['totaltimetofindsoln'] for fbs in validsolutions])
#    print("Solver spent a total of %.4f seconds finding valid, i.e. flow-balanced, solns." % totaltimetofindvalidsolutions)
    print("Solver spent a total of %.4f seconds finding valid solns." % totaltimetofindvalidsolutions)
#    print("The solver failed %d times to find a valid, i.e. flow-balanced, solution." % nonsolutionsfound)
    print("The solver failed %d times to find a valid solution." % nonsolutionsfound)
#    print("Solver spent a total of %.4f seconds finding invalid, (flow-unbalanced) solns." % totaltimetofindnonsolutions)
    print("Solver spent a total of %.4f seconds finding invalid solns." % totaltimetofindnonsolutions)
#    assert sum([fbs['numtimesfound'] for fbs in flowbalancedsolutions]) + nonsolutionsfound == NUMGLOBALOPTIMISERRUNS
    assert sum([fbs['numtimesfound'] for fbs in validsolutions]) + nonsolutionsfound == NUMGLOBALOPTIMISERRUNS
    if NEWTONRAPHSON:
      print("The maximum number of Newton-Raphson iterations that were used to find any particular valid solution, over all starting-points, is %d." % maxNRiterationsUsedToFindGoodSoln)
#    for fbs in flowbalancedsolutions:
    for fbs in validsolutions:
      print("Mean time spent finding a solution close to %s was %.4f" % (fbs['hInjectionNodes'], fbs['totaltimetofindsoln']/fbs['numtimesfound']))
    if nonsolutionsfound > 0:
#      print("Mean time spent finding an invalid, i.e. flow-unbalanced, solution was %.4f" % (totaltimetofindnonsolutions/nonsolutionsfound))
      print("Mean time spent finding an invalid solution was %.4f" % (totaltimetofindnonsolutions/nonsolutionsfound))
    else:
#      print("No invalid, i.e. flow-unbalanced, solutions were found.")
      print("No invalid solutions were found.")


#  elif SEMversion == 'SEM3':  # don't use the solver; instead find shortest-path routes
  elif SEMversion in ('SEM3', 'SEM4'):  # don't use the solver; instead find shortest-path routes
    starttime = time.time()
    print("\nFinding shortest-path routes from all injection-nodes ...")
    vertices = range(N)
    predictedLinkTraversalTime = {}  # traversal-time as might be predicted by a driver
    edges = []
    lK = argsflowfunc['lK']  # lK stores the 'free speed' of each road-link
    linklength = argsflowfunc['linklength']
    for edge in lK.keys():  # lK stores the 'free speed' of each road-link
      (a, b) = list(edge)
      predictedLinkTraversalTime[a, b] = predictedLinkTraversalTime[b, a] = linklength[edge] / lK[edge]
      print("predictedLinkTraversalTime[%d, %d] is %s" % (a, b, predictedLinkTraversalTime[a, b]))
      edges.append( (a,b) )
      edges.append( (b,a) )
    exception = None  # no exception, such as numpy.linalg.LinAlgError, is raised by the global-optimisation solver as it's not being employed
    (pathdistance, nextvertex) = FloydWarshall( edges, vertices, predictedLinkTraversalTime )
    print("pathdistance is", pathdistance)
    print("nextvertex is", nextvertex)

    parentsOnShortestPathsFromInjectionNodes = {}  # for each node, its immediate ancestors on all those shortest paths from injection-nodes that pass through that node
    for j in range(N):
      parentsOnShortestPathsFromInjectionNodes[j] = set()

#    shortestDistanceToAnExit = np.full( shape=len(injectionnodes), fill_value=9e99  )  # initialise to "infinity" each shortest distance from a non-leaf to an exit-node
    shortestDistanceToAnExit = {}  # for each non-leaf node, its shortest distance to an exit-node
    closestExit = {}  # for each non-leaf node, its closest exit-node
    shortestPathToAnExit = {}  # for each non-leaf node, the path to its closest exit-node
#    for i in range(N):  # for each node construct the shortest path to an exit, its length, and the node-ID of that closest exit
    for i in injectionnodes:  # only for injection-nodes do we need to construct the shortest path to an exit, its length, and the node-ID of that closest exit
      for j in range(N):
        if j == i: continue
        (shortestpath, shortestpathlength) = constructShortestPath(i, j, nextvertex, predictedLinkTraversalTime)
        print("Shortest path from %d to %d is %s with length %.4f" % (i, j, shortestpath, shortestpathlength) )
        if j in exitnodes and shortestpathlength < shortestDistanceToAnExit.get(i, 9e99):
#          print("Shortest path from %d to exit-node %d is %s with length %.4f" % (i, j, shortestpath, shortestpathlength) )
          shortestDistanceToAnExit[i] = shortestpathlength
          closestExit[i] = j
          shortestPathToAnExit[i] = shortestpath
#        print("shortestpath is", shortestpath)

    qij = {}  # store initial flows on all arcs, i.e. directed edges, that are in the shortest-path tree
    flowstateijatnodej = {}  # store flow-state, either 'freestate' or 'congestedstate', on all arcs, i.e. directed edges, that are in the shortest-path tree
    for i in injectionnodes:  # proceed down each injection-node's shortest path to an exit, adding for each node after the first one on the path its parent on the path as one of that node's parents on all shortest paths from injection-nodes
      for index, j in enumerate( shortestPathToAnExit[i] ):
        if index > 0:
          parentofj = shortestPathToAnExit[i][index-1]
          parentsOnShortestPathsFromInjectionNodes[j].add( parentofj )
          qij[parentofj, j] = 0.0
          flowstateijatnodej[parentofj, j] = 'freeflowstate'
    print(f"shortestDistanceToAnExit is", shortestDistanceToAnExit)
    print("closestExit is", closestExit)
    print("shortestPathToAnExit is", shortestPathToAnExit)
    print("parentsOnShortestPathsFromInjectionNodes is", parentsOnShortestPathsFromInjectionNodes)
    print("qij is", qij)
#    sys.exit()

    # Every vehicle flows along a shortest path from its injection-node to an exit-node, and the traffic-flow in the entire network is the sum of the flows along these shortest paths. 

    ## (Possible future development) Could use a discretised version (see Treiber & Kesting, pp. 67-75) of the continuity equation to calculate the steady-state traffic-density and out-flow at each injection-point, then at other non-leaf nodes, and then at the exit-nodes. (If a node has in-degree of at least two, then its in-flows and their densities combine, and the discretised continuity equation gives the resultant out-flow and its density at the far end of the adjacent down-flow link.)
    ## The link to the next node is divided into sections; at the beginning of a section the flow and speed are known, which gives the density once traffic has flown into the section, which gives, according to v = v(rho) for the link, the speed at the beginning of the next section. The flow is conserved between adjacent sections.

#    for i in injectionnodes:
#      j = shortestPathToAnExit[i][1]
##      timeinterval = 1/3600  # in hours
##      sectionlength = v0[indexamonginjectionorexitnodes[i]] * timeinterval
#      enteringflow = b0[indexamonginjectionorexitnodes[i]]
#      enteringspeed = v0[indexamonginjectionorexitnodes[i]]
##      sectiondensity = enteringflow * timeinterval / sectionlength
#      sectiondensity = enteringflow / enteringspeed
##      print(f"After injection-flow of {enteringflow:.1f} veh/hr at node {i} has flowed towards node {j} at speed v0={v0[indexamonginjectionorexitnodes[i]]:.1f} km/hr for {timeinterval:.6f} hours, covering a distance of {sectionlength:.3f} km, the density in this section becomes {sectiondensity:.3f} veh/km.")

#    qij = {}  # store flows on all edges
    congestionij = {}  # store congestion-values on all edges
    densityij = {}  # store densities on all edges
    speedij = {}  # store densities on all edges
    for edge in lK.keys():
#      qij[edge[0], edge[1]] = 0.0
      densityij[edge] = 'linkNotUsedByTrafficInAShortestPath'
      speedij[edge] = 'linkNotUsedByTrafficInAShortestPath'
    outflowsByNodeID = {}  # store the out-flows at exit-nodes, as they are calculated by the solution-method (Newton-Raphson, global optimiser, or shortest-path-based method)

    if SEMversion == 'SEM3':
      for i in exitnodes:
#        flow, speed = flowandspeedEnteringNode( i, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnExit, flowijfunc, flowfunc )
        childnode = None  # an exit-node has no child-nodes on any shortest path
#        flow, speed = flowandspeedEnteringNode( i, childnode, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc )  # get rid of injection-speed 'v0' and replace it with the link's free speed
        flow = flowEnteringNode( i, childnode, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes )  # get rid of injection-speed 'v0' and replace it with the link's free speed
        outflowsByNodeID[i] = -flow
    elif SEMversion == 'SEM4':
#      timeaveragedflowsinlinksoverSimulDuration = evolutionofwaveandshockfronts( inflowsandflowperiodsByNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes )
      (cumulNumVehiclesToHaveDepartedLinkDuringSimulation, outflowsByNodeID) = evolutionofwaveandshockfronts( inflowsandflowperiodsByNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qij, densityij, speedij, flowstateijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnExit, flowijfunc, flowfunc, injectionnodes, outflowsByNodeID )
      timeaveragedflowsinlinksoverSimulDuration = {key: cumulNumVehiclesToHaveDepartedLinkDuringSimulation[key] / simulDuration for key in cumulNumVehiclesToHaveDepartedLinkDuringSimulation}

    print("qij is", qij)
    print("outflowsByNodeID is", outflowsByNodeID)
    if SEMversion == 'SEM3':
#      print(f"flow is {flow:.3f}, speed is {speed:.3f}.")
      print(f"flow is {flow:.3f}.")
    elif SEMversion == 'SEM4':
      print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation is {cumulNumVehiclesToHaveDepartedLinkDuringSimulation}.")
      print(f"timeaveragedflowsinlinksoverSimulDuration is {timeaveragedflowsinlinksoverSimulDuration}.")
    endtime = time.time()
    timetorunshortestpathbasedmethod = endtime - starttime
    print(f"Shortest-path-based method took {timetorunshortestpathbasedmethod:.5f} seconds.")
    print("densityij is", densityij)
    print("speedij is", speedij)
    starttime = time.time()
    print("Calculating congestion in edges...")
#    print(f"qij.keys() is {qij.keys()}")
    print(f"densityij.keys() is {densityij.keys()}")
##    for edge in qij.keys():
#    for edge in densityij.keys():
    for edge in lK.keys():
##      humpdensity = argsflowfunc['tau'][frozenset({p,i})]
##      if densityinlink > humpdensity:
##        congestiononlink = linkcapacity - abs(flowleavingparent)
      (r, s) = edge
      if (r, s) in qij:
        flowinedge = qij[r, s]
#      if flowinedge < 0.0:
      elif (s, r) in qij:
        (r, s) = (s, r)
        flowinedge = qij[r, s]
      else:
        congestionij[frozenset({r, s})] = 'linkNotUsedByTrafficInAShortestPath'
        print(f"Congestion on link ({r}, {s}) is 'linkNotUsedByTrafficInAShortestPath'.")
        continue
      print(f"flowinedge qij[{r}, {s}] is {flowinedge}")
      print(f"densityij[{edge}] is {densityij[edge]}")
      congestiononlink = linkcongestion(r, s, flowinedge, densityij[edge], argsflowfunc)
      congestionij[frozenset({r, s})] = congestiononlink
      assert congestiononlink >= 0.0
      linkcapacity = argsflowfunc['linkcapacity'][edge]
      if abs(congestiononlink) > 0.0:  # TODO: is the abs() necessary?
        print(f"CONGESTION of {linkcapacity:.3f}-{flowinedge:.3f}={congestiononlink:.3f} occurs on link ({r}, {s}) which has capacity of {linkcapacity:.3f}.")  # TODO: this measure of congestion gives a negative value when the flow entering the parent is greater than the link-capacity: might need a better measure
      else:
        print(f"There is NO congestion on link ({r}, {s}) which has capacity of {linkcapacity:.3f}.")
    endtime = time.time()
    timetocalculatecongestion = endtime - starttime
    print(f"Calculation of congestion took {timetocalculatecongestion:.5f} seconds.")
    if SEMversion == 'SEM4':
#      origkeys_qij = qij.keys() 
#      for (a,b) in origkeys_qij:
      for (a,b) in qij:
        qij[a,b] = timeaveragedflowsinlinksoverSimulDuration[a,b]
      print(f"qij is now defined to be the time-averaged flows over simulDuration of {simulDuration} hours ({simulDuration*3600:.1f} seconds), qij={qij}")



  # Now that solver has run, save the output in GeoJSON format:
  outputGeoJSON = {}  # fill this with GeoJSON data:
  BBwest = min([pointcoordsbyID[i][0] for i in range(N)])  # bounding-box limits
  BBsouth = min([pointcoordsbyID[i][1] for i in range(N)])
  BBeast = max([pointcoordsbyID[i][0] for i in range(N)])
  BBnorth = max([pointcoordsbyID[i][1] for i in range(N)])
  print(f"BBwest {BBwest}, BBsouth {BBsouth}, BBeast {BBeast}, BBnorth {BBnorth}")
  outputGeoJSON["bbox"] = [BBwest, BBsouth, BBeast, BBnorth]
  outputGeoJSON["features"] = []
  print("inflowsByNodeID is", inflowsByNodeID)
  for nodeID in sorted( pointcoordsbyID.keys() ):
    pointJSON = {}
    pointJSON["geometry"] = {}
    pointJSON["geometry"]["coordinates"] = list(pointcoordsbyID[nodeID])
    pointJSON["geometry"]["type"] = "Point"
    pointJSON["properties"] = {}
    if nodeID in injectionnodes:  # the current node is an injection-node
      pointJSON["properties"]["inflow"] = inflowsByNodeID[nodeID]  # TODO: change this? Pre-existing output-JSON uses "flow" instead of "inflow", but this is easily confused with the "flow" property of links
      if SEMversion == 'SEM2':  # Newton-Raphson or the global optimiser was used
        pointJSON["properties"]["head"] = hInjectionNodes[indexamonginjectionorexitnodes[nodeID]]
    else:  # the current node is an exit-node
      pointJSON["properties"]["inflow"] = outflowsByNodeID[nodeID]
      if SEMversion == 'SEM2':  # Newton-Raphson or the global optimiser was used
        pointJSON["properties"]["head"] = hExitNodes[indexamonginjectionorexitnodes[nodeID]]
    pointJSON["properties"]["id"] = nodeID
    pointJSON["properties"]["type"] = 0  # TODO: is the "type" property always 0 for a Point?
    pointJSON["type"] = "Feature"
    outputGeoJSON["features"].append(pointJSON)

#  for linkID in sorted( linkcoordsbyID.keys() ):
  for (linkcoords, linkproperties) in linestringsWithCoords.items():  # 'linkcoords' is a two-tuple of two-tuples
    linkJSON = {}
    linkJSON["geometry"] = {}
    linkJSON["geometry"]["coordinates"] = [list(tuple) for tuple in linkcoords]
    linkJSON["geometry"]["type"] = "LineString"
    linkJSON["properties"] = {}
#    if nodeID in injectionnodes:  # the current node is an injection-node
#    linkJSON["properties"]["K"] = linkproperties['freespeed']  # might be incorrect - is K a property of the water-flow-based monotonic-flow model?
    linkJSON["properties"]["diameter"] = linkproperties['diameter']
    (a, b) = tuple(linkproperties['endpointNodeIDs'])
    if (a,b) in qij:
      if SEMversion == 'SEM2':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = qij[a, b]
      elif SEMversion == 'SEM3':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = densityij[frozenset({a, b})] * speedij[frozenset({a, b})]
      elif SEMversion == 'SEM4':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = timeaveragedflowsinlinksoverSimulDuration[a, b]
    elif (b,a) in qij:
      if SEMversion == 'SEM2':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -qij[b, a]
      elif SEMversion == 'SEM3':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -densityij[frozenset({b, a})] * speedij[frozenset({b, a})]
      elif SEMversion == 'SEM4':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -timeaveragedflowsinlinksoverSimulDuration[b, a]
    else:
#      try:
#        linkJSON["properties"]["flow"] = -densityij[frozenset({b, a})] * speedij[frozenset({b, a})]
#      except:
#        print(f"densityij[frozenset({{{b}, {a}}})] is {densityij[frozenset({b, a})]}, speedij[frozenset({{{b}, {a}}})] is {speedij[frozenset({b, a})]}")
#        raise
      linkJSON["properties"]["flow"] = 'linkNotUsedByTrafficInAShortestPath'
#    linkJSON["properties"]["k"] =  # TODO: necessary for SEM2?
    linkJSON["properties"]["length"] = linkproperties['length']
    linkJSON["properties"]["matsim_linkID"] = linkproperties['matsim_linkID']
    linkJSON["properties"]["type"] = 1  # TODO: is the "type" property always 0 for a LineString?
    linkJSON["type"] = "Feature"
    outputGeoJSON["features"].append(linkJSON)

  outputGeoJSON["type"] = "FeatureCollection"

  print("outputGeoJSON is", outputGeoJSON)



#   Draw the network as a directed graph:
#  networktype = 'linear'
#  networktype = 'highway'
#  networktype = 'star'
  networktype = 'nontree'
  print("networktype is", networktype)
  print("pointcoordsbyID is", pointcoordsbyID)
  if exception == None:
    print()
    print("qij is %s" % qij)
    print("congestionij is %s" % congestionij)
    if SEMversion in ('SEM3', 'SEM4'):  # found shortest-path routes, rather than using the global optimiser
      print("densityij is %s" % densityij)
      print("speedij is %s" % speedij)
    G = nx.DiGraph()  #G = nx.Graph()
#    DECIMALPLACES = 4  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 3  # number of decimal places to show in values on the graph
    DECIMALPLACES = 2  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 1  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 0  # number of decimal places to show in values on the graph

#    G.add_nodes_from(range(4))
##    G.add_edges_from( [e for e in G.edges()] )
##    G.add_edges_from( qij.keys() )
#    positiveflowedges = [(i,j) if qij[i,j] >= 0.0 else (j,i) for (i,j) in qij]
    positiveflowedges = []
    elabels = {}
    for (i,j) in qij:
#      print("(i,j) (%d, %d)" % (i, j))
#      print("qij %f" % qij[i,j])
      flowrounded = round(qij[i,j], DECIMALPLACES)
      congestionrounded = round(congestionij[frozenset(set((i,j)))], DECIMALPLACES)
      print("flowrounded is %s" % flowrounded)
      if SEMversion == 'SEM2':  # used the solver, rather than finding shortest-path routes
#        hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
#        hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
        if i in exitnodes:
          hi = terminatorPressure
        else:
          hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
        if j in exitnodes:
          hj = terminatorPressure
        else:
          hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
#        print("hi %f, hj %f" % (hi, hj))
        densityonlink = abs(hi - hj)
      else:
        densityonlink = densityij[frozenset({i, j})]
      print("densityonlink is %s" % densityonlink)
      densityrounded = round(densityonlink, DECIMALPLACES)
      if SEMversion == 'SEM2':  # used the solver, rather than finding shortest-path routes
        speedonlink = qij[i,j] / densityonlink
#      else:
      elif SEMversion in ('SEM3', 'SEM4'):
        speedonlink = speedij[frozenset({i, j})]
      print("speedonlink is %s" % speedonlink)
      speedrounded = round(speedonlink, DECIMALPLACES)
#      print("speedrounded is %s" % speedrounded)
      if qij[i,j] >= 0.0:
        positiveflowedges.append( (i,j) )
        flowvaluetodisplay = flowrounded
#        elabels[i,j] = "q %.*f, \u03b2 %.*f, v %.*f" % (DECIMALPLACES, flowrounded, DECIMALPLACES, congestionrounded, DECIMALPLACES, -speedrounded)  # display q, \beta_{ij} (congestion), and speed on each edge {i,j}
      else:
        positiveflowedges.append( (j,i) )
        flowvaluetodisplay = -flowrounded
#        elabels[i,j] = "q %.*f, \u03b2 %.*f, v %.*f" % (DECIMALPLACES, -flowrounded, DECIMALPLACES, congestionrounded, DECIMALPLACES, -speedrounded)  # display q, \beta_{ij} (congestion), and speed on each edge {i,j}
      if abs(densityrounded) < 1e4:
        elabels[i,j] = "q %.*f, \u03c1 %.*f, v %.*f, \u03b2 %.*f" % (DECIMALPLACES, flowvaluetodisplay, DECIMALPLACES, densityrounded, DECIMALPLACES, speedrounded, DECIMALPLACES, congestionrounded)  # display q, rho (density), speed, and beta (congestion) on each edge
      else:
        elabels[i,j] = "q %.*f, \u03c1 %.1e, v %.*f, \u03b2 %.*f" % (DECIMALPLACES, flowvaluetodisplay, densityrounded, DECIMALPLACES, speedrounded, DECIMALPLACES, congestionrounded)  # display q, rho (density), speed, and beta (congestion) on each edge
    G.add_edges_from( positiveflowedges )
    debug("G.nodes() is %s" % G.nodes())
    debug("G.edges() is %s" % G.edges())
    debug("positiveflowedges is %s" % positiveflowedges)
#    print("elabels is %s" % elabels)
##    colourvalmap = {0: 1.0, 2: 0.5714285714285714, 3: 0.0}
#    colourvalmap = {}
#    for t in exitnodes:
#      colourvalmap[t] = 0.55  # colour terminal nodes red
#    colourvalues = [colourvalmap.get(node, 0.45) for node in G.nodes()]  # colour non-terminal nodes light blue
    colourvalues = ['lightblue' if node in exitnodes else 'yellow' for node in G.nodes()]  # colour non-terminal nodes light blue
    debug("colourvalues is %s" % colourvalues)

#    red_edges = [(0, 1), (1, 2)]
#    black_edges = [edge for edge in G.edges() if edge not in red_edges]
    inflowlabels = {}
    hvaluelabels = {}
    for i in range(N):  # node N-1 won't have any neighbours not already considered
      if i in exitnodes:
#        inflowlabels[i] = terminatorPressure  # not used at exit-nodes
#        hvaluelabels[i] = "h = %.*f" % (DECIMALPLACES, terminatorPressure)
        hvaluelabels[i] = "h = %.0f" % terminatorPressure
      else:
        inflowlabels[i] = b0[indexamonginjectionorexitnodes[i]]
        if SEMversion == 'SEM2':  # used the solver, rather than finding shortest-path routes
#          hvaluerounded = round(hInjectionNodes[indexamonginjectionorexitnodes[i]], DECIMALPLACES)
##          hvaluelabels[i] = "%.*f" % (DECIMALPLACES, hvaluerounded)
          hvaluelabels[i] = "h = %.*f" % (DECIMALPLACES, hInjectionNodes[indexamonginjectionorexitnodes[i]])
#          hvaluelabels[i] = "h = %e" % hInjectionNodes[indexamonginjectionorexitnodes[i]]
    debug("inflowlabels is %s" % inflowlabels)
    if SEMversion == 'SEM2':  # used the solver, rather than finding shortest-path routes
      debug("hvaluelabels is %s" % hvaluelabels)

#     Need to create a layout when making separate calls to draw nodes and edges.
    if networktype == 'linear':
      pos = nx.bipartite_layout(G, G.nodes(), align='horizontal')
#      pos = nx.spectral_layout(G)
#      pos = nx.planar_layout(G)
#      pos = nx.kamada_kawai_layout(G)
#      pos = nx.shell_layout(G)
#      pos = nx.spiral_layout(G)
#      pos = nx.circular_layout(G)
    elif networktype == 'highway':  # for plotting as bipartite set, put highway nodes in first subset, remaining nodes in second subset
##      print("G.nodes()[:-1] is %s" % repr(G.nodes()[:-1]) )  # incorrect: type(G.nodes) is <class 'networkx.classes.reportviews.NodeView'>
#      print("type(G.nodes.items()) is %s" % type(G.nodes.items()))
#      print("type(G.nodes()) is %s" % type(G.nodes()))
#      print("list(G.nodes(data=True)) is %s" % list(G.nodes(data=True)))
      nodesinfirstset = [n for n in G.nodes() if n < len(G.nodes())-1]
      debug("nodesinfirstset is %s" % nodesinfirstset)
#      pos = nx.bipartite_layout(G, G.nodes()[:-1], align='horizontal')
      pos = nx.bipartite_layout(G, nodesinfirstset, align='horizontal')
      pos[4][0] += .75  # ad-hoc line for five-node ``highway'' network
    elif networktype == 'star':
#      pos = nx.spring_layout(G, seed=0)  # set 'seed' for deterministic node-layout
      pos = nx.spring_layout(G)  # set 'seed' for deterministic node-layout
    elif networktype == 'nontree':
      pos = nx.spring_layout(G)  # set 'seed' for deterministic node-layout
#    print("pos is %s" % pos)
    assert len(pos) == len(pointcoordsbyID)
    assert set(pos.keys()) == set(pointcoordsbyID.keys())
    for nodeID in pointcoordsbyID:  # set the nodes' positions
      pos[nodeID] = np.array( pointcoordsbyID[nodeID] )
    print("pos is %s" % pos)
###    plt.figure(1, figsize=(9,5))
##    plt.figure(figsize=(9,5))
#    fig, ax = plt.subplots(figsize=(11,5))
#    fig, ax = plt.subplots(figsize=(8,5))
#    fig, ax = plt.subplots(figsize=(9,5))  # for report
    fig, ax = plt.subplots(figsize=(11,5))  # for presentation
    xmin, xmax, ymin, ymax = plt.axis()
#    print("xmin %.2f, xmax %.2f, ymin %.2f, ymax  %.2f" % (xmin, xmax, ymin, ymax))
#    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), node_color = colourvalues, node_size = 500)
    nx.draw_networkx_nodes(G, pos, node_color = colourvalues, node_size = 500)
#    nx.draw_networkx_labels(G, pos)
    nodelabels = {v:'v'+str(v) for v in G.nodes()}
    debug("nodelabels is %s" % nodelabels)
    nx.draw_networkx_labels(G, pos, nodelabels)

#    if networktype == 'linear':
#      pos_hvalues = {v: 1*(coords+np.array([0.0, -.005])) for v, coords in pos.items()}
#    elif networktype == 'highway':
#      pos_hvalues = {v: 1*(coords+np.array([0.0, -.105])) for v, coords in pos.items()}
#    elif networktype == 'star':
##      pos_hvalues = {v: coords+np.array([0.195, -.105]) for v, coords in pos.items()}
#      pos_hvalues = {v: coords+np.array([0.155, -.105]) for v, coords in pos.items()}
#    elif networktype == 'nontree':
###      pos_hvalues = {v: coords+np.array([0.055, -.105]) for v, coords in pos.items()}
##      pos_hvalues = {v: coords+np.array([0.002, -.010]) for v, coords in pos.items()}
#      pos_hvalues = {v: coords+np.array([0.002, -.005]) for v, coords in pos.items()}
#    print("pos_hvalues is %s" % pos_hvalues)
#    BBwest = 144.2235565185547
#    BBeast = 144.24372497558594
#    BBnorth = -36.92287826538086
#    BBsouth = -37.1002311706543
    BBwest = min([pos[i][0] for i in range(N)])  # bounding-box limits
    BBeast = max([pos[i][0] for i in range(N)])
    BBnorth = max([pos[i][1] for i in range(N)])
    BBsouth = min([pos[i][1] for i in range(N)])
    print(f"BBwest {BBwest}, BBeast {BBeast}, BBnorth {BBnorth}, BBsouth {BBsouth}")
#    print("BBeast-BBwest is", BBeast-BBwest)
    pos_hvalues = {}
#    for i in range(N):
    for v, coords in pos.items():
      if BBeast-BBwest > 0:
#        x = (coords[0]-BBwest)/(BBeast-BBwest)*(.74-.2)+.2
#        x = (coords[0]-BBwest)/(BBeast-BBwest)*(.74-.2)+BBwest
#        x = (coords[0]-BBwest)/(BBeast-BBwest)*.01+BBwest
        x = (coords[0]-BBwest) + BBwest
      else:
#        x = .5
        x = BBwest -.003
      if BBnorth-BBsouth > 0:
#        y = (coords[1]-BBsouth)/(BBnorth-BBsouth)*(.8-.2)+.24
#        y = (coords[1]-BBsouth)/(BBnorth-BBsouth)*(.8-.2)+BBsouth
        y = (coords[1]-BBsouth) + BBsouth
      else:
#        y = .5
        y = BBsouth - .003
      pos_hvalues[v] = np.array([x, y])
#      print("np.array([x, y] is", np.array([x, y]))
    print("pos_hvalues is %s" % pos_hvalues)
    if SEMversion == 'SEM2':  # used the solver, rather than finding shortest-path routes
      nx.draw_networkx_labels(G, pos_hvalues, hvaluelabels, font_color='k')

#    nx.draw_networkx_edge_labels(G, pos, edge_labels=elabels)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=elabels, font_color='g')
#    nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', arrows=True)
    nx.draw_networkx_edges(G, pos, edgelist=positiveflowedges, edge_color='b', arrows=True, arrowsize=25)
#    nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=False)

#    BBwest = min([pos[i][0] for i in range(N)])  # bounding-box limits
#    BBeast = max([pos[i][0] for i in range(N)])
#    BBnorth = max([pos[i][1] for i in range(N)])
#    BBsouth = min([pos[i][1] for i in range(N)])
#    print(f"BBwest {BBwest}, BBeast {BBeast}, BBnorth {BBnorth}, BBsouth {BBsouth}")
    for i in range(N):
      if i not in exitnodes:
#        if networktype == 'linear':
#          x = (pos[i][0]+1.58)/3.12  # for ``linear'' network
#          y = pos[i][1]+.57
#        elif networktype == 'highway':
#          x = (pos[i][0]+.86)/2.03  # for five-node ``highway'' network
##          y = pos[i][1]+.60
##          y = pos[i][1]+.10
#          y = pos[i][1]+.57
#        elif networktype == 'star':
#          x = (pos[i][0]+.86)/2.03
#          y = pos[i][1]+.57
#        elif networktype == 'nontree':
##          x = (pos[i][0]+.86)/2.03
#          x = (pos[i][0]+1.06)/2.03 * .8 + .10
#          y = (pos[i][1]+.85) / 1.85 * .8 + .12
        if BBeast-BBwest > 0:
          x = (pos[i][0]-BBwest)/(BBeast-BBwest)*(.74-.2)+.2
        else:
          x = .5
        if BBnorth-BBsouth > 0:
          y = (pos[i][1]-BBsouth)/(BBnorth-BBsouth)*(.8-.2)+.24
        else:
          y = .5
#        print("%d: x %.5f, y %.5f, inflowlabels[i] %.4f" % (i, x, y, inflowlabels[i]))
        fig.text(x, y, str(inflowlabels[i]), bbox=dict(facecolor='red', alpha=.5))
##    fig.text(.3, .8, 'Testing', size=10, transform=ax.transAxes)
#    fig.text(.3, .8, 'Testing', size=10)
    plt.show()

  return outputGeoJSON


#############################################################################
if __name__ == '__main__':  # Needed only if another file shall import this one as a module
#  DEBUG = False  # flag whether to print debugging-messages

  SEMversion = 'SEM2'  # use the basin-hopping global-optimisation solver
#  SEMversion = 'SEM3'  # find shortest-path routes: this is a static (link-volumes, travel-costs, and origin-destination demands are constant over time) deterministic (each vehicle takes with probability one its shortest path, rather than each vehicle having an estimate for each path's cost that differs from the actual cost by a random error) equilibrium (the travel-costs on all paths used from any origin to any destination are equal, while all unused paths have equal or greater costs) method.
  print("__name__ == '__main__': SEMversion is", SEMversion)

#  starttime = time.time()

  # TODO: 'JSONnetworkfilename' should be a parameter of the main() function
  JSONnetworkfilename = 'cmr_1s1d1r_network.links-finalnodedeleted.geojson'
  networktype = 'linear'

##  JSONnetworkfilename = 'sixnodenontree.geojson'
#  JSONnetworkfilename = 'sixnodenontree-node5deleted.geojson'  # node 5 deleted so as to make it an exit-node
##  b0 = np.array([407, 0, 0, 100, 0])  # the inflow at each injection-node is proportional to the population there, and is measured in vehicles per hour
#  b0 = np.array([500, 0, 0, 100, 0])
##  b0 = np.array([1900, 0, 0, 100, 0])
#  networktype = 'nontree'  # Setting 'networktype' determines how the graph/network will be displayed by the matplotlib library NetworkX; if networktype = 'linear' all nodes will be shown on the same horizontal line in accordance with the nx.bipartite_layout layout, while if networktype = 'star' then the nodes will be randomly placed in accordance with the nx.spring_layout layout.

#  JSONnetworkfilename = 'twonodenetwork.geojson'

#  JSONnetworkfilename = 'twonodenetwork-duplicateEndpoint.geojson'

#  JSONnetworkfilename = 'twonodenetwork-duplicateFirstEndpoint.geojson'

#  JSONnetworkfilename = 'twonodenetwork-twoduplicateEndpoints.geojson'

  JSONnetworkfilename = 'twonodenetwork-missingEndpoint.geojson'
#  b0 = np.array([400.0])  # the inflow at each injection-node is proportional to the population there, and is measured in vehicles per hour
#  b0 = np.array([443.0])
#  inflowsByNodeID = {inflowNodeID: inflowAtSingleNode}  # in general, will contain all non-zero (positive) inflows and the corresponding node-IDs
  inflowsByNodeID = {0: 443.0}  # in general, will contain all non-zero (positive) inflows and the corresponding node-IDs
#  b0 = np.array([506.111])  # capacity of the single link is about 505.605, and 506.111 is the lowest value to 3 d.p. that causes the global minimiser to fail to find a flow-balanced solution
#  b0 = np.array([600.0])
  exitnodes = {1}
  networktype = 'linear'

  print(f"JSONnetworkfilename is {JSONnetworkfilename}")

#  outputGeoJSON = runSEM2or3or4( SEMversion, inputfilename, inflowsByNodeID)
  outputGeoJSON = runSEM2or3or4( SEMversion, JSONnetworkfilename, inflowsByNodeID, exitnodes)
