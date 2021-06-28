#!/usr/bin/python3

# TODO: parameter-validation for all functions.

import scipy.optimize
#print("scipy.__version__ is", scipy.__version__)  # need a later version such as 1.4.1 in order to use multivariate scipy.optimize.newton(), i.e. call it on a vector-valued function and a vector-valued initial estimate
import sys
import numpy as np
import math
import scipy.linalg
import random
import networkx as nx
import matplotlib.pyplot as plt
import time
import scipy.stats
import json
import copy

#meanVehicleLength_Km = .005  # these two shouldn't be global variables, so their definitions have been moved inside function 'setuptrafficflowproblem()'
#meanMinGapinStoppedTraffic_Km = .004
#EPS = 1e-3  # this shouldn't be a global variable, so its definition has been moved inside each of the functions propagationSpeedOfShockFront(), propagateSEM4congestionUpstreamToParents(), flowEnteringNodeSEM3(), and evolutionofwaveandshockfrontsSEM4()


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
#  currentPointID = 0
  currentIntIDofPoint = 0  # map each matsim_nodeID (a string) to a non-negative integer, for use as a row- or column-index in the adjacency matrix
  currentLinestringID = 0
  for feature in jsonobj["features"]:
    assert feature["type"] == 'Feature'
    geometry = feature["geometry"]
    if geometry['type'] == "Point":
      key = tuple( geometry['coordinates'] )
      if key not in pointsWithCoords:  # add each non-duplicate Point
####        pointsWithCoords[ key ] = currentPointID
####        currentPointID += 1
###        pointsWithCoords[ key ] = {'matsim_nodeID': feature['properties']['matsim_nodeID']}
##        pointsWithCoords[ key ] = feature['properties']['matsim_nodeID']
#        pointsWithCoords[key]['matsim_nodeID'] = feature['properties']['matsim_nodeID']
        pointsWithCoords[ key ] = {'matsimID': feature['properties']['matsim_nodeID'], 'index': currentIntIDofPoint}
        currentIntIDofPoint += 1
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
##        linestringsWithCoords[key]['diameter'] = int( feature['properties']['diameter'] )  # isn't used so far
#        linestringsWithCoords[key]['diameter'] = float( feature['properties']['diameter'] )  # isn't used so far
        linestringsWithCoords[key]['freespeed'] = float( feature['properties']['freespeed'] )
        linestringsWithCoords[key]['length'] = float( feature['properties']['length'] )
        linestringsWithCoords[key]['matsim_linkID'] = feature['properties']['matsim_linkID']
        linestringsWithCoords[key]['oneway'] = bool(int( feature['properties']['oneway'] ))  # convert value of feature['properties']['oneway'] to a Boolean: "1" becomes 'True' (and also "2", "3", "-1", etc.) and "0" becomes 'False', but "0.0" and any other string not expressing an integer will raise a ValueError. TODO: decide whether and how the 'oneway' property should affect which solutions are allowed
#        linestringsWithCoords[key]['permlanes'] = int( feature['properties']['permlanes'] )
        linestringsWithCoords[key]['permlanes'] = float( feature['properties']['permlanes'] )  # TODO: ask Leorey what a non-integer 'permlanes' value such as 1.5 means
      else:
#        print(f"A LineString with coordinates {key} already exists in linestringsWithCoords")
        numDuplicateLinestrings += 1
#  print("pointsWithCoords is", pointsWithCoords)
  print(f"pointsWithCoords has {len(pointsWithCoords)} members.")
  if numDuplicatePoints > 0:
    pluralIndicator = "s" if numDuplicatePoints > 1 else ""
    print(f"WARNING: {numDuplicatePoints} duplicate point{pluralIndicator} ignored.")
#    raise Exception(f"WARNING: {numDuplicatePoints} duplicate point{pluralIndicator} ignored.")
#  print("linestringsWithCoords is", linestringsWithCoords)
  if numDuplicateLinestrings > 0:
    pluralIndicator = "s" if numDuplicateLinestrings > 1 else ""
    print(f"WARNING: {numDuplicateLinestrings} duplicate LineStrings ignored.")
#    raise Exception(f"WARNING: {numDuplicateLinestrings} duplicate Linestring{pluralIndicator} ignored.")  # TODO: should probably comment out this line, but will wait for more testing to see whether the exception is ever raised
##  exitnodes = []  # as of 16th April 2021, the previously-defined exit-nodes are an input-parameter to this function ('geoJSONtoAdjacencyMatrix')
  for l in linestringsWithCoords:  # add any missing Points at the ends of LineStrings
    for xy in l:
      if xy not in pointsWithCoords:
        print(f"Endpoint {xy} of LineString {l} is not in pointsWithCoords: adding it")
        pointsWithCoords[ xy ] = currentPointID  # no longer works, so cannot have a missing Point that appears in a LineString
#        exitnodes.append( currentPointID )  # include the missing Point as an exit-node: all non-missing-point exit-nodes are already specified as a parameter of this function ('geoJSONtoAdjacencyMatrix')  # TODO: confirm that this is what is generally desired
        exitnodes.add( currentPointID )  # include the missing Point as an exit-node: all non-missing-point exit-nodes are already specified as a parameter of this function ('geoJSONtoAdjacencyMatrix')  # TODO: confirm that this is what is generally desired
        print(f"WARNING: including originally-missing endpoint (PointID {currentPointID}) {xy} of LineString {l} as an exit-node: exitnodes is now {exitnodes}")
        currentPointID += 1
      else:
#        print(f"Point {xy} is in pointsWithCoords")
        pass
#  print("After addition of any missing points, pointsWithCoords is", pointsWithCoords)
  N = len(pointsWithCoords)
  print(f"After addition of any missing points, pointsWithCoords has {N} members.")
  adjacency = np.zeros( shape=(N,N) )  # construct adjacency matrix (used only by SEM2)
#  pointcoordsbyID = {}
  pointcoordsandindexbyID = {}
  pointNodeIDbyIndex = {}
#  indexofnode = 0  # map each matsim_nodeID (a string) to a non-negative integer, for use as a row- or column-index in the adjacency matrix
#pointsWithCoords[ key ] = {'matsimID': feature['properties']['matsim_nodeID'], 'index': currentIntIDofPoint}
#  for (xy, pointID) in pointsWithCoords.items():
  for (xy, pointDict) in pointsWithCoords.items():
###    pointcoordsbyID[pointID] = xy
##    pointcoordsbyID[pointDict['matsimID']] = xy
#    pointcoordsbyID[pointDict['matsimID']] = {'coords': xy, 'index': indexofnode}
#    pointcoordsandindexbyID[pointDict['matsimID']] = {'coords': xy, 'index': indexofnode}
    pointcoordsandindexbyID[pointDict['matsimID']] = {'coords': xy, 'index': pointDict['index']}
    pointNodeIDbyIndex[pointDict['index']] = pointDict['matsimID']
#    indexofnode += 1
#  linkcoordsbyID = {}  # probably don't need, as 'linestringsWithCoords' contains all information about links
#  if flowfunc == 'xExpnegx':  # TODO: must eventually generalise to allow any (non-monotonic in density) flow-function
  freespeed = {}  # free speed per link, in kilometres per hour
  permlanes = {}  # number of lanes per link
  linklength = {}  # length of each link, in kilometres
  linkcapacity = {}  # capacity of each link, in vehicles per hour
  matsimlinkID = {}  # MATSim link-ID of each link
  oneway = {}  # one-way status of each link, either 'True' or 'False'
  for l in linestringsWithCoords:
#    a = pointsWithCoords[ l[0] ]
#    b = pointsWithCoords[ l[1] ]
    a = pointsWithCoords[ l[0] ]['index']
    b = pointsWithCoords[ l[1] ]['index']
#    oneway[ frozenset({a,b}) ] = linestringsWithCoords[l]['oneway']
    oneway[(a,b)] = linestringsWithCoords[l]['oneway']
    adjacency[a,b] = 1
    if not oneway[(a,b)]:  # if {a,b} is not a one-way edge, then also set the (b,a) arc's 'adjacency' element to 1
      adjacency[b,a] = 1
    freespeed[ frozenset({a,b}) ] = linestringsWithCoords[l]['freespeed'] * 3600/1000  # convert from metres per second to kilometres per hour
    permlanes[ frozenset({a,b}) ] = linestringsWithCoords[l]['permlanes']
    linklength[ frozenset({a,b}) ] = linestringsWithCoords[l]['length'] / 1000  # convert from metres to kilometres
    linkcapacity[ frozenset({a,b}) ] = linestringsWithCoords[l]['linkcapacity']
    matsimlinkID[ frozenset({a,b}) ] = linestringsWithCoords[l]['matsim_linkID']
    linestringsWithCoords[l]['endpointNodeIDs'] = [a, b]  # add to each link in 'linestringsWithCoords' a list of its Points' node-IDs
#  print("linestringsWithCoords is", linestringsWithCoords)
  print(f"linestringsWithCoords has {len(linestringsWithCoords)} members")
#  return (adjacency, exitnodes, pointcoordsbyID, linestringsWithCoords, freespeed, permlanes, linklength, linkcapacity, matsimlinkID, oneway)
  return (adjacency, exitnodes, pointcoordsandindexbyID, pointNodeIDbyIndex, linestringsWithCoords, freespeed, permlanes, linklength, linkcapacity, matsimlinkID, oneway)


#def setuptrafficflowproblem( JSONnetworkfilename, b0 ):
def setuptrafficflowproblem( JSONnetworkfilename, exitnodes, SEMversion ):
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

#  if SEMversion in ('SEM1', 'SEM2'):  # will use Newton-Raphson or global optimiser, rather than finding shortest-path routes
  flowijfunc = {'HazenWilliams': HazenWilliamsflowfunc, 'concaveQuadratic': concaveQuadraticflowfunc, 'xExpnegx': xExpnegxflowfunc, 'triangular': triangularflowfunc}
  derivijfunc = {'HazenWilliams': HazenWilliamsderivfunc, 'concaveQuadratic': concaveQuadraticderivfunc, 'xExpnegx': xExpnegxderivfunc, 'triangular': triangularderivfunc}

#  linklength = {}  # length of each link in kilometres; entries will be added if needed

  meanVehicleLength_Km = .005
  meanMinGapinStoppedTraffic_Km = .004

#  (adjacency, exitnodes, pointcoordsbyID, linestringsWithCoords, lK, numlanes, linklength, linkcapacity, matsimlinkID, oneway) = geoJSONtoAdjacencyMatrix( JSONnetworkfilename, exitnodes )  # 'lK' is free speed
  (adjacency, exitnodes, pointcoordsandindexbyID, pointNodeIDbyIndex, linestringsWithCoords, lK, numlanes, linklength, linkcapacity, matsimlinkID, oneway) = geoJSONtoAdjacencyMatrix( JSONnetworkfilename, exitnodes )  # 'lK' is free speed
#  print("adjacency is", adjacency)
##  leafnodes = np.sum(adjacency, axis=1) == 1  # np.sum(adjacency, axis=1) is row-wise sums of 'adjacency'
##  print("leafnodes is", leafnodes)
##  leafnodeindices = [i for (i, x) in enumerate(leafnodes) if x]
##  print("leafnodeindices is", leafnodeindices)
  print("exitnodes is", exitnodes)
##  print("pointcoordsbyID is", pointcoordsbyID)
#  print("pointcoordsandindexbyID is", pointcoordsandindexbyID)
##  b0 = np.array([400])  # the inflow at each (non-leaf) node is proportional to the population there, and is measured in vehicles per hour
#  print("lK is", lK)  # 'lK' is free speed
#  print("numlanes is", numlanes)
#  print("linklength is", linklength)
#  print("linkcapacity is", linkcapacity)
#  print("matsimlinkID is", matsimlinkID)
#  print("oneway is", oneway)
  if flowfunc == 'xExpnegx':
#    CALCULATEPHI = True  # use the original method that applies when links' capacities are not known in advance
    CALCULATEPHI = False  # when the links' capacities are known in advance, use the capacity to calculate that link's value for phi
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
#      print("link %s: phi is %.3f, tau is %.3f, flow-capacity is %.3f" % (edge, phi[edge], tau[edge], linkcapacity[edge]))
    argsflowfunc = {'lK': lK, 'numlanes': numlanes, 'phi': phi, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength, 'oneway': oneway}

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
    argsflowfunc = {'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength, 'numlanes': numlanes, 'maximumdensityperlane': maximumdensityperlane, 'oneway': oneway}

  elif flowfunc == 'triangular':
#    rhoC = {}  # density per lane at capacity
    tau = {}  # in terms of the notation on p. 94 of Treiber and Kesting, tau is numlanes * rho_C
#    maximumdensityperlane = {}
    maximumdensityperlane = 120.0  # in vehicles per kilometre per lane; a constant over all edges - this value taken from Table 8.1, p. 93, 'Traffic Flow Dynamics' (2013) by Treiber and Kesting
    print(f"maximumdensityperlane is {maximumdensityperlane:.4f}")
#    T = 1.4 / 3600  # in hours; a constant over all edges - this value taken from p. 93 of Treiber and Kesting
#    print(f"T is {T:.9f}")
    T = {}  # time-gap between consecutive vehicles
    for edge in lK.keys():
      # N.B. 'linkcapacity' is capacity flow over all lanes, NOT per lane. Similarly, 'tau' is the density over all lanes at capacity flow.
      T[edge] = numlanes[edge] / linkcapacity[edge] - 1/(lK[edge]*maximumdensityperlane)  # in hours; from eq. (8.15), p. 93 of Treiber and Kesting
#      rhoC[edge] = 1 / (lK[edge] * T + 1/maximumdensityperlane)  # density per lane at capacity; from p. 94 of Treiber and Kesting
#      tau[edge] = numlanes[edge] / (lK[edge] * T[edge] + 1/maximumdensityperlane)  # density over entire link-width at capacity; from eq. (8.16), p. 94 of Treiber and Kesting
      tau[edge] = linkcapacity[edge] / lK[edge]  # density over entire link-width at capacity; from eq. (8.13), p. 93 of Treiber and Kesting, substituting Q=C/I and rho_{free} = tau/I
#      print(f"For edge {edge}, slope of free-flow part of triangular fundamental diagram is {linkcapacity[edge]/tau[edge]}; slope of congestion part of triangular fundamental diagram is {-linkcapacity[edge]/(maximumdensityperlane*numlanes[edge]-tau[edge])}.")
##    argsflowfunc = {'hvalmax': hvalmax, 'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength}
    argsflowfunc = {'lK': lK, 'tau': tau, 'linkcapacity': linkcapacity, 'linklength': linklength, 'numlanes': numlanes, 'maximumdensityperlane': maximumdensityperlane, 'T': T, 'oneway': oneway}



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



#  for key in argsflowfunc.keys():
#    print("argsflowfunc['%s'] is %s" % (key, argsflowfunc[key]))
  N = adjacency.shape[0]  # number of nodes in road-network
  print("N is %d" % N)
##  h = np.zeros(N)  #h = np.zeros(shape=(N,1))
##  print("np.sum(adjacency, axis=1) == 1 is", np.sum(adjacency, axis=1) == 1)
##  print("np.where( np.sum(adjacency, axis=1) == 1 ) is", np.where( np.sum(adjacency, axis=1) == 1 ))
#  leafnodesBool = np.sum(adjacency, axis=1) == 1  # np.sum(adjacency, axis=1) is row-wise sums of 'adjacency'
#  print("leafnodesBool is", leafnodesBool)

###  exitnodes = np.where( np.sum(adjacency, axis=1) == 1 ).tolist()  # np.sum(adjacency, axis=1) gives row-wise sum of 'adjacency'
##  exitnodes = [i for i in range(N) if leafnodesBool[i]]
#  print("exitnodes is", exitnodes)
##  injectionnodes = [i for i in range(N) if not leafnodesBool[i]]
#  injectionnodes = [i for i in range(N) if i not in exitnodes]
  injectionnodes = [nodeID for nodeID in pointcoordsandindexbyID if nodeID not in exitnodes]
#  print("injectionnodes is", injectionnodes)
  print(f"len(injectionnodes) is {len(injectionnodes)}")
  indexamonginjectionorexitnodes = np.zeros(N, dtype=int)
#  for k, i in enumerate(list(injectionnodes)):
##    print("k %d, i %d" % (k, i))
#    indexamonginjectionorexitnodes[i] = k
  k = 0
  for nodeID in injectionnodes:
    indexamonginjectionorexitnodes[ pointcoordsandindexbyID[nodeID]['index'] ] = k
    k += 1
#  for k, i in enumerate(list(exitnodes)):
##    print("k %d, i %d" % (k, i))
#    indexamonginjectionorexitnodes[i] = k
  k = 0
  for nodeID in exitnodes:
    indexamonginjectionorexitnodes[ pointcoordsandindexbyID[nodeID]['index'] ] = k
    k += 1
#  print("indexamonginjectionorexitnodes is %s" % indexamonginjectionorexitnodes)
  if SEMversion == 'SEM2':
    neighboursThatCanFlowIntoj = {}  # 'neighboursThatCanFlowIntoj' is used only by SEM2
    neighboursThatCanReceiveFlowFromi = {}  # 'neighboursThatCanReceiveFlowFromi' is used only by SEM2
    for j in range(N):
      neighboursThatCanFlowIntoj[j] = [i for i in range(N) if adjacency[i,j] == 1]
#      print("Node %d has neighbours that can flow into it %s" % (j, neighboursThatCanFlowIntoj[j]))
#    print("neighboursThatCanFlowIntoj is %s" % neighboursThatCanFlowIntoj)
    for i in range(N):
      neighboursThatCanReceiveFlowFromi[i] = [j for j in range(N) if adjacency[i,j] == 1]
#    print("neighboursThatCanReceiveFlowFromi is %s" % neighboursThatCanReceiveFlowFromi)
  else:
    neighboursThatCanFlowIntoj = None
    neighboursThatCanReceiveFlowFromi = None
  if SEMversion in ('SEM1', 'SEM2'):  # will use Newton-Raphson or global optimiser, rather than finding shortest-path routes
#    h[exitnodes] = terminatorPressure
    hExitNodes = terminatorPressure*np.ones(len(exitnodes))
#    print("h is", h)
    print("hExitNodes is", hExitNodes)
  else:
    hExitNodes = None
  assert len(injectionnodes) > 0
  assert len(injectionnodes) < N
##  iteration = 0
#  numcallsoffuncf = 0
#####  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, networktype, adjacency, b0, N, leafnodesBool, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0)
####  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, networktype, adjacency, b0, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0)
###  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0, pointcoordsbyID)
##  return (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, argsflowfunc, pointcoordsbyID)
#  return (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, argsflowfunc, pointcoordsbyID, linestringsWithCoords)
  return (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, hExitNodes, argsflowfunc, pointcoordsandindexbyID, pointNodeIDbyIndex, linestringsWithCoords)


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

#def densityasfunctionofheadloss( hval, i, j, argsflowfunc ):
#  return hval
##  L = argsflowfunc['linklength'][ frozenset({i,j}) ]
##  return hval / L  # try, instead of density = headloss, density = headloss/L; this doesn't appear to help in finding valid solutions of SEM2 (or finding invalid ones more quickly), but (TODO) could try an initial h-vector that isn't all zeros.

def triangularflowfunc(i, j, hiMinushj, argsflowfunc):
# Calculate flow across all lanes.
#  hval = abs(hiMinushj)
##  densityval = densityasfunctionofheadloss( hval, i, j, argsflowfunc )
#  densityval = hval
  densityval = abs(hiMinushj)
  numlanes = argsflowfunc['numlanes']
  numlanes_current = numlanes[frozenset({i,j})]
  maximumdensityperlane = argsflowfunc['maximumdensityperlane']  # a constant over all edges
#  maximumdensityperlane_current = maximumdensityperlane[frozenset({i,j})]
#  if hval > numlanes_current*maximumdensityperlane:
  if densityval > numlanes_current*maximumdensityperlane:
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
##  if hval <= numlanes_current * rhoC_current:
#  if hval <= tau_current:
  if densityval <= tau_current:
##    return hsign*lK_current*hval*numlanes_current  # maximum |q| is numlanes_current*lK/(lK*T+1/maximumdensityperlane)  # from p. 93 of Treiber and Kesting
#    return hsign*lK_current*hval  # maximum |q| is numlanes_current*lK/(lK*T+1/maximumdensityperlane)  # from p. 93 of Treiber and Kesting
    return hsign*lK_current*densityval  # maximum |q| is numlanes_current*lK/(lK*T+1/maximumdensityperlane)  # from p. 93 of Treiber and Kesting
#  else:  # hval > numlanes_current * rhoC_current
#  else:  # hval > tau_current
  else:  # densityval > tau_current
####    return numlanes_current/T*(1-hsign*hval/maximumdensityperlane)
###    return (1/T)*(1-hsign*hval/(numlanes_current*maximumdensityperlane))
##    return (1/T_current)*(1-hsign*hval/(numlanes_current*maximumdensityperlane))
#    return (1/T_current)*(1-hsign*densityval/(numlanes_current*maximumdensityperlane))
    return (numlanes_current/T_current)*(1-hsign*densityval/(numlanes_current*maximumdensityperlane))

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

##def f(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
#def f(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
def f(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  debug("f(): h is %s, type(h) is %s" % (h, type(h)))
  fvect = []
##  print(f"injectionnodeindices is {injectionnodeindices}")
##  print(f"neighboursThatCanFlowIntoj is {neighboursThatCanFlowIntoj}")
#  for j in injectionnodeindices:  # construct Sum_i q_{ij} = b_j for injection-node indices j
#    totalinflowj = 0.0
#    hj = h[indexamonginjectionorexitnodes[j]]
#    for i in neighboursThatCanFlowIntoj[j]:
#      if i in exitnodes:
#        hi = hExitNodes[indexamonginjectionorexitnodes[i]]
#      else:
##        print("h is %s" % h)
#        hi = h[indexamonginjectionorexitnodes[i]]
##      debug("i %d, j %d, hi %.10f, hj %.10f" % (i, j, hi, hj))
#      debug("i %d, j %d" % (i, j))
#      debug("hi %.10f, hj %.10f" % (hi, hj))
##      qij = flowijfunc[flowfunc](i, j, hi, hj, argsflowfunc)
#      qij = flowijfunc[flowfunc](i, j, hi-hj, argsflowfunc)
#      debug("qij is %.5f" % qij)
#      print(f"i {i}, j {j}, hi {hi:.6f}, hj {hj:.6f}, qij {qij:.6f}")
#      totalinflowj += qij
#    totalinflowj += b0[indexamonginjectionorexitnodes[j]]
#    fvect.append( totalinflowj )
##    print("fvect[%d] is" % j, fvect[j])
#    print("f(): totalinflowj at (injection-)node %d is %.10f" % (j, totalinflowj))
  for i in injectionnodeindices:  # construct Sum_i q_{ij} = b_j for injection-node indices i
    totaloutflowi = 0.0
    hi = h[indexamonginjectionorexitnodes[i]]
    for j in neighboursThatCanReceiveFlowFromi[i]:
      if j in exitnodes:
        hj = hExitNodes[indexamonginjectionorexitnodes[j]]
      else:
#        print("h is %s" % h)
        hj = h[indexamonginjectionorexitnodes[j]]
#      debug("j %d, i %d, hi %.10f, hj %.10f" % (j, i, hi, hj))
      debug("j %d, i %d" % (j, i))
      debug("hi %.10f, hj %.10f" % (hi, hj))
      qij = flowijfunc[flowfunc](i, j, hi-hj, argsflowfunc)
      debug("qij is %.5f" % qij)
#      print(f"i {i}, j {j}, hi {hi:.6f}, hj {hj:.6f}, qij {qij:.6f}")
      totaloutflowi += qij
    totaloutflowi -= b0[indexamonginjectionorexitnodes[i]]
    fvect.append( totaloutflowi )
#    print("fvect[%d] is" % i, fvect[i])
    debug("f(): totaloutflowi at (injection-)node %d is %.10f" % (i, totaloutflowi))
  debug("fvect is %s" % fvect)
  return np.array(fvect)

#def fnormsquared(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
def fnormsquared(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
#  normsquared = scipy.linalg.norm( f(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )**2
  normsquared = scipy.linalg.norm( f(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )**2
#  debug("normsquared is %.10f" % normsquared)
  return normsquared

#def fnorm(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
def fnorm(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
#  norm = scipy.linalg.norm( f(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )
  norm = scipy.linalg.norm( f(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )
#  debug("norm is %.10f" % norm)
  return norm

#  congestiononlink = linkcongestion(r, s, qrs, abs(hr-hs), argsflowfunc)
#  congestiononlink = linkcongestion(p, j, flowenteringparent, densityinlink, argsflowfunc)
def linkcongestion(r, s, flow, density, argsflowfunc):
  # The current measure of congestion, which is zero if the link's density is <= density at capacity and is the link's capacity less its actual flow otherwise, is problematic; for example, a link into which a congestion-zone has flowed upsttream might have very low flow only because its head-node has very low demand (i.e. injection-flow), and this will be measured to be high congestion. A better measure might follow from assessing congestion not on the links but at the injection-nodes, taking it to be the demand (injection-flow) less the total inflow into links of the network (noting that total inflow will be less than demand where one or more of the adjacent links has had its flow reduced by an upstream-propagating congestion-zone).
#  print(f"linkcongestion(): calculating congestion on link {frozenset({r,s})} from flow {flow} and density {density}...")
  humpdensity = argsflowfunc['tau'][frozenset({r,s})]
##  print("tau[%s] is %.4f" % ((r,s), argsflowfunc['tau'][(r,s)]))
#  print("tau[%s] is %.4f" % ((r,s), argsflowfunc['tau'][frozenset({r,s})]))
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

def congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
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

def totalcongestioninnetwork(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  (congestionij, totalcongestion) = congestioninnetwork(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
  return totalcongestion

#def fnormlesstotalcongestion(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
def fnormlesstotalcongestion(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  totalcongestionfactor = .2  # TODO: what should this value be? - consider using Lagrange multipliers
#  return fnorm(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) - totalcongestionfactor*totalcongestioninnetwork(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
  return fnorm(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) - totalcongestionfactor*totalcongestioninnetwork(h, injectionnodeindices, indexamonginjectionorexitnodes, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)

#def fnormplustotalcongestion(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
def fnormplustotalcongestion(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc):
  totalcongestionfactor = .2  # TODO: what should this value be? - consider using Lagrange multipliers
#  return fnorm(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) + totalcongestionfactor*totalcongestioninnetwork(h, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
  return fnorm(h, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) + totalcongestionfactor*totalcongestioninnetwork(h, injectionnodeindices, indexamonginjectionorexitnodes, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)

eps_GlobalMinconvergence = 7e-8  # empirically obtained value that seems (so far) to yield solutions comparable to those given by Newton-Raphson with eps_NRconvergence = 1e-3 (James Hilton's value)
def printlocalminimum(h, fAtMinimum, accepted):
#  print("Local minimum %.4e at h %s: accepted is %s" % (fAtMinimum, h, accepted))
  # Stop the basinhopping routine if the function-value at this minimum is smaller than a threshold:
  if fAtMinimum < eps_GlobalMinconvergence and accepted: return True

#def solveforhvector(NEWTONRAPHSON, N, injectionnodes, exitnodes, indexamonginjectionorexitnodes, terminatorPressure, hExitNodes, neighboursThatCanFlowIntoj, flowfunction, derivfunction, hJmin, hJmax, hinitInjectionNodes, argsflowfunc, ENERGYFUNCTION, flowijfunc, flowfunc, b0):
def solveforhvector(NEWTONRAPHSON, N, injectionnodeindices, exitnodeindices, indexamonginjectionorexitnodes, terminatorPressure, hExitNodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, flowfunction, derivfunction, hJmin, hJmax, hinitInjectionNodes, argsflowfunc, ENERGYFUNCTION, flowijfunc, flowfunc, b0):
  if not NEWTONRAPHSON:  # use a minimum-finding method to find the minimum of ||f(h)||^2; the chosen method seeks either a local or a global minimum, and either uses the Jacobian or is a Jacobian-less method such as the secant method or a hybrid method, perhaps approximating the Jacobian:
#    hInjectionNodes = scipy.optimize.newton(f, np.zeros(len(injectionnodes)), tol=1e-9, maxiter=100)  # uses the secant method if "derivative" 'fprime' of f isn't given, otherwise uses the Newton-Raphson method; however, it's not clear from the documentation of 'newton' that 'fprime' is a Jacobian or that this method can handle multivariate problems. (Untried: if the "second-order derivative" 'fprime2' of f is also given, then uses Halley’s method.) Using 'tol' and 'maxiter', as for example in 'tol=1e-9, maxiter=100', allows a more precise solution to be found, at the cost of running more than the default number (fifty) of iterations.
#    hInjectionNodes = scipy.optimize.fsolve(f, np.zeros(len(injectionnodes)))  # when a function 'fprime' to compute the Jacobian of f isn't specified, the Jacobian is calculated by a forward-difference approximation ('fsolve' is a wrapper around MINPACK’s 'hybrd' algorithm, which finds a zero of a system of N non-linear functions in N variables by a modification of the Powell hybrid method). This 'fsolve' method requires fewer iterations than 'newton' when neither is supplied with an 'fprime' function for the Jacobian/"derivative".

##    hInjectionNodes = scipy.optimize.fmin_bfgs(fnormsquared, np.zeros(len(injectionnodes)))  # getting trapped in local minimum?
#    res = scipy.optimize.minimize(fnormsquared, np.zeros(len(injectionnodes)), tol=1.0)

##    hJ0 = random.choices(range(hJmin, hJmax+1), k=3)  # initial guess at solution
#    hJ0 = [0] * len(injectionnodes)  # default (zero-vector) initial guess at solution
    hJ0 = [0] * len(injectionnodeindices)  # default (zero-vector) initial guess at solution

    hJ0 = hinitInjectionNodes
    print("Initial guess for global optimiser is h = %s" % repr(hJ0))

#    # Specify bounds on hJ (hInjectionNodes), as required by L-BFGS-B:
    # Specify bounds on hJ (hInjectionNodeIndices), as required by L-BFGS-B:
#    hJmin = int(-4)
#    hJmax = int(4)

#    hJmin = int(-5)  # values used to obtain problematic solution in which h_1 is negative despite h_2 and b_2 being positive
#    hJmax = int(4)
#    hJmin = int(-18)  # values used to obtain a second problematic solution in which h_1 is large and negative; implies h_1 could be arbitrarily large and negative and still yield a solution
#    hJmax = int(5)

##    hJminvals = [hJmin] * 2
#    hJminvals = [hJmin] * len(injectionnodes)
    hJminvals = [hJmin] * len(injectionnodeindices)
#    hJmaxvals = [hJmax] * len(injectionnodes)
    hJmaxvals = [hJmax] * len(injectionnodeindices)
    bounds = [(low, high) for low, high in zip(hJminvals, hJmaxvals)]
    # Method L-BFGS-B is for bound-constrained minimisation; works when (check this) the problem is smooth and bounded (is it?):
#    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, tol=.001)
#    minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds, tol=.0001)
#    minimizer_kwargs = dict(method="L-BFGS-B", tol=.0001)  # drop the bounds-constraint: seems to make the global optimiser run in less time.
#    minimizer_kwargs = dict(args=(injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0), method="L-BFGS-B", tol=.0001)
    if ENERGYFUNCTION == 'conservationOfMassResidual':
      tol=.0001
    elif ENERGYFUNCTION == 'conservationOfMassResidualPlusTotalCongestion':
      tol=.01
    elif ENERGYFUNCTION == 'conservationOfMassResidualLessTotalCongestion':
#      tol=.01  # TODO: establish a good value
      pass
##    minimizer_kwargs = dict(args=(injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc), method="L-BFGS-B", tol=tol)  # 'gtol' option of L-BFGS-B is set to value of 'tol'
    minimizer_kwargs = dict(args=(injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc), method="L-BFGS-B", tol=tol)  # 'gtol' option of L-BFGS-B is set to value of 'tol'
#    minimizer_kwargs = dict(args=(injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc), method="L-BFGS-B", tol=tol)  # 'gtol' option of L-BFGS-B is set to value of 'tol'
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

#    hInjectionNodes = res.x  # the solution-array returned
    hInjectionNodeIndices = res.x  # the solution-array returned
    print("res is %s" % res)
#    print(f"hInjectionNodeIndices is {hInjectionNodeIndices}")


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
    h[exitnodeindices] = terminatorPressure
#    hInjectionNodes = terminatorPressure*np.ones( len(injectionnodes) )
#    h[injectionnodes] = 0.0  # shouldn't need as h is set to np.zeros(N) above
###    h[1] = 2.5  # specify a non-zero initial guess for h
###    h[2] = 0.0
##    h[1], h[2], h[3] = (0, 1.0, 4.0)
##    h[1], h[2] = (1.0, 4.0)
#    for ind, jn in enumerate(injectionnodes):
    for ind, jn in enumerate(injectionnodeindices):
      h[jn] = hinitInjectionNodes[ind]
    print("Initial guess for Newton-Raphson is h = %s" % h)
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
        if k in exitnodeindices:
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
        if i in exitnodeindices:
          hi = hExitNodes[indexamonginjectionorexitnodes[i]]
        else:
#          hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
          hi = h[i]
#        for j in [k for k in range(i+1,N) if k in neighboursThatCanFlowIntoj[i]]:
#        for j in neighboursThatCanFlowIntoj[i]:
        for j in neighboursThatCanReceiveFlowFromi[i]:
          if j in exitnodeindices:
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
#          if i in injectionnodes and j in injectionnodes:
          if i in injectionnodeindices and j in injectionnodeindices:  # we set A[i,j] and A[j,i] only for an edge {i,j} with non-leaf endpoint-nodes. On the other hand, for an edge {i,j} adjacent to a exit-node we don't set A[i,j] or A[j,i], the (i,j)th and (j,i)th elements of J_F(h); this is because, if w.l.o.g. i is the exit-node, then F_i = h_i and so A[i,j] = \partial F_i / \partial h_j = 0, and in the system of equations -A*dh = F(h), A[j,i] will be multiplied by dh(i), which is zero (as h(i) is fixed at 'terminatorPressure'), meaning that A[j,i] might as well be zero.
            A[i, j] = -dval
            A[j, i] = -dval
            debug("A[%d,%d] and A[%d,%d] set to -dval of %.5f" % (i, j, j, i, -dval))
#          if i in injectionnodes:
          if i in injectionnodeindices:
            d[i] += dval
            debug("d[%d] += dval of %.5f" % (i, dval))
            b[i] -= qval
            debug("b[%d] -= qval of %.5f" % (i, qval))
#          if j in injectionnodes:
          if j in injectionnodeindices:
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
#    hInjectionNodes = h[injectionnodes]
    hInjectionNodeIndices = h[injectionnodeindices]
#    print("fnormsquared(hInjectionNodes) is %.6e" % fnormsquared(hInjectionNodes))
####    print("fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0) is %.6e" % fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, hvalmax, lK, b0))
###    print("fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) is %.6e" % fnormsquared(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
##    print("fnorm(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) is %.6e" % fnorm(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
#    print("fnorm(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) is %.6e" % fnorm(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
    print("fnorm(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) is %.6e" % fnorm(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
    normsratio = dh_dot/dh0_dot
#    if dh_dot/dh0_dot > eps_NRconvergence and iteration >= MAXITERS:
    if normsratio > eps_NRconvergence and iteration >= MAXITERS:
      print("=========== ALERT: Newton-Raphson method has not converged: has reached maximum number (%d) of iterations without ||dh||^2 (%.12f) falling to %.4e * ||dh0||^2 (||dh0||^2 = %.5f); ||dh||^2/||dh0||^2 = %.6e" % (MAXITERS, dh_dot, eps_NRconvergence, dh0_dot, normsratio))
#    if dh_dot/dh0_dot <= eps_NRconvergence:
    if normsratio <= eps_NRconvergence:
      print("Newton-Raphson method has converged, in that ||dh||^2 (%.12f) has fallen to no more than %.4e * ||dh0||^2 (||dh0||^2 = %.5f); ||dh||^2/||dh0||^2 = %.6e" % (dh_dot, eps_NRconvergence, dh0_dot, normsratio))


#  print("Solution for injection-nodes is hInjectionNodes =", hInjectionNodes)
  print("Solution for injectionnodeindices is hInjectionNodeIndices =", hInjectionNodeIndices)
###  print("Inflows at nodes for this solution is %s" % f(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
##  print("Inflows at nodes for this solution is %s" % f(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
#  print("Inflows at nodes for this solution is %s" % f(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
  print("Outflows at nodes for this solution is %s" % f(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc))
  if NEWTONRAPHSON:
    print("Flow network took %d Newton-Raphson iterations." % iteration)
#  print("Function f was called %d times." % numcallsoffuncf)
##  print("The solution-finding method took %.5f seconds." % (endtime-starttime))

##   Calculate Qmax (maximum flow) and totaloutflow
#   Calculate Qabsmax (maximum absolute flow) and totaloutflow
  Qabsmax = 0.0
  qij = {}  # store flows on all edges
#  congestionij = {}  # store congestion-values on all edges
#  for i in range(N-1):
  for i in range(N):
    if i in exitnodeindices:
      hi = hExitNodes[indexamonginjectionorexitnodes[i]]
    else:
#      hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
      hi = hInjectionNodeIndices[indexamonginjectionorexitnodes[i]]
#    for j in [k for k in range(i+1,N) if k in neighboursThatCanFlowIntoj[i]]:
    for j in neighboursThatCanReceiveFlowFromi[i]:
      if j in exitnodeindices:
        hj = hExitNodes[indexamonginjectionorexitnodes[j]]
      else:
#        hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
        hj = hInjectionNodeIndices[indexamonginjectionorexitnodes[j]]
#      hval = abs(hi-hj)
      print("i %d, j %d, hi %.3f, hj %.3f, hval %.3f" % (i, j, hi, hj, abs(hi-hj)))
#      qij[i, j] = flowijfunc[flowfunc](i, j, hi, hj, argsflowfunc)
      qij[i, j] = flowijfunc[flowfunc](i, j, hi-hj, argsflowfunc)
      print("i %d, j %d, hi %.5f, hj %.5f, qij[%d,%d] is %.5f" % (i, j, hi, hj, i, j, qij[i,j]))
      Qabsmax = max(Qabsmax, abs(qij[i,j]))
#  for j in range(N):
#    if j in exitnodeindices:
#      hi = hExitNodes[indexamonginjectionorexitnodes[j]]
#    else:
##      hi = hInjectionNodes[indexamonginjectionorexitnodes[j]]
#      hi = hInjectionNodeIndices[indexamonginjectionorexitnodes[j]]
##    for i in [k for k in range(j+1,N) if k in neighboursThatCanFlowIntoj[j]]:
#    for i in neighboursThatCanFlowIntoj[j]:
#      if i in exitnodeindices:
#        hj = hExitNodes[indexamonginjectionorexitnodes[i]]
#      else:
##        hj = hInjectionNodes[indexamonginjectionorexitnodes[i]]
#        hj = hInjectionNodeIndices[indexamonginjectionorexitnodes[i]]
##      hval = abs(hi-hj)
#      print("i %d, j %d, hi %.3f, hj %.3f, hval %.3f" % (i, j, hi, hj, abs(hi-hj)))
##      qij[j, i] = flowijfunc[flowfunc](j, i, hi, hj, argsflowfunc)
#      qij[i, j] = flowijfunc[flowfunc](i, j, hi-hj, argsflowfunc)
#      print("i %d, j %d, hi %.5f, hj %.5f, qij[%d,%d] is %.5f" % (i, j, hi, hj, i, j, qij[i,j]))
#      Qabsmax = max(Qabsmax, abs(qij[i,j]))
  print("Qabsmax is", Qabsmax)

  outflowsByNodeID = {}  # outflow at an exit-node is defined to be the negative of its total inflow
  totaloutflow = 0.0
  for t in exitnodeindices:
    outflowsByNodeID[t] = 0.0
#    totaloutflowj = 0.0
#    ht = h[indexamonginjectionorexitnodes[t]]
    ht = hExitNodes[indexamonginjectionorexitnodes[t]]
#    print("Node %d has neighbours %s" % (t, neighboursThatCanFlowIntoj))
#    assert len(neighboursThatCanFlowIntoj[t]) == 1  # no longer valid, as can now choose a non-leaf node to be a exit- (i.e. exit-) node
    for i in neighboursThatCanFlowIntoj[t]:
##      assert i in injectionnodes  # this is invalid now that any node can be an exit-node: can have adjacent exit-nodes
###      if i in exitnodes:
####        hi = hExitNodes[indexamonginjectionorexitnodes[i]]
###       pass
###      else:
###        hi = h[indexamonginjectionorexitnodes[i]]
##      hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
#      if i in injectionnodes:
      if i in injectionnodeindices:
#        hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
        hi = hInjectionNodeIndices[indexamonginjectionorexitnodes[i]]
      else:
        hi = hExitNodes[indexamonginjectionorexitnodes[i]]
      hval = abs(hi-ht)
      print(f"qij is {qij}")
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
#  (congestionij, totalnetworkcongestion) = congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
  (congestionij, totalnetworkcongestion) = congestioninnetwork(hInjectionNodeIndices, injectionnodeindices, indexamonginjectionorexitnodes, exitnodeindices, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
#  return (exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, iteration, outflowsByNodeID)
  return (exception, hInjectionNodeIndices, qij, flowbalanced, congestionij, totalnetworkcongestion, iteration, outflowsByNodeID)




def generateAllIndices( n, numpointsperdimension ):  # return a list of n-tuples, each tuple containing n integers between zero and (numpointsperdimension-1) inclusive, with each possible such tuple included in the list exactly once
  if n == 1:
#    return [(val,) for val in range(maxindex+1)]
    return [(val,) for val in range(numpointsperdimension)]
  elif n > 1:
#    return [ (val,) + t for val in range(maxindex+1) for t in generateAllIndices(n-1, maxindex) ]
    return [ (val,) + t for val in range(numpointsperdimension) for t in generateAllIndices(n-1, numpointsperdimension) ]



def FloydWarshall( arcs, vertices, arcweights ):
  """
  A dynamic-programming algorithm: find shortest path between every pair of nodes, given the distance between each pair of arc-adjacent nodes as stored in 'arcweights'. In the returned arrays, pathdistance[i,j] is the shortest-path distance from i to j, and nextNodeInShortestPath[i,j] is the next node after i in the shortest path from i to j.
  """
  N = len(vertices)
#  pathdistance = 9e99 * np.ones( shape=(N,N) )  # set each path-distance to "infinity"
  pathdistance = np.full( shape=(N,N), fill_value=9e99  )  # set each path-distance to "infinity"
  nextNodeInShortestPath = np.full( shape=(N,N), fill_value=None )
  for arc in arcs:
#    arcAsList = list(arc)
##    print("arcAsList is", arcAsList)
#    (a, b) = arcAsList
    (a, b) = arc
#    pathdistance[a, b] = edgeweights[e]
    pathdistance[a, b] = arcweights[a, b]
    nextNodeInShortestPath[a, b] = b
  for v in vertices:
    pathdistance[v, v] = 0
    nextNodeInShortestPath[v, v] = v
#  print("FloydWarshall(): pathdistance is", pathdistance)
#  print("FloydWarshall(): nextNodeInShortestPath is", nextNodeInShortestPath)
  for k in vertices:
    for i in vertices:
      for j in vertices:
        if pathdistance[i, j] > pathdistance[i, k] + pathdistance[k, j]:
          pathdistance[i, j] = pathdistance[i, k] + pathdistance[k, j]
          nextNodeInShortestPath[i, j] = nextNodeInShortestPath[i, k]
#          print("k %d: pathdistance[%d, %d] is %s" % (k, i, j, pathdistance[i,j]))
#          print("nextNodeInShortestPath[%d, %d] is %s" % (i, j, nextNodeInShortestPath[i,j]))
  return (pathdistance, nextNodeInShortestPath)


def constructShortestPath(i, j, nextNodeInShortestPath, edgeweights):
  """
  Given the distance between each pair of adjacent vertices as stored in 'edgeweights', and 'nextNodeInShortestPath' in which nextNodeInShortestPath[i,j] is the next vertex after i in the shortest path from i to j, construct the shortest path from i to j.
  """
  print(f"Entering constructShortestPath({i}, {j}, nextNodeInShortestPath, edgeweights)")
  if nextNodeInShortestPath[i, j] == None:
##    return []
#    return ([], f"infinity: can't go from node {i} to {j}") 
    return ([], 9e99) 
  path = [i]
#  while i != j:
#    i = nextNodeInShortestPath[i, j]
#    path.append( i )
  pathlength = 0.0
  a = i
  while a != j:
    b = nextNodeInShortestPath[a, j]
    pathlength += edgeweights[a, b]
    path.append( b )
    a = b
  return (path, pathlength)


#def adjustSpeedandDensityToLink( flowfunction, p, j, densityenteringparent, argsflowfunc, flowenteringparent, speedenteringparent ):
#  EPS = 1e-3
##  EPS = 1e-0  # might be all the precision that's required
#  densityinlink = densityenteringparent
#  previousspeedinlink = 1e9
#  speedinlink = speedenteringparent
#  numiterations = 0
##  if flowfunc in ('concaveQuadratic', 'xExpnegx', 'triangular'):
##    humpdensity = argsflowfunc['tau'][frozenset({p,j})]
##  elif flowfunc == 'triangular':
##    humpdensity = numlanes*argsflowfunc['rhoC'][frozenset({p,j})]
#  humpdensity = argsflowfunc['tau'][frozenset({p,j})]
##  if densityinlink >= humpdensity:
#  if densityinlink > humpdensity:
#    print(f"WARNING: density in link {{{p},{j}}} of {densityinlink:.6f} is over the hump {humpdensity:.6f} by a gap of {(densityinlink-humpdensity):.6e}.")
#  # Iterate until speedinlink converges:
#  CALCULATESPEEDBEFOREDENSITY = True
##  CALCULATESPEEDBEFOREDENSITY = False
#  print("CALCULATESPEEDBEFOREDENSITY is", CALCULATESPEEDBEFOREDENSITY)
#  while abs(speedinlink - previousspeedinlink) > EPS:
#    previousspeedinlink = speedinlink
#    # First, as speed is a function of density, the speed adjusts in accordance with the density and the conditions (free speed, number of lanes) on the link
#    theoreticalflowinlink = flowfunction(p, j, densityinlink, argsflowfunc)  # flow that, for the given density, would result from the conditions (free speed, number of lanes) on this link
#    print(f"iteration {numiterations}: given density of {densityinlink:.3f}, theoreticalflowinlink is {theoreticalflowinlink:.3f}")
#    if CALCULATESPEEDBEFOREDENSITY:
##      if abs(densityinlink) >= EPS:
#      if abs(densityinlink) >= EPS or abs(theoreticalflowinlink) >= EPS:
#        speedinlink = theoreticalflowinlink / densityinlink  # v = q/rho: the speed adjusts in accordance with the density and the conditions (free speed, number of lanes) on the link
##        print(f"adjustSpeedandDensityToLink(): speedinlink adjusted to {speedinlink:.3E}")
#        print(f"speedinlink adjusted to {theoreticalflowinlink:.3f}/{densityinlink:.3f} = {speedinlink:.3f}")
#      else:  # both flow and density in link are close to zero: speedinlink retains its existing value
#        pass
##      assert speedinlink >= 0.0
#      assert speedinlink >= -EPS, f"speedinlink is {speedinlink:.3f}, numiterations is {numiterations}"
#      # Second, the density adjusts in accordance with the original entering flow and the new speed on the link (density = flow/speed)
#      if abs(speedinlink) >= 1e-18:
#        densityinlink = flowenteringparent / speedinlink  # the density adjusts in accordance with the new adjusted speed and the constraint that flow is conserved
#        print(f"densityinlink adjusted to {flowenteringparent:.3f}/{speedinlink:.3f} = {densityinlink:.3f}")
#      else:
#        densityinlink = 9e99  # TODO: is it valid to set this to a high value, or should it be set to infinity/inf/NAN?
##        densityinlink = 1e18
##        densityinlink = 9.99999999999e17
#        print(f"densityinlink adjusted to {densityinlink:.3f}")
#      print(f"densityinlink is {densityinlink:.3f}")
#    else:  # CALCULATESPEEDBEFOREDENSITY == False
#      if abs(densityinlink) >= EPS or abs(theoreticalflowinlink) >= EPS:
#        densityinlink = theoreticalflowinlink / speedinlink  # rho = q/v: the density adjusts in accordance with the speed and the conditions (free speed, number of lanes) on the link
#        print(f"densityinlink adjusted to {theoreticalflowinlink:.3f}/{speedinlink:.3f} = {densityinlink:.3f}")
#      else:  # both flow and density in link are close to zero: densityinlink retains its existing value
#        pass
#      assert densityinlink >= -EPS, f"densityinlink is {densityinlink:.3f}, numiterations is {numiterations}"
#      # Second, the speed adjusts in accordance with the original entering flow and the new density on the link (speed = flow/density)
#      if abs(densityinlink) >= 1e-18:
#        speedinlink = flowenteringparent / densityinlink  # the speed adjusts in accordance with the new adjusted density and the constraint that flow is conserved
#        print(f"speedinlink adjusted to {flowenteringparent:.3f}/{densityinlink:.3f} = {speedinlink:.3f}")
#      else:
#        speedinlink = 9e99  # TODO: is it valid to set this to a high value, or should it be set to infinity/inf/NAN?
##        speedinlink = 1e18
#        print(f"speedinlink adjusted to {speedinlink:.3f}")
#      print(f"speedinlink is {speedinlink:.3f}")
#    numiterations += 1
#  print(f"abs(theoreticalflowinlink - flowenteringparent) is {abs(theoreticalflowinlink - flowenteringparent):.3E}")
#  EPS_FLOWINLINK = 2e-2
#  assert abs(theoreticalflowinlink - flowenteringparent) <= EPS_FLOWINLINK
#  print(f"adjustSpeedandDensityToLink(): numiterations is {numiterations}")
#  print(f"adjustSpeedandDensityToLink(): {densityinlink:.6f}*{speedinlink:.6f} = {densityinlink*speedinlink:.6f}")
#  return (speedinlink, densityinlink, numiterations)


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


def propagateSEM3congestionUpstreamToParents( i, totalflowleavingnodeiDownstreamOfCongestion, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes ):
  print(f"totalflowleavingnodeiDownstreamOfCongestion at node {i} is {totalflowleavingnodeiDownstreamOfCongestion}")
  initialtotalflowEnteringNodeiFromItsParents = sum([qijatnodej[p, i] for p in parentsOnShortestPathsFromInjectionNodes[i]])
  print(f"Initial total flow entering node {i} from its parent-nodes is {initialtotalflowEnteringNodeiFromItsParents:.3f}")
  totalflowEnteringNodei = 0.0
  if i in injectionnodes:
    injectionflowati = b0[indexamonginjectionorexitnodes[i]]
    if injectionflowati > 0.0:
      print(f"Adding to totalflowEnteringNodei the injectionflow at node {i} of {injectionflowati:.5f}")
    totalflowEnteringNodei += injectionflowati
  for p in parentsOnShortestPathsFromInjectionNodes[i]:
##    qijatnodej[p, i] *= totalflowleavingnodeiDownstreamOfCongestion / initialtotalflowEnteringNodei  # divide congested flow between the links entering node i: the flow in each decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
#    qijatnodej[p, i] *= (totalflowleavingnodeiDownstreamOfCongestion - injectionflowati) / initialtotalflowEnteringNodeiFromItsParents  # divide congested flow between the links entering node i: the flow in each decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
#    congestedflowpi = qijatnodej[p, i]
    initialflowpi = qijatnodej[p, i]
    print(f"initial flow[{p}, {i}] is {initialflowpi}")
    if initialflowpi == 0.0:  # no need to consider a parent-link that has no incoming flow, as congestion can't affect it  # FIXME: should check for closeness to zero of flow's absolute value
      print(f"initialflowpi of {initialflowpi} is zero, which cannot be reduced by congestion: there is no real congestion in link ({p}, {i}).")
      continue
#    congestedflowpi = qijatnodej[p, i] * (totalflowleavingnodeiDownstreamOfCongestion - injectionflowati) / initialtotalflowEnteringNodeiFromItsParents  # divide congested flow between the links entering node i: the flow in each decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
    congestedflowpi = qijatnodej[p, i] * totalflowleavingnodeiDownstreamOfCongestion / (initialtotalflowEnteringNodeiFromItsParents + injectionflowati)  # divide congested flow between the links entering node i: the flow in each link, and also by implication the injection-flow, decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
    print(f"congested flow[{p}, {i}] is {congestedflowpi}")
    print(f"As congestion-zone propagates upstream, flow in link [{p}, {i}] decreases from {initialflowpi} to {congestedflowpi}")
    qijatnodej[p, i] = congestedflowpi
    totalflowEnteringNodei += congestedflowpi
    edge = frozenset({p, i})
#    densityijatnodej[frozenset({p, i})] = densityoncongestedbranchoftriangularfundamentaldiagram( qijatnodej[p, i], edge, argsflowfunc )
    congesteddensitypi = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowpi, edge, argsflowfunc )
    densityijatnodej[edge] = congesteddensitypi
    print(f"densityijatnodej is now {densityijatnodej}")
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
#    speedijatnodej[frozenset({p, i})] = qijatnodej[p, i] / densityijatnodej[edge]
    speedijatnodej[edge] = congestedflowpi / congesteddensitypi
    print(f"Density in link {edge} is now {densityijatnodej[edge]}, speed is now {speedijatnodej[edge]}.")
    propagateSEM3congestionUpstreamToParents( p, congestedflowpi, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )

  if totalflowleavingnodeiDownstreamOfCongestion < totalflowEnteringNodei:
#    print(f"ALERT: totalflowEnteringNodei of {totalflowEnteringNodei} not equal to totalflowleavingnodeiDownstreamOfCongestion of {totalflowleavingnodeiDownstreamOfCongestion}.")
    print(f"ALERT: at node {i}, totalflowleavingnodeiDownstreamOfCongestion {totalflowleavingnodeiDownstreamOfCongestion} < totalflowEnteringNodei {totalflowEnteringNodei}: node {i} must be an injection-node, and there is a shock where traffic-flow enters the network there.")
    assert i in injectionnodes
#  assert totalflowEnteringNodei == totalflowleavingnodeiDownstreamOfCongestion
  assert totalflowleavingnodeiDownstreamOfCongestion <= totalflowEnteringNodei


def propagateSEM4congestionUpstreamToParents( upstreamPropagatingShockFront, a, timeThatShockFrontReachedNodea, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, densityijatnodej, speedijatnodej, flowstateijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes, propagatingwavefronts, nexttimestepduration ):
  EPS = 1e-3
  congestionAlertMessage = None  # will change if congestion is detected at an injection-node
  congestedFlowThatHasPropagatedUpstreamToNodea = upstreamPropagatingShockFront['flow']
  print(f"congestedFlowThatHasPropagatedUpstreamToNodea at node {a} is {congestedFlowThatHasPropagatedUpstreamToNodea} veh/hr.")
  initialtotalflowEnteringNodeaFromItsParents = sum([qijatnodej[p,a] for p in parentsOnShortestPathsFromInjectionNodes[a]])
  print(f"Initial total flow entering node {a} from its parent-nodes is {initialtotalflowEnteringNodeaFromItsParents:.3f} veh/hr.")
  newtotalflowEnteringNodea = 0.0  # will include the injection-flow at a, if any
  if a in injectionnodes:
    injectionflowata = b0[indexamonginjectionorexitnodes[a]]
#    if injectionflowata > 0.0:
#      print(f"Adding to newtotalflowEnteringNodea the injectionflow at node {a} of {injectionflowata:.5f}")
#    newtotalflowEnteringNodea += injectionflowata
  else:
    injectionflowata = 0.0
  print(f"injectionflow at node {a} is {injectionflowata:.5f} veh/hr.")
  flowReductionFactor = congestedFlowThatHasPropagatedUpstreamToNodea / (initialtotalflowEnteringNodeaFromItsParents + injectionflowata)  # shall divide congested flow between the links entering node a: the flow in each link, and also by implication the injection-flow, decreases proportionally  # TODO: improve this if not realistic, or find a better method in the literature
  print(f"flowReductionFactor is congestedFlowThatHasPropagatedUpstreamToNodea / (initialtotalflowEnteringNodeaFromItsParents + injectionflowata) = {congestedFlowThatHasPropagatedUpstreamToNodea} / ({initialtotalflowEnteringNodeaFromItsParents} + {injectionflowata}) = {flowReductionFactor}")
  congestedinjectionflowata = flowReductionFactor * injectionflowata
  if congestedinjectionflowata < injectionflowata:
    congestionAlertMessage = f"ALERT: at injection-node {a}, congestedinjectionflowata {congestedinjectionflowata} < injectionflowata {injectionflowata} veh/hr: there is a shock where traffic-flow enters the network there."
    print(congestionAlertMessage)
    assert a in injectionnodes
  print(f"Adding to newtotalflowEnteringNodea the congestedinjectionflowata (possibly zero) at node {a} of {congestedinjectionflowata:.5f} veh/hr")
  newtotalflowEnteringNodea += congestedinjectionflowata
##  assert flowReductionFactor < 1.0  # not generally true: flowReductionFactor is sometimes equal to 1
#  assert flowReductionFactor <= 1.0  # TODO: perhaps is sometimes greater than 1, but could be due to a bug - must be careful not to propagate a shock-front upstream to a "parent"-node on a shortest path coming from an injection-node that didn't original supply the flow trailed by the shock-front
##  if flowReductionFactor < 1.0:  # if flow is reduced by the congestion, then congestion affects the flow on every link entering node a:
  for p in parentsOnShortestPathsFromInjectionNodes[a]:
    initialflowpa = qijatnodej[p,a]
    print(f"initial flow[{p},{a}] is {initialflowpa}")
    if initialflowpa == 0.0:  # no need to consider a parent-link that has no incoming flow, as congestion can't affect it  # FIXME: should check for closeness to zero of flow's absolute value
      print(f"initialflowpa of {initialflowpa} is zero, which cannot be reduced by congestion: there is no real congestion in link ({p},{a}).")
      continue
    congestedflowpa = flowReductionFactor * initialflowpa
    print(f"congested flow[{p},{a}] is {congestedflowpa}")
    if flowReductionFactor < 1.0:
      print(f"As congestion-zone propagates upstream, flow in link [{p},{a}] decreases from {initialflowpa} to {congestedflowpa}")
      qijatnodej[p,a] = congestedflowpa
      print(f"qijatnodej[{p},{a}] is now {qijatnodej[p,a]}")
    print(f"Adding to newtotalflowEnteringNodea the congestedflowpa entering node {a} from its parent {p} of {congestedflowpa:.5f} veh/hr")
    newtotalflowEnteringNodea += congestedflowpa
    edge = frozenset({p,a})
    congesteddensitypa = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowpa, edge, argsflowfunc )
#    assert densityijatnodej[edge] != congesteddensitypa
    densityijatnodej[edge] = congesteddensitypa
    print(f"densityijatnodej[{edge}] is now {densityijatnodej[edge]}")
    upstreamflowpa = initialflowpa  # needed to calculate shock-front's new propagation-speed on link (p,a)
    if flowstateijatnodej[p,a] == 'freeflowstate':
      upstreamdensitypa = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamflowpa, edge, argsflowfunc )
    elif flowstateijatnodej[p,a] == 'congestedflowstate':
      upstreamdensitypa = densityoncongestedbranchoftriangularfundamentaldiagram( upstreamflowpa, edge, argsflowfunc )
    if flowReductionFactor == 1.0:  # TODO: perhaps replace with 'if abs(flowReductionFactor - 1.0) <= EPS:'
      assert abs(upstreamflowpa - congestedflowpa) <= EPS
      print(f"upstreamdensitypa is {upstreamdensitypa}, congesteddensitypa is {congesteddensitypa}.")
      assert abs(upstreamdensitypa - congesteddensitypa) <= EPS
    elif flowReductionFactor < 1.0:
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
#      speedijatnodej[frozenset({p,a})] = qijatnodej[p,a] / densityijatnodej[edge]
      speedijatnodej[edge] = congestedflowpa / congesteddensitypa
      flowstateijatnodej[p,a] = 'congestedflowstate'
      print(f"densityijatnodej{edge} is now {densityijatnodej[edge]}, speedijatnodej{edge} is now {speedijatnodej[edge]}, flowstateijatnodej[{p},{a}] is {flowstateijatnodej[p,a]}.")
      lengthoflinkpa = argsflowfunc['linklength'][ frozenset({p,a}) ]
      timetillupstreamnodeonlink = -lengthoflinkpa / shockFrontPropagationSpeed
##      propagateSEM3congestionUpstreamToParents( p, congestedflowpa, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )
      propagatingwavefronts.append( {'time': timeThatShockFrontReachedNodea, 'currentlink': frozenset({p,a}), 'sourcenodeonlink': a, 'destinationnodeonlink': p, 'distancereacheddownlink': 0.0, 'propagationspeed': -shockFrontPropagationSpeed, 'flow': congestedflowpa, 'density': congesteddensitypa, 'timetilldestinationnodeonlink': lengthoflinkpa/(-shockFrontPropagationSpeed), 'nodesonpathtoanexit': 'upstream-propagating shock-front, so not applicable', 'flowaheadofwavefront': upstreamflowpa, 'densityaheadofwavefront': upstreamdensitypa} )  # 'propagationspeed' of a propagating wavefront is in the direction of its travel, so is always positive

  print(f"upstreamPropagatingShockFront is {upstreamPropagatingShockFront}")
  upstreamPropagatingShockFront['timetilldestinationnodeonlink'] = 0.0
  upstreamPropagatingShockFront['distancereacheddownlink'] += upstreamPropagatingShockFront['propagationspeed'] * nexttimestepduration
  sourceNodeOfShockFront = upstreamPropagatingShockFront['sourcenodeonlink']
  lengthofcurrentlink = argsflowfunc['linklength'][ frozenset({sourceNodeOfShockFront, a}) ]
  try:
    assert abs(upstreamPropagatingShockFront['distancereacheddownlink'] - lengthofcurrentlink) <= EPS
  except:
    print(f"upstreamPropagatingShockFront['distancereacheddownlink'] is {upstreamPropagatingShockFront['distancereacheddownlink']} km, lengthofcurrentlink is {lengthofcurrentlink}")
    raise
  print(f"As it propagates upstream and reaches node {a}, wavefront {upstreamPropagatingShockFront} is removed from propagatingwavefronts.")
  propagatingwavefronts.remove( upstreamPropagatingShockFront )
  print(f"propagatingwavefronts is now {propagatingwavefronts}")

  print(f"newtotalflowEnteringNodea is {newtotalflowEnteringNodea} veh/hr, congestedFlowThatHasPropagatedUpstreamToNodea is {congestedFlowThatHasPropagatedUpstreamToNodea} veh/hr.")
#  if congestedFlowThatHasPropagatedUpstreamToNodea < newtotalflowEnteringNodea:
#    print(f"ALERT: at node {a}, congestedFlowThatHasPropagatedUpstreamToNodea {congestedFlowThatHasPropagatedUpstreamToNodea} < newtotalflowEnteringNodea {newtotalflowEnteringNodea} veh/hr: node {a} must be an injection-node, and there is a shock where traffic-flow enters the network there.")
#    assert a in injectionnodes
  if abs(congestedFlowThatHasPropagatedUpstreamToNodea - newtotalflowEnteringNodea) > EPS:
    print(f"BUG: newtotalflowEnteringNodea of {newtotalflowEnteringNodea} not equal to congestedFlowThatHasPropagatedUpstreamToNodea of {congestedFlowThatHasPropagatedUpstreamToNodea} veh/hr.")
  assert abs(newtotalflowEnteringNodea - congestedFlowThatHasPropagatedUpstreamToNodea) <= EPS
#  assert congestedFlowThatHasPropagatedUpstreamToNodea <= newtotalflowEnteringNodea
  return congestionAlertMessage


#def flowandspeedEnteringNode( i, parentsOnShortestPathsFromInjectionNodes, qijatnodej, congestionij, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnyExit, flowijfunc, flowfunc ):
def flowEnteringNodeSEM3( i, childofi, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, flowijfunc, flowfunc, injectionnodes ):
  # Traffic injected at a node flows towards the next node on its shortest path towards an exit.
#  # At the beginning of the first link the flow and speed are known, which gives the density; for traffic flowing into the link this density gives, according to v = v(rho) for the link, the speed in the link. The flow is conserved between each injection-node and the end of the first link on its shortest path to an exit; if a link is divided for numerical reasons into sub-sections, then the flow is conserved between adjacent sections.
  # At an injection-node the inflow is known; this divides among any adjacent links that are used by shortest paths, and the portion of flow on each such link gives its density, on the free-flow side of its triangular fundamental diagram; the speed is the link's free speed V0 (see Figure 8.9, p. 92 of Treiber and Kesting).
  EPS = 1e-3
  print(f"Entering flowEnteringNodeSEM3( {i}, ... ):")
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
##    flowenteringparent, speedenteringparent = flowandspeedEnteringNode( p, parentsOnShortestPathsFromInjectionNodes, qijatnodej, congestionij, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnyExit, flowijfunc, flowfunc )
#    flowenteringparent, speedenteringparent = flowandspeedEnteringNode( p, i, parentsOnShortestPathsFromInjectionNodes, qijatnodej, congestionij, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, flowijfunc, flowfunc )
    flowenteringparent = flowEnteringNodeSEM3( p, i, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, flowijfunc, flowfunc, injectionnodes )
#    densityenteringparent = flowenteringparent / speedenteringparent
#    print(f"Node {p} has closest exit-node {closestExit[p]}, at distance {shortestDistanceToAnyExit[p]:.4f} with shortest path {shortestPathToAnyExit[p]}")
#    print(f"p is {p}, i is {i}")
    assert i == shortestPathToAnyExit[p][1]  # i is the next node in the shortest path from p to an exit-node
    linkcapacity = argsflowfunc['linkcapacity'][frozenset({p,i})]
    if flowenteringparent > linkcapacity:
      # TODO: is it physically realistic to have traffic-flow in excess of the link-capacity use a different route to an exit, perhaps a second-shortest path? Note there might be several different paths all with the shortest length.
#      print(f"WARNING: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {i}) which has a capacity of {linkcapacity:.3f}; high congestion and density, and low speed, will result.")
      print(f"WARNING: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {i}) which has a capacity of {linkcapacity:.3f}: the link will have flow at this capacity, while a congestion-zone will propagate up the upstream link(s); the congestion-zone will have density to the right of the hump in its link's flow-function, and correspondingly low speed.")
      flowleavingparent = linkcapacity
      densityleavingparent = argsflowfunc['tau'][frozenset({p,i})]
      speedleavingparent = flowleavingparent / densityleavingparent
#      for gp in parentsOnShortestPathsFromInjectionNodes[p]:
      propagateSEM3congestionUpstreamToParents( p, flowleavingparent, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )
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
    qijatnodej[p, i] = flowleavingparent
    densityijatnodej[frozenset({p, i})] = densityinlink
    print(f"AFTER DEFINING densityijatnodej[frozenset({{{p}, {i}}})], densityijatnodej is {densityijatnodej}")
    speedijatnodej[frozenset({p, i})] = speedinlink
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


def examineCandidateWavefrontsForIntersectionWithThisWavefront(thiswf, candidateWavefrontsToIntersect, nextwavefronttopropagate, wavefrontIntersected, L, speed_thiswf, distreacheddownlink_thiswf, minTimeForAWavefrontToIntersectWithAnotherOnSameLink, minTimeForAWavefrontToReachNextNodeOnPath, pointReachedByNextPropagatingWaveFront):  # assume that 'thiswf' is an upstream-propagating wavefront
  for cwf in candidateWavefrontsToIntersect:
    # Calculate remaining time (might be infinity) till candidate wavefront would be intersected by 'thiswf':
    speed_cwf = cwf['propagationspeed']  # propagation-speed of candidate wavefront
    distreacheddownlink_cwf = cwf['distancereacheddownlink']
    print(f"speed_cwf is {speed_cwf} km/hr, distreacheddownlink_cwf is {distreacheddownlink_cwf} km")
    if cwf['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':
      timeTillIntersected = (L - distreacheddownlink_thiswf - distreacheddownlink_cwf) / (speed_thiswf + speed_cwf)
    else:  # the intersecting wavefronts are both propagating upstream
      timeTillIntersected = (distreacheddownlink_cwf - distreacheddownlink_thiswf) / (speed_thiswf - speed_cwf)
    print(f"timeTillIntersected is {timeTillIntersected} hours")
    assert timeTillIntersected > 0.0  # TODO: isn't true if both wavefronts are propagating upstream and cwf is behind thiswf, i.e. hasn't propagated as great a distance up the link
    if timeTillIntersected < minTimeForAWavefrontToIntersectWithAnotherOnSameLink and timeTillIntersected <= minTimeForAWavefrontToReachNextNodeOnPath:  # for the intersection to be such, it must occur by the time (before?) any propagating wavefront reaches its next node (the node at the end of its link)
      assert timeTillIntersected < minTimeForAWavefrontToReachNextNodeOnPath  # TODO: it would be convenient if this were always true, for otherwise nextwavefronttopropagate might reach its destination-node simultaneously with intersecting the other wavefront on the link, although that case might not be too hard to sort out
      minTimeForAWavefrontToIntersectWithAnotherOnSameLink = timeTillIntersected
#      print(f"minTimeForThisWavefrontToIntersectWithAnotherOnSameLink is {minTimeForThisWavefrontToIntersectWithAnotherOnSameLink} hours")
      print(f"minTimeForAWavefrontToIntersectWithAnotherOnSameLink is {minTimeForAWavefrontToIntersectWithAnotherOnSameLink} hours")
      assert cwf['currentlink'] == thiswf['currentlink']
      nextwavefronttopropagate = thiswf
      wavefrontIntersected = cwf
      print(f"wavefrontIntersected is {wavefrontIntersected}")
      print(f"L of {thiswf['currentlink']} is {L} km.")
      if cwf['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':
        pointReachedByNextPropagatingWaveFront = 'intersect with wavefront propagating in opposite direction (downstream) on same link'
      else:  # the intersecting wavefronts are both propagating upstream
        pointReachedByNextPropagatingWaveFront = 'catch up to slower-propagating shockwave on same link'
    elif timeTillIntersected >= minTimeForAWavefrontToIntersectWithAnotherOnSameLink:
      print(f"timeTillIntersected of {timeTillIntersected} is no less than currently-best minTimeForAWavefrontToIntersectWithAnotherOnSameLink of {minTimeForAWavefrontToIntersectWithAnotherOnSameLink} hours, so is rejected.")
    elif timeTillIntersected > minTimeForAWavefrontToReachNextNodeOnPath:
      print(f"timeTillIntersected of {timeTillIntersected} is less than currently-best minTimeForAWavefrontToIntersectWithAnotherOnSameLink of {minTimeForAWavefrontToIntersectWithAnotherOnSameLink} hours, but exceeds minTimeForAWavefrontToReachNextNodeOnPath of {minTimeForAWavefrontToReachNextNodeOnPath} hours, so the intersection cannot occur before the soonest arrival of a wavefront at the end-node of its link.")
  return (minTimeForAWavefrontToIntersectWithAnotherOnSameLink, nextwavefronttopropagate, wavefrontIntersected, pointReachedByNextPropagatingWaveFront)


def calculateNextTimestepDuration(propagatingwavefronts, argsflowfunc):
  LARGEVAL = 9e99
  print(f"propagatingwavefronts is {propagatingwavefronts}")
  minTimeForAWavefrontToReachNextNodeOnPath = min([wv['timetilldestinationnodeonlink'] for wv in propagatingwavefronts])  # initial candidate for the next time-step's duration, during which all wavefronts will propagate; at the end of the minimum time just calculated, at least one of the wavefronts will intersect with the next node on its shortest path to an exit; TODO: finding the minimum time till destination-node takes O(W) = O(N^2) time, where W is number of wavefronts, if the wavefronts aren't maintained in order of 'timetilldestinationnodeonlink', so consider keeping them in that order so that the minimum time till destination-node can be found in O(1) time.
  print(f"minTimeForAWavefrontToReachNextNodeOnPath is {minTimeForAWavefrontToReachNextNodeOnPath} hours.")
  wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes = [wv for wv in propagatingwavefronts if wv['timetilldestinationnodeonlink'] == minTimeForAWavefrontToReachNextNodeOnPath]  # every wavefront that will reach, at the end of the minimum time just calculated, the next node on its shortest path to an exit
#  assert len(wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes) == 1  # TODO: modify this assertion if it turns out to be not always true, which it probably will eventually; it's likely that len(wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes) can be any integer >= 1
##  nextwavefronttopropagate = wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes[0]  # this wavefront, in the next time-step, will either reach the next node on its shortest path to an exit or will first intersect with another wavefront propagating on the same link  #(in which case both wavefronts are probably upstream-propagating shock-fronts; TODO: is probably not generally true, but will assert it as a check)
  nextwavefronttopropagate = wavefrontsThatTakeExactlyThatTimeToReachTheirNextNodes[0]  # if there turned out to be no intersection before then between a pair of wavefronts propagating on the same link, then this wavefront, in the next time-step, would be the first to reach the next node on its shortest path to an exit; the value of 'nextwavefronttopropagate' will change below if an imminent intersection between wavefronts is found
#  print(f"nextwavefronttopropagate, which in the next time-step will either reach the next node on its shortest path to an exit or will first intersect with another wavefront propagating on the same link, is {nextwavefronttopropagate}")
  print(f"Candidate nextwavefronttopropagate, which in the absence of an intersection between same-link wavefronts would, in the next time-step, be the first to reach its next node (on its shortest path to an exit), is {nextwavefronttopropagate}")
#  edgePropagatingOn = nextwavefronttopropagate['currentlink']
#  otherWavefrontsOnSameLink = [wv for wv in propagatingwavefronts if wv['currentlink'] == edgePropagatingOn and wv != nextwavefronttopropagate]  # TODO: finding the other wavefronts on the same link takes O(W) = O(N^2) time, where W is number of wavefronts, if the wavefronts aren't maintained in order of 'currentlink', so consider keeping them in that order so that the other wavefronts can be found in O(1) time.
#  print(f"otherWavefrontsOnSameLink is {otherWavefrontsOnSameLink}")
#  speed_nwftp = nextwavefronttopropagate['propagationspeed']  # propagation-speed of nextwavefronttopropagate
#  distreacheddownlink_nwftp = nextwavefronttopropagate['distancereacheddownlink']
#  print(f"speed_nwftp is {speed_nwftp} km/hr, distreacheddownlink_nwftp is {distreacheddownlink_nwftp} km/hr")

  # The previous method of choosing nextwavefronttopropagate was to choose a wavefront that would reach its next node in the shortest time; this was incorrect, because there might be another wavefront that's far from its next node but is just about to intersect a wavefront on its link.
  # Find next wavefront to propagate, i.e. the one that will be the soonest to either intersect with another on its link or reach its next node at the end of its link:
  wavefrontIntersected = None  # this and the two variables below will change if we find two wavefronts on the same link that will intersect
  pointReachedByNextPropagatingWaveFront = None
  minTimeForAWavefrontToIntersectWithAnotherOnSameLink = LARGEVAL
  for edgePropagatingOn in argsflowfunc['lK'].keys():
    upstreamPropagatingWavefrontsOnLink = [wv for wv in propagatingwavefronts if wv['currentlink'] == edgePropagatingOn and wv['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front']
    if upstreamPropagatingWavefrontsOnLink:
      print(f"upstreamPropagatingWavefrontsOnLink {edgePropagatingOn} is {upstreamPropagatingWavefrontsOnLink}")
    if not upstreamPropagatingWavefrontsOnLink: continue
    wavefrontsOnLink = [wv for wv in propagatingwavefronts if wv['currentlink'] == edgePropagatingOn]
    print(f"Wavefronts on link {edgePropagatingOn} are {wavefrontsOnLink}.")
#    if not wavefrontsOnLink: continue
    L = argsflowfunc['linklength'][ edgePropagatingOn ]
#    minTimeForAWavefrontOnThisLinkToReachNextNode = min([wv['timetilldestinationnodeonlink'] for wv in wavefrontsOnLink])  # an upper bound on the time it can take until the next intersection of wavefronts, if any, on this link; if it took longer than this then a wavefront would have already reached its next node
##    timeTillIntersected = LARGEVAL  # this will change if we find two wavefronts on the same link that shall soon intersect
#    upstreamPropagatingWavefrontsOnLink = [wv for wv in wavefrontsOnLink if wv['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front']
#    print(f"upstreamPropagatingWavefrontsOnLink is {upstreamPropagatingWavefrontsOnLink}")
    # Generate each pair of wavefronts that will potentially intersect (no pair will be generated more than once by the following procedure):
#    for thiswf in wavefrontsOnLink:
    for thiswf in upstreamPropagatingWavefrontsOnLink:  # any pair of intersecting wavefronts will contain at least one that's propagating upstream
      assert L == argsflowfunc['linklength'][ thiswf['currentlink'] ]
      speed_thiswf = thiswf['propagationspeed']  # propagation-speed of 'thiswf'
      distreacheddownlink_thiswf = thiswf['distancereacheddownlink']
#      print(f"speed_thiswf is {speed_thiswf} km/hr, distreacheddownlink_thiswf is {distreacheddownlink_thiswf} km/hr")
      timeTillThisWFwouldReachNextNodeOnPath = thiswf['timetilldestinationnodeonlink']
      assert thiswf['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front'  # TODO: as this should always be true, can delete the 'if' clause below and its code-block
      if thiswf['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':  # thiswf is downstream-propagating; can potentially intersect only with upstream-propagating wavefronts on same link
        candidateWavefrontsToIntersect = [wv for wv in wavefrontsOnLink if wv['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front']
        print(f"candidateWavefrontsToIntersect is {candidateWavefrontsToIntersect}")
        (minTimeForAWavefrontToIntersectWithAnotherOnSameLink, nextwavefronttopropagate, wavefrontIntersected, pointReachedByNextPropagatingWaveFront) = examineCandidateWavefrontsForIntersectionWithThisWavefront(thiswf, candidateWavefrontsToIntersect, nextwavefronttopropagate, wavefrontIntersected, L, speed_thiswf, distreacheddownlink_thiswf, minTimeForAWavefrontToIntersectWithAnotherOnSameLink, minTimeForAWavefrontToReachNextNodeOnPath, pointReachedByNextPropagatingWaveFront)
      else:  # thiswf is upstream-propagating; can potentially intersect only with slower upstream-propagating wavefronts on same link that are a greater distance up the link (can intersect with downstream-propagating wavefronts, but such an intersection will be found by the case above)
#        candidateWavefrontsToIntersect = [wv for wv in wavefrontsOnLink if wv['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front' or wv['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front' and wv['propagationspeed'] < thiswf['propagationspeed']]  # 'wv != thiswf' is implied by the condition on wavefronts' speeds
        candidateWavefrontsToIntersect = [wv for wv in wavefrontsOnLink if wv['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front' or wv['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front' and wv['propagationspeed'] < thiswf['propagationspeed'] and wv['distancereacheddownlink'] > thiswf['distancereacheddownlink']]  # 'wv != thiswf' is implied by the condition on wavefronts' speeds
        print(f"candidateWavefrontsToIntersect is {candidateWavefrontsToIntersect}")
        (minTimeForAWavefrontToIntersectWithAnotherOnSameLink, nextwavefronttopropagate, wavefrontIntersected, pointReachedByNextPropagatingWaveFront) = examineCandidateWavefrontsForIntersectionWithThisWavefront(thiswf, candidateWavefrontsToIntersect, nextwavefronttopropagate, wavefrontIntersected, L, speed_thiswf, distreacheddownlink_thiswf, minTimeForAWavefrontToIntersectWithAnotherOnSameLink, minTimeForAWavefrontToReachNextNodeOnPath, pointReachedByNextPropagatingWaveFront)
#        for cwf in candidateWavefrontsToIntersect:
#          # Calculate remaining time till candidate wavefront would be caught up by 'thiswf':
#          speed_cwf = cwf['propagationspeed']  # propagation-speed of candidate wavefront
#          distreacheddownlink_cwf = cwf['distancereacheddownlink']
#          print(f"speed_cwf is {speed_cwf} km/hr, distreacheddownlink_cwf is {distreacheddownlink_cwf} km")
#          timeTillIntersected = (distreacheddownlink_cwf - distreacheddownlink_thiswf) / (speed_thiswf - speed_cwf)
#          assert timeTillIntersected > 0.0
#          print(f"timeTillIntersected is {timeTillIntersected} hours")
#          if timeTillIntersected < minTimeForAWavefrontToIntersectWithAnotherOnSameLink and timeTillIntersected <= minTimeForAWavefrontToReachNextNodeOnPath:  # for the intersection to be such, it must occur by the time (before??) nextwavefronttopropagate reaches its destination-node
#            assert timeTillIntersected < minTimeForAWavefrontToReachNextNodeOnPath  # TODO: it would be convenient if this were always true, for otherwise nextwavefronttopropagate might reach its destination-node simultaneously with intersecting the other wavefront on the link, although that case might not be too hard to sort out
#            minTimeForAWavefrontToIntersectWithAnotherOnSameLink = timeTillIntersected
#            nextwavefronttopropagate = thiswf
#            wavefrontIntersected = cwf
#            pointReachedByNextPropagatingWaveFront = 'catch up to slower-propagating shockwave on same link'

      assert minTimeForAWavefrontToIntersectWithAnotherOnSameLink > 0.0

#      print(f"Wavefront {thiswf} has candidateWavefrontsToIntersect {candidateWavefrontsToIntersect}")

#  if minTimeForThisWavefrontToIntersectWithAnotherOnSameLink < LARGEVAL:
  if minTimeForAWavefrontToIntersectWithAnotherOnSameLink < LARGEVAL:
#    nexttimestepduration = minTimeForThisWavefrontToIntersectWithAnotherOnSameLink
    nexttimestepduration = minTimeForAWavefrontToIntersectWithAnotherOnSameLink
  else:  # no pair of wavefronts that shall soon intersect was found
    print(f"No pair of soon-to-intersect same-link wavefronts was found.")
#    nexttimestepduration = minTimeForAWavefrontOnThisLinkToReachNextNode
    nexttimestepduration = minTimeForAWavefrontToReachNextNodeOnPath
    pointReachedByNextPropagatingWaveFront = f'reach end-node {nextwavefronttopropagate["destinationnodeonlink"]} at end of link'

  print(f"nexttimestepduration is {nexttimestepduration} hours, pointReachedByNextPropagatingWaveFront is '{pointReachedByNextPropagatingWaveFront}'")
##  assert nexttimestepduration > 0.0
#  assert nexttimestepduration > 1e-5
  assert nexttimestepduration > 1e-6

#  print(f"nextwavefronttopropagate, which in the next time-step will either reach the next node on its shortest path to an exit or will first intersect with another wavefront propagating on the same link, is {nextwavefronttopropagate}")
  print(f"nextwavefronttopropagate, which in the next time-step will {pointReachedByNextPropagatingWaveFront}, is {nextwavefronttopropagate}")
#  sys.exit()  # for debugging

##  catchuptime = {}  # for each other wavefront on the same link, the time it would be caught up to by 'nextwavefronttopropagate'
#  speedsEqual = False  # TODO: debugging-variable: delete once working
#  wavefrontIntersected = None  # this and the variable below will change if we find a wavefront on the same link that will be intersected by nextwavefronttopropagate
#  minTimeForNextWavefrontToIntersectWithAnotherOnSameLink = LARGEVAL
#  # Out of the other wavefronts on the same link as nextwavefronttopropagate, find the one, if any, that nextwavefronttopropagate will catch up soonest:
#  for owf in otherWavefrontsOnSameLink:
#    # Calculate remaining time (might be never) till this wavefront would be intersected by 'nextwavefronttopropagate':
#    speed_owf = owf['propagationspeed']  # propagation-speed of owf (other wavefront on same link)
#    distreacheddownlink_owf = owf['distancereacheddownlink']
#    print(f"speed_owf is {speed_owf} km/hr, distreacheddownlink_owf is {distreacheddownlink_owf} km")
#    if nextwavefronttopropagate['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':  # is downstream-propagating
##      assert owf['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front'  # TODO: not generally true: a downstream-propagating wavefront could meet an upstream-propagating one on the same link
#      if owf['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':  # if nextwavefronttopropagate is downstream-propagating and so is the other wavefront, then both wavefronts must be propagating at the same speed
#        assert speed_owf == speed_nwftp
##    assert speed_owf < speed_nwftp  # the two wavefronts' speeds should never be equal, unless they are both already at the same distance down the link, which shouldn't occur  # TODO: isn't always true if both wavefronts are downstream-propagating (and perhaps isn't always true even if they're both upstream-propagating), but that should be tested for above
#    assert speed_owf <= speed_nwftp  # the two wavefronts' speeds might be equal, for example if they're both downstream-propagating  # TODO: is probably not always true: for example, the other wavefront might be propagating downstream at a faster rate than nextwavefronttopropagate is propagating upstream
##    assert distreacheddownlink_owf > distreacheddownlink_nwftp  # not always true if if both wavefronts are downstream-propagating, and in fact isn't always true even if they're both upstream-propagating
###    assert speed_owf != speed_nwftp  # the two wavefronts' speeds should never be equal, unless they are both already at the same distance down the link, which shouldn't occur
#    if nextwavefronttopropagate['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front' and owf['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front' or nextwavefronttopropagate['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front' and owf['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':  # if nextwavefronttopropagate and the other wavefront are both propagating in the same direction, either downstream or upstream
#      if speed_nwftp == speed_owf:
#        timeTillIntersected = LARGEVAL  # will never catch up
#      else:
#        timeTillIntersected = (distreacheddownlink_owf - distreacheddownlink_nwftp) / (speed_nwftp - speed_owf)
#    else:  # if nextwavefronttopropagate and the other wavefront are propagating in opposite directions on the link
#      timeTillIntersected = (L - distreacheddownlink_nwftp - distreacheddownlink_owf) / (speed_nwftp + speed_owf)
#    print(f"timeTillIntersected is {timeTillIntersected} hours")
#    assert timeTillIntersected > 0.0
###    catchuptime[owf] = timeTillIntersected
##    if timeTillIntersected < minTimeForNextWavefrontToIntersectWithAnotherOnSameLink and timeTillIntersected <= minTimeForAWavefrontToReachNextNodeOnPath:  # for the intersection to be such, it must occur by the time (before??) nextwavefronttopropagate reaches its destination-node
##      assert timeTillIntersected < minTimeForAWavefrontToReachNextNodeOnPath  # TODO: it would be convenient if this were always true, for otherwise nextwavefronttopropagate might reach its destination-node simultaneously with intersecting the other wavefront on the link, although that case might not be too hard to sort out
#    if timeTillIntersected < minTimeForNextWavefrontToIntersectWithAnotherOnSameLink and timeTillIntersected <= timeTillThisWFwouldReachNextNodeOnPath:  # for the intersection to be such, it must occur by the time (before??) nextwavefronttopropagate reaches its destination-node
#      assert timeTillIntersected < timeTillThisWFwouldReachNextNodeOnPath  # TODO: it would be convenient if this were always true, for otherwise nextwavefronttopropagate might reach its destination-node simultaneously with intersecting the other wavefront on the link, although that case might not be too hard to sort out
#      minTimeForNextWavefrontToIntersectWithAnotherOnSameLink = timeTillIntersected
#      wavefrontIntersected = owf
###    if catchuptime:
##    print("catchuptime is", catchuptime)
##    wavefrontIntersected = min(catchuptime, catchuptime.get)  # obtain key in dictionary 'catchuptime' having the minimum value; the default for 'min' is to give the minimum key, but here the keys are ordered by the function catchuptime.get, i.e. ordered by the keys' dictionary-values (see https://docs.python.org/3/library/functions.html#min)
##    minTimeForNextWavefrontToIntersectWithAnotherOnSameLink = catchuptime[wavefrontIntersected]
#    print(f"minTimeForNextWavefrontToIntersectWithAnotherOnSameLink is {minTimeForNextWavefrontToIntersectWithAnotherOnSameLink} hours")
#    print(f"wavefrontIntersected is {wavefrontIntersected}")
#    assert minTimeForNextWavefrontToIntersectWithAnotherOnSameLink > 0.0
##    assert nextwavefronttopropagate['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable'  # TODO: this assertion might prove to be not generally correct
##    assert wavefrontIntersected['nodesonpathtoanexit'] == 'upstream-propagating shock-front, so not applicable'  # TODO: this assertion might prove to be not generally correct
#    if minTimeForNextWavefrontToIntersectWithAnotherOnSameLink < LARGEVAL:
#      nexttimestepduration = minTimeForNextWavefrontToIntersectWithAnotherOnSameLink
#      pointReachedByNextPropagatingWaveFront = 'catch up to slower-propagating shockwave on same link'
#      print(f"nexttimestepduration is {nexttimestepduration} hours, pointReachedByNextPropagatingWaveFront is '{pointReachedByNextPropagatingWaveFront}'")
##    else:
##      print(f"nexttimestepduration is {nexttimestepduration} hours, pointReachedByNextPropagatingWaveFront is '{pointReachedByNextPropagatingWaveFront}'")
#    elif speed_nwftp == speed_owf:  # these two or three lines are debugging-code: delete eventually
#      speedsEqual = True

##  else:  # there are no other wavefronts propagating on the same link as 'nextwavefronttopropagate'
#  if wavefrontIntersected == None:  # if no other wavefront propagating on the same link as 'nextwavefronttopropagate' will be intersected by it
##    nexttimestepduration = min([wv['timetilldestinationnodeonlink'] for wv in propagatingwavefronts])  # the next time-step's duration, during which all wavefronts will propagate; at the end of the next time-step one (or more) of the wavefronts will intersect with the next node on its shortest path to an exit; TODO: finding the next time-step's duration takes O(W) = O(N^2) time, where W is number of wavefronts, if the wavefronts aren't maintained in order of 'timetilldestinationnodeonlink', so consider keeping them in that order so that the next time-step's duration can be found in O(1) time.
##    nexttimestepduration = timeTillChosenWavefrontWouldReachNextNodeOnPath
#    nexttimestepduration = minTimeForAWavefrontToReachNextNodeOnPath
#    pointReachedByNextPropagatingWaveFront = f'reach end-node {nextwavefronttopropagate["destinationnodeonlink"]} at end of link'
#    if speedsEqual and speed_nwftp == speed_owf:  # these four lines are debugging-code: delete eventually
#      print(f"speed_nwftp {speed_nwftp} and speed_owf {speed_owf} are equal.")
#      print(f"nexttimestepduration is {nexttimestepduration} hours, pointReachedByNextPropagatingWaveFront is '{pointReachedByNextPropagatingWaveFront}'")
  return (nexttimestepduration, pointReachedByNextPropagatingWaveFront, nextwavefronttopropagate, wavefrontIntersected)


def printSimulationClockTime(timestepcount, simulationclocktime):
  print(f"==============Time-step {timestepcount}: simulationclocktime is {simulationclocktime} hours ({simulationclocktime*3600} seconds).")


def verifyOnEveryLinkThatPropagatingFlowMatchesqijValues(propagatingwavefronts, argsflowfunc, qijatnodei, qijatnodej ):
  for edgePropagatingOn in argsflowfunc['lK'].keys():
    downstreamPropagatingWavefrontsOnLinkDict = {wv['distancereacheddownlink']: wv for wv in propagatingwavefronts if wv['currentlink'] == edgePropagatingOn and wv['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front'}
    print(f"downstreamPropagatingWavefrontsOnLinkDict {edgePropagatingOn} is {downstreamPropagatingWavefrontsOnLinkDict}")
##    wavefrontIntersected = min(catchuptime, catchuptime.get)  # obtain key in dictionary 'catchuptime' having the minimum value; the default for 'min' is to give the minimum key, but here the keys are ordered by the function catchuptime.get, i.e. ordered by the keys' dictionary-values (see https://docs.python.org/3/library/functions.html#min)
    if downstreamPropagatingWavefrontsOnLinkDict:
      leastDownstreamWavefrontDistance = min(downstreamPropagatingWavefrontsOnLinkDict)  # 'min' obtains key in dictionary 'downstreamPropagatingWavefrontsOnLinkDict' having the minimum value, i.e. distancereacheddownlink
      farthestDownstreamWavefrontDistance = max(downstreamPropagatingWavefrontsOnLinkDict)
      leastDownstreamWavefront = downstreamPropagatingWavefrontsOnLinkDict[ leastDownstreamWavefrontDistance ]
      farthestDownstreamWavefront = downstreamPropagatingWavefrontsOnLinkDict[ farthestDownstreamWavefrontDistance ]
#      print(f"farthestDownstreamWavefront is {farthestDownstreamWavefront}")
      print(f"leastDownstreamWavefrontDistance is {leastDownstreamWavefrontDistance}")
      print(f"farthestDownstreamWavefrontDistance is {farthestDownstreamWavefrontDistance}")
      a = leastDownstreamWavefront['sourcenodeonlink']
      b = leastDownstreamWavefront['destinationnodeonlink']
      print(f"leastDownstreamWavefront['flow'] is {leastDownstreamWavefront['flow']}, qijatnodei[{a},{b}] is {qijatnodei[a,b]}")
      assert leastDownstreamWavefront['flow'] == qijatnodei[a,b]
      flowAheadOfFarthestDownstreamWavefront = farthestDownstreamWavefront['flowaheadofwavefront']
    else:
      leastDownstreamWavefront = None
      farthestDownstreamWavefront = None
      flowAheadOfFarthestDownstreamWavefront = None
    upstreamPropagatingWavefrontsOnLinkDict = {wv['distancereacheddownlink']: wv for wv in propagatingwavefronts if wv['currentlink'] == edgePropagatingOn and wv['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front'}
    if upstreamPropagatingWavefrontsOnLinkDict:
      print(f"upstreamPropagatingWavefrontsOnLinkDict {edgePropagatingOn} is {upstreamPropagatingWavefrontsOnLinkDict}")
    if upstreamPropagatingWavefrontsOnLinkDict:
      leastUpstreamShockfrontDistance = min(upstreamPropagatingWavefrontsOnLinkDict)  # 'min' obtains key in dictionary 'upstreamPropagatingWavefrontsOnLinkDict' having the minimum value, i.e. distancereacheddownlink
      farthestUpstreamShockfrontDistance = max(upstreamPropagatingWavefrontsOnLinkDict)
      leastUpstreamShockfront = upstreamPropagatingWavefrontsOnLinkDict[ leastUpstreamShockfrontDistance ]
      farthestUpstreamShockfront = upstreamPropagatingWavefrontsOnLinkDict[ farthestUpstreamShockfrontDistance ]
#      print(f"farthestUpstreamShockfront is {farthestUpstreamShockfront}")
      print(f"leastUpstreamShockfrontDistance is {leastUpstreamShockfrontDistance}")
      print(f"farthestUpstreamShockfrontDistance is {farthestUpstreamShockfrontDistance}")
      a = leastUpstreamShockfront['destinationnodeonlink']
      b = leastUpstreamShockfront['sourcenodeonlink']
      print(f"leastUpstreamShockfront['flow'] is {leastUpstreamShockfront['flow']}, qijatnodej[{a},{b}] is {qijatnodej[a,b]}")
      assert leastUpstreamShockfront['flow'] == qijatnodej[a,b]
      flowAheadOfFarthestUpstreamShockfront = farthestUpstreamShockfront['flowaheadofwavefront']
    else:
      leastUpstreamShockfront = None
      farthestUpstreamShockfront = None
      flowAheadOfFarthestUpstreamShockfront = None
#    if not upstreamPropagatingWavefrontsOnLinkDict: continue
#    L = argsflowfunc['linklength'][ edgePropagatingOn ]
    if flowAheadOfFarthestDownstreamWavefront and flowAheadOfFarthestUpstreamShockfront:
      assert flowAheadOfFarthestDownstreamWavefront == flowAheadOfFarthestUpstreamShockfront


#def evolutionofwaveandshockfrontsSEM4( inflowsandflowperiodsByNodeID, subflowsToAssignedExitNodesByInjectionNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, densityijatnodej, speedijatnodej, flowstateijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, shortestPathToExit, flowijfunc, flowfunc, injectionnodes, outflowsByNodeID, pointcoordsandindexbyID, pointNodeIDbyIndex, nextNodeInShortestPath, directedNeighboursInShortestPathTreeOfNode, assignedExitNodesReachedDownstreamOfArc ):
def evolutionofwaveandshockfrontsSEM4( inflowsandflowperiodsByNodeID, subflowsToAssignedExitNodesByInjectionNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, densityijatnodej, speedijatnodej, flowstateijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, shortestPathToExit, flowijfunc, flowfunc, injectionnodes, outflowsByNodeID, pointcoordsandindexbyID, pointNodeIDbyIndex, nextNodeInShortestPath, successorNodeOnShortestPathFromNodeToExitNode, assignedExitNodesReachedDownstreamOfArc ):
  EPS = 1e-6
  INITIALTIMESTAMP = 0.0
  congestionAlertMessages = []  # will change if congestion is detected at an injection-node
#  propagatingwavefronts = {}
#  for edge in argsflowfunc['lK']:
#    propagatingwavefronts[edge] = []
  propagatingwavefronts = []
#  print(f"inflowsandflowperiodsByNodeID.values() is {inflowsandflowperiodsByNodeID.values()}")
  totalinflows = sum(u['inflow'] for u in inflowsandflowperiodsByNodeID.values())
  print(f"totalinflows is {totalinflows} veh/hr.")
#  sys.exit()
  print(f"inflowsandflowperiodsByNodeID is {inflowsandflowperiodsByNodeID}")
  # Create at time t=0 the initial wavefronts at injection-nodes:
  for nodeIDa in inflowsandflowperiodsByNodeID:
    print(f"For node {nodeIDa}, inflowsandflowperiodsByNodeID['{nodeIDa}'] is {inflowsandflowperiodsByNodeID[nodeIDa]}, subflowsToAssignedExitNodesByInjectionNodeID[{nodeIDa}] is {subflowsToAssignedExitNodesByInjectionNodeID[nodeIDa]}")
#  totalinflows = sum(d['inflow'] for d in inflowsandflowperiodsByNodeID.values())
    inflowata = inflowsandflowperiodsByNodeID[nodeIDa]['inflow']
    assert sum(d for d in subflowsToAssignedExitNodesByInjectionNodeID[nodeIDa].values()) == inflowata  # the subflows assigned to exit-nodes sum to the injection-flow
    a = pointcoordsandindexbyID[nodeIDa]['index'] 
    if inflowata <= EPS:
      print(f"WARNING: inflow of {inflowata} veh/hr at injection-node {a} is too close to zero: no downstream-propagating wavefront is created at {a}.")
      continue
    if a not in shortestPathToAnyExit:
      print(f"WARNING: there is no path from injection-node {a} along directed arcs to an exit-node: no downstream-propagating wavefront is created at {a}.")  # TODO: should traffic at node a proceed as far into the network as it can go? (I think not, as traffic won't be sent into the network if it can't reach an exit-node; if it were sent into the network, it might cause other traffic to become congested.)
      congestionAlertMessage = f"ALERT: at injection-node {a}, inflow of {inflowata} veh/hr has no path to an exit-node, so cannot be assigned to one."
      print(congestionAlertMessage)
      congestionAlertMessages.append( congestionAlertMessage )
      assert inflowata > EPS
      continue
#    shortestpathtoanexit_froma = shortestPathToAnyExit[a]
    # For each successor-node of a on a shortest-path, construct a list of the exit-nodes that traffic at a is assigned to that are reached by going next to b:
    exitNodesReachedViaSuccessorNode = {}
    for exitnode in subflowsToAssignedExitNodesByInjectionNodeID[nodeIDa]:
      e = pointcoordsandindexbyID[exitnode]['index']
      successorNodeToReachThisExitNode = nextNodeInShortestPath[a,e]
      if successorNodeToReachThisExitNode in exitNodesReachedViaSuccessorNode:
        exitNodesReachedViaSuccessorNode[successorNodeToReachThisExitNode].append( e )
      else:
        exitNodesReachedViaSuccessorNode[successorNodeToReachThisExitNode] = [e]
    print(f"exitNodesReachedViaSuccessorNode for injection-node {nodeIDa} is {exitNodesReachedViaSuccessorNode}")
#    b = shortestpathtoanexit_froma[1]  # b is the next node in the shortest path from a to an exit-node
    for b in exitNodesReachedViaSuccessorNode.keys():  # b is the next node in the shortest path from a to an exit-node
      initialflowab = 0.0
      print(f"initial flow[{a}, {b}] is {initialflowab} veh/hr")
      assert initialflowab == 0.0  # FIXME: should check for closeness to zero of flow's absolute value
      subflowsOnPathsToExits = {e: {'subflow': subflowsToAssignedExitNodesByInjectionNodeID[nodeIDa][pointNodeIDbyIndex[e]], 'shortestpath': shortestPathToExit[a,e]} for e in exitNodesReachedViaSuccessorNode[b]}
      print(f"subflowsOnPathsToExits is {subflowsOnPathsToExits}")
      linkcapacity = argsflowfunc['linkcapacity'][frozenset({a,b})]
      print(f"linkcapacity is {linkcapacity} veh/hr")
      sumOfSubflowsOnCurrentLink = sum(subflowsOnPathsToExits[e]['subflow'] for e in subflowsOnPathsToExits)
      print(f"sumOfSubflowsOnCurrentLink is {sumOfSubflowsOnCurrentLink} veh/hr.")
      if sumOfSubflowsOnCurrentLink <= linkcapacity:
        upstreamfreeflowab = sumOfSubflowsOnCurrentLink
      else:
        upstreamfreeflowab = linkcapacity
        # If inflow exceeds capacity, report an upstream-propagating shockwave created at the injection-node that affects only the hypothetical injection-"link" by which injection-traffic might be supposed to enter the network:
#        congestionAlertMessage = f"ALERT: at injection-node node {a} at time {INITIALTIMESTAMP} hours, linkcapacity {linkcapacity} < injection-flow {inflowsandflowperiodsByNodeID[nodeIDa]['inflow']}: there is a shock where traffic-flow enters the network there."
        flowReductionFactor = linkcapacity / sumOfSubflowsOnCurrentLink  # subflows on link are reduced by the same factor
        print(f"flowReductionFactor is linkcapacity / sumOfSubflowsOnCurrentLink = {linkcapacity} / {sumOfSubflowsOnCurrentLink} = {flowReductionFactor}")
        for e in exitNodesReachedViaSuccessorNode[b]:
          subflowsOnPathsToExits[e]['subflow'] *= flowReductionFactor
        print(f"subflowsOnPathsToExits is now {subflowsOnPathsToExits}")
#          congestedinjectionflowata = flowReductionFactor * injectionflowata
        congestionAlertMessage = f"ALERT: at injection-node {a} ('{pointNodeIDbyIndex[a]}') into arc ({a},{b}) at time {INITIALTIMESTAMP} hours, link-capacity {linkcapacity} veh/hr < arc's assigned inflow {sumOfSubflowsOnCurrentLink}: there is a shock where traffic-flow enters the network at node {a} into arc ({a},{b}), and the assigned inflow in that arc is reduced by a factor of {linkcapacity}/{sumOfSubflowsOnCurrentLink} = {flowReductionFactor}."
        print(congestionAlertMessage)
        congestionAlertMessages.append( congestionAlertMessage )
      print(f"Upstream free flow[{a}, {b}] is {upstreamfreeflowab} veh/hr")
#      print(f"Inflow of {inflowsandflowperiodsByNodeID[nodeIDa]['inflow']} veh/hr at node {a} causes a downstream-propagating wavefront to form, with the inflow - perhaps restricted by capacity {linkcapacity} veh/hr of link ({a},{b}) - freely flowing upstream of the wavefront, and zero flow downstream of it.")
      print(f"Inflow of {upstreamfreeflowab} veh/hr into arc ({a},{b}) at node {a} causes a downstream-propagating wavefront to form, with the inflow - perhaps restricted by capacity {linkcapacity} veh/hr of link ({a},{b}) - freely flowing upstream of the wavefront, and zero flow downstream of it.")
      downstreamflowab = initialflowab
      print(f"Downstream flow[{a}, {b}] is {downstreamflowab} veh/hr")
      assert sum(subflowsOnPathsToExits[e]['subflow'] for e in subflowsOnPathsToExits) + downstreamflowab == upstreamfreeflowab  # the subflows assigned to exit-nodes reached via the current link sum to the upstream free flow on the link
#      if inflowsandflowperiodsByNodeID[nodeIDa]['inflow'] <= linkcapacity:
#        upstreamfreeflowab = inflowsandflowperiodsByNodeID[nodeIDa]['inflow']
#      else:
#        upstreamfreeflowab = linkcapacity
      edge = frozenset({a, b})
      upstreamfreedensityab = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamfreeflowab, edge, argsflowfunc )
      downstreamdensityab = densityonfreeflowbranchoftriangularfundamentaldiagram( downstreamflowab, edge, argsflowfunc )
      print(f"downstreamdensityab is {downstreamdensityab} veh/km")
      assert downstreamdensityab == 0.0
#      if abs(downstreamflowab) <= EPS:
##        assert (a, b) not in qijatnodej
##        assert (b, a) not in qijatnodej
#        qijatnodej[a, b] = downstreamflowab
#        print(f"qijatnodej[{a}, {b}] is now {qijatnodej[a, b]}")
##      densityijatnodej[edge] = downstreamdensityab
##      print(f"densityijatnodej is now {densityijatnodej} veh/km")
      qijatnodei[a,b] = upstreamfreeflowab
      densityijatnodei[a,b] = upstreamfreedensityab
      print(f"qijatnodei[{a},{b}] is {qijatnodei[a,b]} veh/hr.")
      print(f"As wavefront propagates downstream, flow in link [{a}, {b}] will increase from {initialflowab} to {upstreamfreeflowab} veh/hr and density increase from {downstreamdensityab} to {upstreamfreedensityab} veh/km.")
      print(f"Calculating propagation-speed: upstreamfreeflowab is {upstreamfreeflowab}, upstreamfreedensityab is {upstreamfreedensityab}, downstreamflowab is {downstreamflowab}, downstreamdensityab is {downstreamdensityab}")
      wavefrontPropagationSpeed = propagationSpeedOfShockFront( upstreamfreeflowab, upstreamfreedensityab, downstreamflowab, downstreamdensityab )
      if wavefrontPropagationSpeed < 0.0:
        upordownstreamText = 'will travel upstream'
      elif wavefrontPropagationSpeed > 0.0:
        upordownstreamText = 'will travel downstream'
      else:
        upordownstreamText = 'is stationary'
      print(f"wavefrontPropagationSpeed is {wavefrontPropagationSpeed} km/hr: the front {upordownstreamText}.")
      assert wavefrontPropagationSpeed >= 0.0
      if wavefrontPropagationSpeed == 0.0:  # can be zero if, for example, the inflow from node a into link (a,b) (upstream of the "wavefront") is zero, because the downstream flow is of course also zero
        print(f"upstreamfreeflowab {upstreamfreeflowab} = downstreamflowab {downstreamflowab} (veh/hr): there is no real propagating wavefront.")
      assert wavefrontPropagationSpeed >= EPS
      assert wavefrontPropagationSpeed == argsflowfunc['lK'][ frozenset({a,b}) ]
#      speedijatnodej[edge] = upstreamfreeflowab / upstreamfreedensityab
      upstreamfreespeedinedgeab = upstreamfreeflowab / upstreamfreedensityab
      assert upstreamfreespeedinedgeab == argsflowfunc['lK'][frozenset({a,b})]
#      print(f"Density in link {edge} is now {upstreamfreedensityab} veh/km, speed is now {upstreamfreespeedinedgeab} km/hr.")
      print(f"Traffic-speed in link {edge} will increase from zero to {upstreamfreespeedinedgeab} km/hr.")
      lengthofcurrentlink = argsflowfunc['linklength'][ frozenset({a,b}) ]
      timetilldestinationnodeonlink = lengthofcurrentlink / wavefrontPropagationSpeed
####      if abs(wavefrontPropagationSpeed) > 0.0:  # record only the wavefronts with non-zero propagation-speed
###      propagatingwavefronts.append( {'time': INITIALTIMESTAMP, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': a, 'destinationnodeonlink': b, 'distancereacheddownlink': 0.0, 'propagationspeed': wavefrontPropagationSpeed, 'flow': upstreamfreeflowab, 'density': upstreamfreedensityab, 'timetilldestinationnodeonlink': timetilldestinationnodeonlink, 'nodesonpathtoanexit': shortestpathtoanexit_froma, 'flowaheadofwavefront': 0.0, 'densityaheadofwavefront': 0.0} )
##      for e in exitNodesReachedViaSuccessorNode[b]:  # b is the next node in the shortest path from a to the currently-considered exit-node e
##        shortestpathtocurrentexit_froma = shortestPathToExit[a,e]
#      propagatingwavefronts.append( {'time': INITIALTIMESTAMP, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': a, 'destinationnodeonlink': b, 'distancereacheddownlink': 0.0, 'propagationspeed': wavefrontPropagationSpeed, 'flow': upstreamfreeflowab, 'density': upstreamfreedensityab, 'timetilldestinationnodeonlink': timetilldestinationnodeonlink, 'subflowsonpathstoexits': subflowsOnPathsToExits, 'flowaheadofwavefront': 0.0, 'densityaheadofwavefront': 0.0} )
      propagatingwavefronts.append( {'time': INITIALTIMESTAMP, 'injectionNodeOfOrigin': a, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': a, 'destinationnodeonlink': b, 'distancereacheddownlink': 0.0, 'propagationspeed': wavefrontPropagationSpeed, 'flow': upstreamfreeflowab, 'density': upstreamfreedensityab, 'timetilldestinationnodeonlink': timetilldestinationnodeonlink, 'nodesonpathtoanexit': 'downstream-propagating wavefront', 'subflowsonpathstoexits': subflowsOnPathsToExits, 'flowaheadofwavefront': 0.0, 'densityaheadofwavefront': 0.0} )  # N.B. Each wavefront trails a flow behind it equal to the flow ahead plus the wavefront's (additive) subflows.

  print(f"propagatingwavefronts is {propagatingwavefronts}")
#  sys.exit()

  simulationclocktime = INITIALTIMESTAMP
  timestepcount = 0
  printSimulationClockTime(timestepcount, simulationclocktime)
  cumulNumVehiclesToHaveDepartedLinkDuringSimulation = {}
  for (i, j) in qijatnodej:
    cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j] = 0

  while simulationclocktime < simulDuration and propagatingwavefronts:
    for pwf in propagatingwavefronts:  # confirm that no wavefront, whether upstream- or downstream-propagating, has the same density behind and ahead of it
      assert abs(pwf['densityaheadofwavefront'] - pwf['density']) > EPS
#    downstreamPropagatingWavefronts = [wv for wv in propagatingwavefronts if wv['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front']
#    for dpwf in downstreamPropagatingWavefronts:  # confirm that every downstream-propagating wavefront has zero flow and density ahead of it - but this is NOT generally true
#      assert abs(dpwf['flowaheadofwavefront']) <= EPS  # not generally true
#      assert abs(dpwf['densityaheadofwavefront']) <= EPS  # not generally true
    (nexttimestepduration, pointReachedByNextPropagatingWaveFront, nextwavefronttopropagate, wavefrontIntersected) = calculateNextTimestepDuration( propagatingwavefronts, argsflowfunc )  # see whether, in the next time-step, nextwavefronttopropagate will reach the next node on its shortest path to an exit or will first intersect with another wavefront propagating on the same link
    if nextwavefronttopropagate['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':
      injectionNodeWhereFlowOriginated = nextwavefronttopropagate['injectionNodeOfOrigin']
#    print(f"nexttimestepduration is {nexttimestepduration} hours, pointReachedByNextPropagatingWaveFront is '{pointReachedByNextPropagatingWaveFront}'")
#    print(f"nextwavefronttopropagate is {nextwavefronttopropagate}")
    if simulationclocktime + nexttimestepduration > simulDuration:
      nexttimestepduration = simulDuration - simulationclocktime
      pointReachedByNextPropagatingWaveFront = f'travel until the simulation clock-time reaches simulDuration of {simulDuration} hours'
      print(f"nexttimestepduration, adjusted to not exceed simulDuration of {simulDuration} hours, is {nexttimestepduration} hours")
    # Update the cumulative number of vehicles to have departed each link so far during the simulation:
    for (i, j) in qijatnodej:
      print(f"Current flow on link ({i}, {j}) between most downstream wavefront in that link and end of link at node {j} is qijatnodej[{i}, {j}]={qijatnodej[i, j]}.")
      numVehiclesDepartingLinkDuringThisTimestep = qijatnodej[i,j] * nexttimestepduration
      print(f"numVehiclesDepartingLinkDuringThisTimestep is {numVehiclesDepartingLinkDuringThisTimestep}")
      cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j] += numVehiclesDepartingLinkDuringThisTimestep
      print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation[{i}, {j}] is now {cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j]}")
#      print(f"qijatnodei[{i},{j}] is {qijatnodei[i,j]}")
#      assert qijatnodei[i,j] <= totalinflows
#      print(f"qijatnodej[{i},{j}] is {qijatnodej[i,j]}")
#      assert qijatnodej[i,j] <= totalinflows

    print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation is {cumulNumVehiclesToHaveDepartedLinkDuringSimulation}")
    print(f"qijatnodei is {qijatnodei}")
    print(f"qijatnodej is {qijatnodej}")

    assert nextwavefronttopropagate['propagationspeed'] > 0.0  # wavefront's propagation-speed in direction of its travel
    a = nextwavefronttopropagate['sourcenodeonlink']
    b = nextwavefronttopropagate['destinationnodeonlink']
    edge = frozenset({a, b})
    currentlink = nextwavefronttopropagate['currentlink']
#    print(f"Wavefront propagates away from node {a} on link {currentlink} at speed {nextwavefronttopropagate['propagationspeed']} km/hr, taking {nexttimestepduration} hours to reach node {b}.")
    print(f"Wavefront propagates away from node {a} on link {currentlink} at speed {nextwavefronttopropagate['propagationspeed']} km/hr, taking {nexttimestepduration} hours to {pointReachedByNextPropagatingWaveFront}.")
    wavefrontflowab = nextwavefronttopropagate['flow']
    nextwavefronttopropagate['time'] = simulationclocktime + nexttimestepduration
    lengthofcurrentlink = argsflowfunc['linklength'][ currentlink ]
    print(f"lengthofcurrentlink {currentlink} is {lengthofcurrentlink} km.")
#    if pointReachedByNextPropagatingWaveFront[:22] == 'reach destination-node':
    if pointReachedByNextPropagatingWaveFront[:14] == 'reach end-node':
      # The chosen next wavefront to propagate reaches its current link's end-node at the end of this time-step, so propagate it further, either up- or downstream, into the shortest-path tree:
      if nextwavefronttopropagate['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front':
##        numnodesremainingtobetraversedonpath = 0
#        numnodesremainingtobetraversedonpath = 'upstream-propagating shock-front, so not applicable'
        nodesRemainToBeTraversedOnPath = False
      else:  # nextwavefronttopropagate is downstream-propagating
###        print(f"nextwavefronttopropagate['nodesonpathtoanexit'] is {nextwavefronttopropagate['nodesonpathtoanexit']}")
##        indexInPathOfCurrentNode = nextwavefronttopropagate['nodesonpathtoanexit'].index(b)
##        numnodesremainingtobetraversedonpath = len(nextwavefronttopropagate['nodesonpathtoanexit']) - 1 - indexInPathOfCurrentNode
#        print(f"nextwavefronttopropagate['subflowsonpathstoexits'] is {nextwavefronttopropagate['subflowsonpathstoexits']}")
        if len(nextwavefronttopropagate['subflowsonpathstoexits'].keys()) == 1 and b in nextwavefronttopropagate['subflowsonpathstoexits']:  # if there's only exit-node headed for by flow behind nextwavefronttopropagate and that exit-node is b, then there are no remaining nodes to traverse on the path  # TODO: we assume that no exit-node in 'subflowsonpathstoexits' occurs on the shortest path to another exit-node, as otherwise the traffic carrying on to the second exit-node might as well leave the network at the first exit-node; this assumption perhaps isn't always true.
          nodesRemainToBeTraversedOnPath = False
        else:
          nodesRemainToBeTraversedOnPath = True
#      if isinstance(numnodesremainingtobetraversedonpath, int) and numnodesremainingtobetraversedonpath > 0:  # update downstream-propagating wavefront as it passes through node b to the next link on its path
      if nextwavefronttopropagate['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front' and nodesRemainToBeTraversedOnPath:  # update downstream-propagating wavefront as it passes through node b to the next link on its path
        assert nextwavefronttopropagate['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front'
#        wavefrontflowdifferenceab = nextwavefronttopropagate['flow'] - nextwavefronttopropagate['flowaheadofwavefront']  # amount by which, when nextwavefronttopropagate reaches node b, flow leaving node b towards c will increase
#        c = nextwavefronttopropagate['nodesonpathtoanexit'][ indexInPathOfCurrentNode + 1]
        sumOfSubflowsPropagatingDownstreamFromNodeb = 0.0  # must calculate the total flow that will propagate downstream from node b, taking into account capacities on onwards links from b
##        for e in nextwavefronttopropagate['subflowsonpathstoexits']:
##          indexInPathOfCurrentNode = nextwavefronttopropagate['subflowsonpathstoexits'][e]['shortestpath'].index(b)
##          c = nextwavefronttopropagate['subflowsonpathstoexits'][e]['shortestpath'][ indexInPathOfCurrentNode + 1]
##          assert c == nextNodeInShortestPath[b,e]  # TODO: if this assertion always holds, then can delete this line and the two above it (might save substantial execution-time)
##          c = nextNodeInShortestPath[b,e]
##          print(f"On shortest path towards exit-node {e}, next onward node is c={c}:")
#        print(f"Directed neighbours of node {b} in shortest-path tree are {directedNeighboursInShortestPathTreeOfNode[b]}.")
        assignedExitNodesOnNextWavefront = nextwavefronttopropagate['subflowsonpathstoexits'].keys()
#        successorsOfNodebOnPathsToAssignedExitNodesOfNextWavefront = [successorNodeOnShortestPathFromNodeToExitNode[b,e] for e in assignedExitNodesOnNextWavefront]
        successorsOfNodebOnPathsToAssignedExitNodesOfNextWavefront = set()
        assignedExitNodesOnNextWavefrontReachedDownstreamOfArc = {}  # assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c] is the set of nextwavefronttopropagate's assigned exit-nodes that are downstream of arc (b,c)
        for e in assignedExitNodesOnNextWavefront:
          c = successorNodeOnShortestPathFromNodeToExitNode[b,e]
          successorsOfNodebOnPathsToAssignedExitNodesOfNextWavefront.add( c )
          if (b,c) in assignedExitNodesOnNextWavefrontReachedDownstreamOfArc:
            assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c].add( e )
          else:
            assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c] = {e}
        print(f"Successors of node {b} on shortest paths to nextwavefront's assigned exit-nodes {list(assignedExitNodesOnNextWavefront)} are {successorsOfNodebOnPathsToAssignedExitNodesOfNextWavefront}.")
#        print(f"nextwavefronttopropagate is {nextwavefronttopropagate}")
#        for c in directedNeighboursInShortestPathTreeOfNode[b]:
        for c in successorsOfNodebOnPathsToAssignedExitNodesOfNextWavefront:
#          if c == a: continue  # traffic won't reach node b then immediately turn and head back towards node a
          assert c != a  # node a cannot be a successor of node b on a shortest path to one of nextwavefront's assigned exit-nodes
#          print(f"Exit-nodes assigned to injection-node {injectionNodeWhereFlowOriginated} that are reached downstream of arc ({b},{c}) on shortest paths are {assignedExitNodesReachedDownstreamOfArc[injectionNodeWhereFlowOriginated,b,c]}.")
          print(f"Assigned exit-nodes on nextwavefronttopropagate that are reached downstream of arc ({b},{c}) on shortest paths are {assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c]}.")
          subflowabtowardsc = 0.0  # subflow which, after nextwavefronttopropagate reaches node b, is intended to flow from node b towards c bound for downstream exit-nodes, although this subflow might be reduced by insufficient capacity on arc (b,c) combined with existing flow on it
###          for e in assignedExitNodesReachedDownstreamFromNode[c]:
##          for e in assignedExitNodesReachedDownstreamOfArc[b,c]:
#          for e in assignedExitNodesReachedDownstreamOfArc[injectionNodeWhereFlowOriginated,b,c]:
          for e in assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c]:
            try:
              assert nextNodeInShortestPath[b,e] == c
            except:
              print(f"nextNodeInShortestPath[{b},{e}] {nextNodeInShortestPath[b,e]} does not equal c={c}.")
              raise
            subflowabtowardsc += nextwavefronttopropagate['subflowsonpathstoexits'][e]['subflow']  # subflow which, after nextwavefronttopropagate reaches node b, is intended to flow from node b towards c bound for exit-node e, although this subflow might be reduced by insufficient capacity on arc (b,c) combined with existing flow on it
            print(f"The subflow which, after nextwavefronttopropagate reaches node {b}, is intended to flow from node {b} towards {c} bound for assigned exit-node {e}, is {nextwavefronttopropagate['subflowsonpathstoexits'][e]['subflow']} veh/hr.")
          print(f"Total subflow which, after nextwavefronttopropagate reaches node {b}, is intended to flow from node {b} towards {c} bound for assigned exit-nodes, is {subflowabtowardsc} veh/hr.")
          nextlinkonpath = frozenset({b, c})
#          print(f"nextlinkonpath is {nextlinkonpath}")
          nextlinkcapacity = argsflowfunc['linkcapacity'][nextlinkonpath]
#          print(f"At node {b}, wavefront shall encounter capacity on next link {nextlinkonpath} of {nextlinkcapacity} veh/hr.")
          print(f"At node {b} at start of arc ({b},{c}), wavefront's subflow of {subflowabtowardsc} veh/hr shall encounter capacity on next link {nextlinkonpath} of {nextlinkcapacity} veh/hr.")
          # N.B. The wavefront trails a flow behind it equal to the flow ahead plus the wavefront's (additive) subflows.
##          if wavefrontflowab <= nextlinkcapacity - qijatnodei[b,c]:  # (incorrect) flow from a into b towards c cannot exceed edge {b,c}'s capacity less the flow already leaving b towards c
#          if wavefrontflowdifferenceab <= nextlinkcapacity - qijatnodei[b,c]:  # (incorrect) wavefront-flow minus flowahead from a into b towards c cannot exceed edge {b,c}'s capacity less the flow already leaving b towards c
          if subflowabtowardsc <= nextlinkcapacity - qijatnodei[b,c]:  # wavefront-flow from a into b towards c doesn't exceed edge {b,c}'s capacity less the flow already leaving b towards c
#            print(f"The wavefront's flow of {nextwavefronttopropagate['flow']} veh/hr less flowaheadofwavefront {nextwavefronttopropagate['flowaheadofwavefront']} is no greater than the next link's capacity of {nextlinkcapacity} less existing flow qijatnodei[{b},{c}] of {qijatnodei[b,c]} leaving {b} towards {c}, so no shock-front or bottleneck results.")
            print(f"The wavefront's subflow of {subflowabtowardsc} veh/hr on arc ({b},{c}) is no greater than the next link's capacity of {nextlinkcapacity} less existing flow qijatnodei[{b},{c}] of {qijatnodei[b,c]} leaving {b} towards {c}, so no shock-front or bottleneck results.")
            additionalFlowabc = subflowabtowardsc
            print(f"additionalFlowabc from ({a},{b}) into ({b},{c}) is {additionalFlowabc} veh/hr.")
            flowReductionFactor = 1.0  # there is no reduction to subflowabtowardsc
            sumOfSubflowsPropagatingDownstreamFromNodeb += subflowabtowardsc
            flowofwavefrontonextlink = subflowabtowardsc + qijatnodei[b,c]
            # The wavefront traverses node b and enters the next link, perhaps changing its flow, propagation-speed, and density:

##          else:  # wavefrontflowab > nextlinkcapacity - qijatnodei[b,c]
#          else:  # wavefrontflowdifferenceab > nextlinkcapacity - qijatnodei[b,c]
          else:  # subflowabtowardsc > nextlinkcapacity - qijatnodei[b,c]
##            print(f"The wavefront's flow of {nextwavefronttopropagate['flow']} veh/hr exceeds (by {nextwavefronttopropagate['flow']-nextlinkcapacity+qijatnodei[b,c]}) the next link's capacity of {nextlinkcapacity} less existing flow qijatnodei[{b},{c}] of {qijatnodei[b,c]} leaving {b} towards {c}, so a bottleneck results: an upstream-propagating shock-front is created, and the wavefront's flow will be reduced before propagating down the next link.")
#            print(f"The wavefront's flow of {nextwavefronttopropagate['flow']} veh/hr less flowaheadofwavefront {nextwavefronttopropagate['flowaheadofwavefront']} exceeds (by {nextwavefronttopropagate['flow']-nextlinkcapacity+qijatnodei[b,c]}) the next link's capacity of {nextlinkcapacity} less existing flow qijatnodei[{b},{c}] of {qijatnodei[b,c]} leaving {b} towards {c}, so a bottleneck results: an upstream-propagating shock-front is created, and the wavefront's flow will be reduced before propagating down the next link.")
            additionalFlowabc = nextlinkcapacity - qijatnodei[b,c]
            print(f"additionalFlowabc from ({a},{b}) into ({b},{c}) is {additionalFlowabc} veh/hr.")
            print(f"The wavefront's subflow of {subflowabtowardsc} veh/hr on arc ({b},{c}) exceeds (by {subflowabtowardsc-nextlinkcapacity+qijatnodei[b,c]}) the next link's capacity of {nextlinkcapacity} less existing flow qijatnodei[{b},{c}] of {qijatnodei[b,c]} leaving {b} towards {c}, so only a reduced part {additionalFlowabc} of subflow {subflowabtowardsc} veh/hr will be able to flow through from {b} into ({b},{c}); the reduced additional part {additionalFlowabc} will count towards a total congested flow that will form at the downstream end (node {b}) of ({a},{b}), where there will be a bottleneck and where an upstream-propagating shock-front shall be created.")
            flowReductionFactor = additionalFlowabc / subflowabtowardsc  # subflows on link are reduced by the same factor
            print(f"flowReductionFactor is additionalFlowabc / subflowabtowardsc = {additionalFlowabc} / {subflowabtowardsc} = {flowReductionFactor}")
            congestionAlertMessage = f"ALERT: at node {b} ('{pointNodeIDbyIndex[b]}') into arc ({b},{c}) at time {simulationclocktime + nexttimestepduration} hours, link-capacity {nextlinkcapacity} veh/hr minus qijatnodei[{b},{c}] {qijatnodei[b,c]} < arc's assigned inflow {subflowabtowardsc}: there is a shock where traffic-flow enters at node {b} into arc ({b},{c}), and the assigned inflow in that arc is reduced by a factor of ({nextlinkcapacity} - {qijatnodei[b,c]})/{subflowabtowardsc} = {flowReductionFactor}."
            print(congestionAlertMessage)
            congestionAlertMessages.append( congestionAlertMessage )
            sumOfSubflowsPropagatingDownstreamFromNodeb += additionalFlowabc
###            congestedflowab = nextlinkcapacity
##            congestedflowab = nextlinkcapacity - qijatnodei[b,c]
##            print(f"congested flow[{a}, {b}] is {congestedflowab} veh/hr")
##            print(f"As congestion-zone propagates upstream, flow in link [{a}, {b}] will decrease from {wavefrontflowab} to {congestedflowab} veh/hr")
#            nextwavefronttopropagate['flow'] = nextlinkcapacity  # flow from b towards c is now the edge {b,c}'s capacity
            flowofwavefrontonextlink = nextlinkcapacity  # flow from b towards c is now the edge {b,c}'s capacity  # TODO: must copy nextwavefronttopropagate once for each successor-node c; should probably use nextNodeInShortestPath[a,e]
##            congestedflowab = qijatnodej[a, b] * totalflowleavingnodeiDownstreamOfCongestion / (initialtotalflowEnteringNodeiFromItsParents + injectionflowati)  # TODO: divide congested flow between the links entering node b: the flow in each link, and also by implication the injection-flow, decreases proportionally
#            # Now that wavefront has reached the end of its current link, update qijatnodej, densityijatnodej, and speedijatnodej as they are at the link's departing-end:
##            assert abs(qijatnodej[a, b] - congestedflowab) > EPS  # not generally true
#            qijatnodej[a, b] = congestedflowab
##            print(f"edge is {edge}, propagatingwavefronts is {propagatingwavefronts}")
##            assert len( [wv for wv in propagatingwavefronts if wv['currentlink']==edge] ) == 1  # TODO: this assertion might prove to be not generally correct
#            assert edge == frozenset({a, b})
#            congesteddensityab = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowab, edge, argsflowfunc )
#            densityijatnodej[edge] = congesteddensityab
##            print(f"densityijatnodej is now {densityijatnodej} veh/km")
#            speedijatnodej[edge] = congestedflowab / congesteddensityab
#            flowstateijatnodej[a, b] = 'congestedflowstate'
#            print(f"qijatnodej[{a}, {b}] is now {qijatnodej[a, b]} veh/hr, densityijatnodej[{edge}] is {densityijatnodej[edge]}, speedijatnodej[{edge}] is {speedijatnodej[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
#            upstreamfreeflowab = wavefrontflowab
#            upstreamfreedensityab = densityonfreeflowbranchoftriangularfundamentaldiagram( upstreamfreeflowab, edge, argsflowfunc )
#            print(f"Calculating propagation-speed of shock-front: upstreamfreeflowab is {upstreamfreeflowab}, upstreamfreedensityab is {upstreamfreedensityab}, congestedflowab is {congestedflowab}, congesteddensityab is {congesteddensityab}")
#            shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamfreeflowab, upstreamfreedensityab, congestedflowab, congesteddensityab )
#            if shockFrontPropagationSpeed < 0.0:
#              upordownstreamText = 'travels upstream'
#            elif shockFrontPropagationSpeed > 0.0:
#              upordownstreamText = 'travels downstream'
#            else:
#              upordownstreamText = 'is stationary'
#            print(f"shockFrontPropagationSpeed is {shockFrontPropagationSpeed} km/hr: the front {upordownstreamText}")
#            assert shockFrontPropagationSpeed < 0.0
##            if shockFrontPropagationSpeed == 0.0:  # can be zero if, for example, the link's free flow upstream of the shock-front is zero, as then the "congested" flow is also zero
##              print(f"upstreamfreeflowab {upstreamfreeflowab} = congestedflowab {congestedflowab}: there is no real congestion.")
##            print(f"Density in link {edge} is now {densityijatnodej[edge]}, speed is now {speedijatnodej[edge]}.")
#            timetillupstreamnodeonlink = -lengthofcurrentlink / shockFrontPropagationSpeed
#            assert timetillupstreamnodeonlink > 0.0
#            shockfrontdict = {'time': simulationclocktime + nexttimestepduration, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': b, 'destinationnodeonlink': a, 'distancereacheddownlink': 0.0, 'propagationspeed': -shockFrontPropagationSpeed, 'flow': congestedflowab, 'density': congesteddensityab, 'timetilldestinationnodeonlink': timetillupstreamnodeonlink, 'nodesonpathtoanexit': 'upstream-propagating shock-front, so not applicable', 'flowaheadofwavefront': upstreamfreeflowab, 'densityaheadofwavefront': upstreamfreedensityab}
#            print(f"New upstream-propagating shock-front is {shockfrontdict}")
#            assert congestedflowab < upstreamfreeflowab
#            propagatingwavefronts.append( shockfrontdict )
###            print(f"After appending new upstream-propagating shock-front, propagatingwavefronts is {propagatingwavefronts}")

##            nextwavefronttopropagate['density'] = nextwavefronttopropagate['flow'] / nextwavefronttopropagate['propagationspeed']
#            if qijatnodei[b,c] < nextlinkcapacity:
#              print(f"After increasing by {nextlinkcapacity - qijatnodei[b,c]}, qijatnodei[{b},{c}] is now {nextlinkcapacity} veh/hr.")
#              qijatnodei[b,c] = nextlinkcapacity
#            assert qijatnodei[b,c] == nextlinkcapacity

          if abs(flowReductionFactor) < EPS: continue  # don't create a new downstream-propagating wavefront if it would carry no flow; continue on to b's next directed neighbour-node c
          # Update wavefront's properties on the next downstream link:
          copyOfWavefrontPropagatingOnNextArcbc = {}  # a copy of nextwavefronttopropagate is created for each successor-arc (b,c); should probably use nextNodeInShortestPath[a,e]
          copyOfWavefrontPropagatingOnNextArcbc['time'] = simulationclocktime + nexttimestepduration
          copyOfWavefrontPropagatingOnNextArcbc['injectionNodeOfOrigin'] = nextwavefronttopropagate['injectionNodeOfOrigin']
          copyOfWavefrontPropagatingOnNextArcbc['currentlink'] = nextlinkonpath
          copyOfWavefrontPropagatingOnNextArcbc['sourcenodeonlink'] = b
          copyOfWavefrontPropagatingOnNextArcbc['destinationnodeonlink'] = c
          copyOfWavefrontPropagatingOnNextArcbc['distancereacheddownlink'] = 0.0
          copyOfWavefrontPropagatingOnNextArcbc['propagationspeed'] = argsflowfunc['lK'][frozenset({b,c})]
          copyOfWavefrontPropagatingOnNextArcbc['flow'] = flowofwavefrontonextlink
          copyOfWavefrontPropagatingOnNextArcbc['density'] = copyOfWavefrontPropagatingOnNextArcbc['flow'] / copyOfWavefrontPropagatingOnNextArcbc['propagationspeed']
          lengthofnextlink = argsflowfunc['linklength'][ nextlinkonpath ]
          print(f"lengthofnextlink is {lengthofnextlink} km")
          copyOfWavefrontPropagatingOnNextArcbc['timetilldestinationnodeonlink'] = lengthofnextlink / copyOfWavefrontPropagatingOnNextArcbc['propagationspeed']
          copyOfWavefrontPropagatingOnNextArcbc['nodesonpathtoanexit'] = 'downstream-propagating wavefront'
#        for e in nextwavefronttopropagate['subflowsonpathstoexits']:
          copyOfWavefrontPropagatingOnNextArcbc['subflowsonpathstoexits'] = {}
##          for e in assignedExitNodesReachedDownstreamOfArc[injectionNodeWhereFlowOriginated,b,c]:
#          for e in assignedExitNodesOnNextWavefront:
          for e in assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c]:
#            if flowReductionFactor < 1.0:
            copyOfWavefrontPropagatingOnNextArcbc['subflowsonpathstoexits'][e] = {'subflow': flowReductionFactor*nextwavefronttopropagate['subflowsonpathstoexits'][e]['subflow'], 'shortestpath': nextwavefronttopropagate['subflowsonpathstoexits'][e]['shortestpath']}
#            if flowReductionFactor < 1.0:
#              for e in exitNodesReachedViaSuccessorNode[b]:
#                subflowsOnPathsToExits[e]['subflow'] *= flowReductionFactor
#              print(f"subflowsOnPathsToExits is now {subflowsOnPathsToExits}")
##flowofwavefrontonextlink 
##          copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] = 0.0  # the flow and density ahead of copyOfWavefrontPropagatingOnNextArcbc might be greater than zero if there was already a flow out of b towards c before this wavefront arrived at b
          copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] = qijatnodei[b,c]
#          copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] = 0.0
          copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] = densityijatnodei[b,c]
          try:
##            assert sum([copyOfWavefrontPropagatingOnNextArcbc['subflowsonpathstoexits'][e]['subflow'] for e in assignedExitNodesReachedDownstreamOfArc[injectionNodeWhereFlowOriginated,b,c]]) == copyOfWavefrontPropagatingOnNextArcbc['flow']
##            assert qijatnodei[b,c] + sum([copyOfWavefrontPropagatingOnNextArcbc['subflowsonpathstoexits'][e]['subflow'] for e in assignedExitNodesReachedDownstreamOfArc[injectionNodeWhereFlowOriginated,b,c]]) == copyOfWavefrontPropagatingOnNextArcbc['flow']
#            assert qijatnodei[b,c] + sum([copyOfWavefrontPropagatingOnNextArcbc['subflowsonpathstoexits'][e]['subflow'] for e in assignedExitNodesOnNextWavefront) == copyOfWavefrontPropagatingOnNextArcbc['flow']
            assert qijatnodei[b,c] + sum([copyOfWavefrontPropagatingOnNextArcbc['subflowsonpathstoexits'][e]['subflow'] for e in assignedExitNodesOnNextWavefrontReachedDownstreamOfArc[b,c]]) == copyOfWavefrontPropagatingOnNextArcbc['flow']
          except:
            print(f"copyOfWavefrontPropagatingOnNextArcbc is {copyOfWavefrontPropagatingOnNextArcbc}")
            raise
          if copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] != copyOfWavefrontPropagatingOnNextArcbc['density']:  # TODO: I'm not sure whether, in the case that 'densityaheadofwavefront' == 'density', copyOfWavefrontPropagatingOnNextArcbc is later removed: it ought to be
            potentialwavefrontpropagationspeed = (copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) / (copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density'])
            try:
              assert abs(potentialwavefrontpropagationspeed - copyOfWavefrontPropagatingOnNextArcbc['propagationspeed']) <= EPS  # confirm that copyOfWavefrontPropagatingOnNextArcbc['propagationspeed'], which above was defined to be the next link's free speed, is equal to the propagation-speed as calculated by formula (8.9), p. 87 of Treiber (2013)
            except:
              print(f"copyOfWavefrontPropagatingOnNextArcbc is {copyOfWavefrontPropagatingOnNextArcbc}")
              raise
#          if wavefrontflowab <= nextlinkcapacity - qijatnodei[b,c]:  # TODO: flow from a into b towards c cannot exceed edge {b,c}'s capacity less the flow already leaving b towards c
          originalqijatnodei = qijatnodei[b,c]
#          print(f"Before nextwavefronttopropagate reaches node {b}, qijatnodei[{b},{c}] is {qijatnodei[b,c]}, wavefrontflowab is {wavefrontflowab}.")
          print(f"Before nextwavefronttopropagate reaches node {b}, qijatnodei[{b},{c}] is {qijatnodei[b,c]}.")
#          qijatnodei[b,c] += wavefrontflowab
#          if qijatnodei[b,c] >= nextlinkcapacity:
#            qijatnodei[b,c] = nextlinkcapacity
#            print(f"qijatnodei[{b},{c}] is now {nextlinkcapacity} veh/hr.")
#            assert qijatnodei[b,c] == nextlinkcapacity
#          else:
#            print(f"After increasing by {wavefrontflowab}, qijatnodei[{b},{c}] is now {qijatnodei[b,c]} veh/hr.")
#            assert qijatnodei[b,c] < nextlinkcapacity
##          if qijatnodei[b,c] <= nextlinkcapacity - wavefrontflowab:
#          if qijatnodei[b,c] <= nextlinkcapacity - wavefrontflowdifferenceab:
          if qijatnodei[b,c] <= nextlinkcapacity - additionalFlowabc:
##            qijatnodei[b,c] += wavefrontflowab
#            qijatnodei[b,c] += wavefrontflowdifferenceab
            qijatnodei[b,c] += additionalFlowabc
##            print(f"After increasing by wavefrontflowab {wavefrontflowab}, qijatnodei[{b},{c}] is now {qijatnodei[b,c]} veh/hr.")
#            print(f"After increasing by wavefrontflowdifferenceab {wavefrontflowdifferenceab}, qijatnodei[{b},{c}] is now {qijatnodei[b,c]} veh/hr.")
            print(f"After increasing by additionalFlowabc {additionalFlowabc}, qijatnodei[{b},{c}] is now {qijatnodei[b,c]} veh/hr.")
            flowincrease = qijatnodei[b,c] - originalqijatnodei
            try:
#              assert flowincrease == wavefrontflowdifferenceab  # amount by which, when nextwavefronttopropagate reaches node b, flow leaving node b towards c will increase
              assert flowincrease == additionalFlowabc  # amount by which, when nextwavefronttopropagate reaches node b, flow leaving node b towards c will increase
            except:
#              print(f"flowincrease is {flowincrease}, wavefrontflowdifferenceab is {wavefrontflowdifferenceab} veh/hr.")
              print(f"flowincrease is {flowincrease}, is additionalFlowabc {additionalFlowabc} veh/hr.")
              raise
          elif qijatnodei[b,c] < nextlinkcapacity:
            qijatnodei[b,c] = nextlinkcapacity
#            print(f"qijatnodei[{b},{c}] can't increase by wavefrontflowab {wavefrontflowab} without exceeding nextlinkcapacity {nextlinkcapacity}, so qijatnodei[{b},{c}] is now {nextlinkcapacity} veh/hr.")
            print(f"qijatnodei[{b},{c}] can't increase by additionalFlowabc {additionalFlowabc} without exceeding nextlinkcapacity {nextlinkcapacity}, so qijatnodei[{b},{c}] is now {nextlinkcapacity} veh/hr.")
          elif qijatnodei[b,c] == nextlinkcapacity:
#            print(f"WARNING: qijatnodei[{b},{c}] is already at nextlinkcapacity {nextlinkcapacity}, so can't increase by wavefrontflowab {wavefrontflowab} or indeed at all: qijatnodei[{b},{c}] is UNCHANGED from {nextlinkcapacity} veh/hr.")
            print(f"WARNING: qijatnodei[{b},{c}] is already at nextlinkcapacity {nextlinkcapacity}, so can't increase by additionalFlowabc {additionalFlowabc} or indeed at all: qijatnodei[{b},{c}] is UNCHANGED from {nextlinkcapacity} veh/hr.")
          densityijatnodei[b,c] = densityonfreeflowbranchoftriangularfundamentaldiagram( qijatnodei[b,c], frozenset({b,c}), argsflowfunc )
          print(f"densityijatnodei[{b},{c}] is now {densityijatnodei[b,c]} veh/km.")
#          copyOfWavefrontPropagatingOnNextArcbc['flow'] = qijatnodei[b,c]
#          copyOfWavefrontPropagatingOnNextArcbc['density'] = densityijatnodei[b,c]
          assert copyOfWavefrontPropagatingOnNextArcbc['flow'] == qijatnodei[b,c]
          assert copyOfWavefrontPropagatingOnNextArcbc['density'] == densityijatnodei[b,c]

          try:
#            assert qijatnodei[b,c] / densityijatnodei[b,c] == argsflowfunc['lK'][frozenset({b,c})]
            assert abs(qijatnodei[b,c] / densityijatnodei[b,c] - argsflowfunc['lK'][frozenset({b,c})]) <= EPS
          except:
            print(f"qijatnodei[{b},{c}] is {qijatnodei[b,c]}, densityijatnodei[{b},{c}] is {densityijatnodei[b,c]}, argsflowfunc['lK'][frozenset({{{b},{c}}})] is {argsflowfunc['lK'][frozenset({b,c})]}")
            print(f"qijatnodei[{b},{c}] / densityijatnodei[{b},{c}] - argsflowfunc['lK'][frozenset({{{b},{c}}})] is {qijatnodei[b,c] / densityijatnodei[b,c] - argsflowfunc['lK'][frozenset({b,c})]}")
            raise
#          try:
#            assert abs(copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density']) > EPS  # the densities ahead of and behind the downstream-propagating wavefront on link (b,c), having just passed through node b, should not be equal; if they are, and the flows ahead of and behind are also equal (which I expect they will always be as a consequence), then copyOfWavefrontPropagatingOnNextArcbc has no effect and can be removed
#          except:
#            print(f"copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] is {copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront']}, copyOfWavefrontPropagatingOnNextArcbc['density']) is {copyOfWavefrontPropagatingOnNextArcbc['density']}")
#            print(f"copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] is {copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront']}, copyOfWavefrontPropagatingOnNextArcbc['flow']) is {copyOfWavefrontPropagatingOnNextArcbc['flow']}")
#            assert abs(copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) <= EPS  # TODO: I expect this to be true, i.e. that when, having just passed through node b, the densities ahead of and behind the downstream-propagating wavefront on link (b,c) are equal, the flows ahead of and behind it are also equal
#            print(f"The densities ahead of behind the downstream-propagating wavefront copyOfWavefrontPropagatingOnNextArcbc are equal, as are the flows, so wavefront {copyOfWavefrontPropagatingOnNextArcbc} is removed from propagatingwavefronts.")
#            propagatingwavefronts.remove( copyOfWavefrontPropagatingOnNextArcbc )
#            print(f"propagatingwavefronts is now {propagatingwavefronts}")
#            raise
          if abs(copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density']) <= EPS:  # the densities ahead of and behind the downstream-propagating wavefront on link (b,c), having just passed through node b, should not be equal; if they are, and the flows ahead of and behind are also equal (which I expect they will always be as a consequence), then copyOfWavefrontPropagatingOnNextArcbc has no effect and can be removed
            print(f"copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] is {copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront']}, copyOfWavefrontPropagatingOnNextArcbc['density']) is {copyOfWavefrontPropagatingOnNextArcbc['density']}")
            print(f"copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] is {copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront']}, copyOfWavefrontPropagatingOnNextArcbc['flow']) is {copyOfWavefrontPropagatingOnNextArcbc['flow']}")
            assert abs(copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) <= EPS  # TODO: I expect this to be true, i.e. that when, having just passed through node b, the densities ahead of and behind the downstream-propagating wavefront on link (b,c) are equal, the flows ahead of and behind it are also equal
            print(f"The densities ahead of behind the downstream-propagating wavefront copyOfWavefrontPropagatingOnNextArcbc are equal, as are the flows, so wavefront {copyOfWavefrontPropagatingOnNextArcbc} has no effect and is removed from propagatingwavefronts.")
            propagatingwavefronts.remove( copyOfWavefrontPropagatingOnNextArcbc )
            print(f"propagatingwavefronts is now {propagatingwavefronts}")
          else:
            try:
#              assert (copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) / (copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density']) == copyOfWavefrontPropagatingOnNextArcbc['propagationspeed']  # the propagation-speed of the downstream-propagating wavefront, having just passed through node b, should now equal the free speed on link (b,c)
              assert abs( (copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) / (copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density']) - copyOfWavefrontPropagatingOnNextArcbc['propagationspeed'] ) <= EPS  # the propagation-speed of the downstream-propagating wavefront, having just passed through node b, should now equal the free speed on link (b,c)
            except:
              print(f"(copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) / (copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density']) - copyOfWavefrontPropagatingOnNextArcbc['propagationspeed'] is {(copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['flow']) / (copyOfWavefrontPropagatingOnNextArcbc['densityaheadofwavefront'] - copyOfWavefrontPropagatingOnNextArcbc['density']) - copyOfWavefrontPropagatingOnNextArcbc['propagationspeed']}")
              raise
            assert copyOfWavefrontPropagatingOnNextArcbc['propagationspeed'] == argsflowfunc['lK'][frozenset({b,c})]
            print(f"New downstream-propagating wavefront is {copyOfWavefrontPropagatingOnNextArcbc}")
#            assert congestedflowab < upstreamfreeflowab
            assert copyOfWavefrontPropagatingOnNextArcbc['flowaheadofwavefront'] < copyOfWavefrontPropagatingOnNextArcbc['flow']
            propagatingwavefronts.append( copyOfWavefrontPropagatingOnNextArcbc )
#            print(f"After appending new downstream-propagating wavefront, propagatingwavefronts is {propagatingwavefronts}")

        if sumOfSubflowsPropagatingDownstreamFromNodeb < wavefrontflowab:  # an upstream-propagating shockfront is created
          print(f"The wavefront's flow of {wavefrontflowab} veh/hr on arc ({a},{b}) exceeds (by {wavefrontflowab-sumOfSubflowsPropagatingDownstreamFromNodeb}) the total additional flow of {sumOfSubflowsPropagatingDownstreamFromNodeb} that can go from node {b} down succeeding links, so a bottleneck results: an upstream-propagating shock-front is created in ({a},{b}), and the wavefront's flow will be reduced from {wavefrontflowab} to {sumOfSubflowsPropagatingDownstreamFromNodeb} veh/hr.")
#          sumOfSubflowsPropagatingDownstreamFromNodeb += congestedflowab
          congestedflowab = sumOfSubflowsPropagatingDownstreamFromNodeb
          print(f"Total congested flow[{a}, {b}] is {congestedflowab} veh/hr")
          print(f"As congestion-zone propagates upstream, flow in link [{a}, {b}] will decrease from {wavefrontflowab} to {congestedflowab} veh/hr")
          # Now that wavefront has reached the end of its current link, update qijatnodej, densityijatnodej, and speedijatnodej at the link's departing-end:
#          assert abs(qijatnodej[a, b] - congestedflowab) > EPS  # not generally true
          qijatnodej[a, b] = congestedflowab
#          print(f"edge is {edge}, propagatingwavefronts is {propagatingwavefronts}")
          assert edge == frozenset({a, b})
          congesteddensityab = densityoncongestedbranchoftriangularfundamentaldiagram( congestedflowab, edge, argsflowfunc )
          densityijatnodej[edge] = congesteddensityab
#          print(f"densityijatnodej is now {densityijatnodej} veh/km")
          speedijatnodej[edge] = congestedflowab / congesteddensityab
          flowstateijatnodej[a, b] = 'congestedflowstate'
          print(f"qijatnodej[{a}, {b}] is now {qijatnodej[a, b]} veh/hr, densityijatnodej[{edge}] is {densityijatnodej[edge]}, speedijatnodej[{edge}] is {speedijatnodej[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
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
#          print(f"Density in link {edge} is now {densityijatnodej[edge]}, speed is now {speedijatnodej[edge]}.")
          timetillupstreamnodeonlink = -lengthofcurrentlink / shockFrontPropagationSpeed
          assert timetillupstreamnodeonlink > 0.0
          shockfrontdict = {'time': simulationclocktime + nexttimestepduration, 'currentlink': frozenset({a,b}), 'sourcenodeonlink': b, 'destinationnodeonlink': a, 'distancereacheddownlink': 0.0, 'propagationspeed': -shockFrontPropagationSpeed, 'flow': congestedflowab, 'density': congesteddensityab, 'timetilldestinationnodeonlink': timetillupstreamnodeonlink, 'nodesonpathtoanexit': 'upstream-propagating shock-front, so not applicable', 'flowaheadofwavefront': upstreamfreeflowab, 'densityaheadofwavefront': upstreamfreedensityab}
          print(f"New upstream-propagating shock-front is {shockfrontdict}")
          assert congestedflowab < upstreamfreeflowab
          propagatingwavefronts.append( shockfrontdict )
#          print(f"After appending new upstream-propagating shock-front, propagatingwavefronts is {propagatingwavefronts}")

        else:  # there is no bottleneck at node b, so no upstream-propagating shockfront is created
          # Now that wavefront has reached the end of its current link, update qijatnodej, densityijatnodej, and speedijatnodej at the link's departing-end:
          qijatnodej[a,b] = wavefrontflowab  # flow at the end of the current link
          densityijatnodej[edge] = nextwavefronttopropagate['density']
          speedijatnodej[edge] = nextwavefronttopropagate['propagationspeed']
          flowstateijatnodej[a,b] = 'freeflowstate'
          print(f"qijatnodej[{a}, {b}] is now {qijatnodej[a, b]} veh/hr, densityijatnodej[{edge}] is {densityijatnodej[edge]}, speedijatnodej[{edge}] is {speedijatnodej[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
#          qijatnodei[b,c] += wavefrontflowab  # TODO: correct to comment out?
#          print(f"After increasing by {wavefrontflowab}, qijatnodei[{b},{c}] is now {qijatnodei[b,c]} veh/hr.")

        print(f"Wavefront {nextwavefronttopropagate} has been copied on each succeeding arc that begins at node {b} and heads towards an assigned exit-node of nextwavefronttopropagate, and is now removed from propagatingwavefronts.")
        propagatingwavefronts.remove( nextwavefronttopropagate )
        print(f"propagatingwavefronts is now {propagatingwavefronts}")

      else:  # either nextwavefronttopropagate is upstream-propagating, or there are no more untraversed links remaining on this shortest path (i.e. numnodesremainingtobetraversedonpath <= 0)
        if nextwavefronttopropagate['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front':  # the wavefront is upstream-propagating; does not necessarily exit the network at its next node
          parentnodesofb = parentsOnShortestPathsFromInjectionNodes[b]
          print(f"Node {b}'s parent-nodes on shortest paths from injection-nodes are {parentnodesofb}")
#          print(f"qijatnodei[{b},{a}] is currently {qijatnodei[b,a]} veh/hr, densityijatnodei[{frozenset({b,a})}] is currently {densityijatnodei[frozenset({b,a})]} veh/km.")
          originalinjectionflowatb = qijatnodei[b,a]
          print(f"Before shock-front reaches node {b}, qijatnodei[{b},{a}] is {qijatnodei[b,a]} veh/hr.")
          qijatnodei[b,a] = nextwavefronttopropagate['flow']
          densityijatnodei[frozenset({b,a})] = nextwavefronttopropagate['density']
          print(f"qijatnodei[{b},{a}] is now {qijatnodei[b,a]} veh/hr, densityijatnodei[{frozenset({b,a})}] is now {densityijatnodei[frozenset({b,a})]} veh/km.")
          if not parentnodesofb:  # b has no parent-nodes on the shortest paths, so nextwavefronttopropagate exits the network at its next node:
            if qijatnodei[b,a] < originalinjectionflowatb:
              congestionAlertMessage = f"ALERT: at injection-node {pointNodeIDbyIndex[b]}, congested injection-flow {qijatnodei[b,a]} < original injection-flow {originalinjectionflowatb} veh/hr: there is a shock where traffic-flow enters the network there."
              print(congestionAlertMessage)
              congestionAlertMessages.append( congestionAlertMessage )
              assert pointNodeIDbyIndex[b] in injectionnodes
            nextwavefronttopropagate['timetilldestinationnodeonlink'] = 0.0
            nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['propagationspeed'] * nexttimestepduration
#            assert nextwavefronttopropagate['distancereacheddownlink'] == argsflowfunc['linklength'][ frozenset({a, b}) ]
            assert abs(nextwavefronttopropagate['distancereacheddownlink'] - lengthofcurrentlink) <= EPS
            print(f"As it exits the network at exit-node {b}, wavefront {nextwavefronttopropagate} is removed from propagatingwavefronts.")
            propagatingwavefronts.remove( nextwavefronttopropagate )
          else:  # update upstream-propagating shockwave as it travels up node b's parent-links
            timeThatShockFrontReachedNodeb = simulationclocktime + nexttimestepduration
            print(f"timeThatShockFrontReachedNode {b} is {timeThatShockFrontReachedNodeb} hours.")
#            assert len(parentsOnShortestPathsFromInjectionNodes[b]) == 1  # TODO: generalise to the case in which node b has more than one parent-link on the shortest paths from injection-nodes
#            p0 = parentsOnShortestPathsFromInjectionNodes[b][0]  # p0 is the single parent-node of b
            print(f"Before propagateSEM4congestionUpstreamToParents(), nextwavefronttopropagate is {nextwavefronttopropagate}")
#            for p in parentnodesofb:
##              flowenteringparent = flowEnteringNodeSEM3( p, i, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, flowijfunc, flowfunc, injectionnodes )
#              try:
##                assert b == shortestPathToAnyExit[p][1]  # b is the next node in the shortest path from p to an exit-node  # TODO: probably no longer true in the general case in which traffic at an injection-node is assigned to more than one exit-node
#                assert b == shortestPathToExit[p,e][1]  # b is the next node in the shortest path from p to exit-node e  # TODO: we don't know which node e is; could find out, but perhaps not worth the trouble
#              except:
##                print(f"Node b {b} does not equal shortestPathToAnyExit[{p}][{1}] {shortestPathToAnyExit[p][1]}.")
#                print(f"Node b {b} does not equal shortestPathToExit[{p},{e}][{1}] {shortestPathToAnyExit[p,e][1]}.")
#                raise
##              linkcapacity = argsflowfunc['linkcapacity'][frozenset({p,i})]
##              if flowenteringparent > linkcapacity:
##                print(f"WARNING: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {i}) which has a capacity of {linkcapacity:.3f}: the link will have flow at this capacity, while a congestion-zone will propagate up the upstream link(s); the congestion-zone will have density to the right of the hump in its link's flow-function, and correspondingly low speed.")
##              flowpb = qijatnodej[p, b]
##              densitypb = densityijatnodej[frozenset({p, b})]
##              speedpb = flowpb / densitypb
##            propagateSEM4congestionUpstreamToParents( nextwavefronttopropagate, b, p, flowpb, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes )
            congestionAlertMessage = propagateSEM4congestionUpstreamToParents( nextwavefronttopropagate, b, timeThatShockFrontReachedNodeb, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, densityijatnodej, speedijatnodej, flowstateijatnodej, argsflowfunc, injectionnodes, b0, indexamonginjectionorexitnodes, propagatingwavefronts, nexttimestepduration )
            congestionAlertMessages.append( congestionAlertMessage )

        else:  # the wavefront is downstream-propagating, and exits the network at its next node
          print(f"b is {b}, pointNodeIDbyIndex[{b}] is {pointNodeIDbyIndex[b]}, exitnodes is {exitnodes}.")
#          assert b in exitnodes
          assert pointNodeIDbyIndex[b] in exitnodes
          nextwavefronttopropagate['timetilldestinationnodeonlink'] = 0.0
          nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['propagationspeed'] * nexttimestepduration
          print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']}, lengthofcurrentlink is {lengthofcurrentlink}")
#          assert nextwavefronttopropagate['distancereacheddownlink'] == argsflowfunc['linklength'][ frozenset({a, b}) ]
          assert abs(nextwavefronttopropagate['distancereacheddownlink'] - lengthofcurrentlink) <= EPS
          qijatnodej[a, b] = nextwavefronttopropagate['flow']  # wavefront has reached the end of its current link, so update qijatnodej, densityijatnodej, and speedijatnodej
          densityijatnodej[edge] = nextwavefronttopropagate['density']
          speedijatnodej[edge] = nextwavefronttopropagate['propagationspeed']
          flowstateijatnodej[a, b] = 'freeflowstate'
          print(f"qijatnodej[{a}, {b}] is now {qijatnodej[a, b]} veh/hr, densityijatnodej[{edge}] is {densityijatnodej[edge]}, speedijatnodej[{edge}] is {speedijatnodej[edge]}, flowstateijatnodej[{a}, {b}] is {flowstateijatnodej[a, b]}.")
          print(f"As it exits the network at exit-node {b}, wavefront {nextwavefronttopropagate} is removed from propagatingwavefronts.")
          propagatingwavefronts.remove( nextwavefronttopropagate )
#          assert nextwavefronttopropagate['nodesonpathtoanexit'][-1] == b
          assert len(nextwavefronttopropagate['subflowsonpathstoexits'].keys()) == 1 and b in nextwavefronttopropagate['subflowsonpathstoexits'] and nextwavefronttopropagate['subflowsonpathstoexits'][b]['shortestpath'][-1] == b
          outflowsByNodeID[pointNodeIDbyIndex[b]] = -nextwavefronttopropagate['flow']

    elif pointReachedByNextPropagatingWaveFront == 'catch up to slower-propagating shockwave on same link':
      assert nextwavefronttopropagate['flow'] < wavefrontIntersected['flow']  # TODO: perhaps not generally true, but I expect it is: the upstream-propagating shockwave should be trailing (have behind it, downstream) a flow that's lower (more congested) than the flow downstream of the upstream-propagating shockwave it's about to intersect
      assert nextwavefronttopropagate['flowaheadofwavefront'] == wavefrontIntersected['flow']  # TODO: this line and the next one are perhaps not generally true, but I expect they are: before the intersection, there should be no other wavefront between nextwavefronttopropagate and wavefrontIntersected
      assert nextwavefronttopropagate['densityaheadofwavefront'] == wavefrontIntersected['density']
#      nextwavefronttopropagate['timetilldestinationnodeonlink'] -= nexttimestepduration
      nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['propagationspeed'] * nexttimestepduration
      wavefrontIntersected['timetilldestinationnodeonlink'] -= nexttimestepduration
      wavefrontIntersected['distancereacheddownlink'] += wavefrontIntersected['propagationspeed'] * nexttimestepduration
      print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']} km, wavefrontIntersected['distancereacheddownlink'] is {wavefrontIntersected['distancereacheddownlink']} km")
#      assert nextwavefronttopropagate['distancereacheddownlink'] == wavefrontIntersected['distancereacheddownlink']
      assert abs(nextwavefronttopropagate['distancereacheddownlink'] - wavefrontIntersected['distancereacheddownlink']) <= EPS
#      print(f"nextwavefronttopropagate['timetilldestinationnodeonlink'] is {nextwavefronttopropagate['timetilldestinationnodeonlink']}, wavefrontIntersected['timetilldestinationnodeonlink'] is {wavefrontIntersected['timetilldestinationnodeonlink']}")
#      print(f"lengthofcurrentlink is {lengthofcurrentlink} km")
      upstreamflow = wavefrontIntersected['flowaheadofwavefront']
      upstreamdensity = wavefrontIntersected['densityaheadofwavefront']
      congestedflow = nextwavefronttopropagate['flow']
      congesteddensity = nextwavefronttopropagate['density']
      print(f"Calculating new propagation-speed of shock-front: upstreamflow is {upstreamflow}, upstreamdensity is {upstreamdensity}, congestedflow is {congestedflow}, congesteddensity is {congesteddensity}")
      assert nextwavefronttopropagate['flow'] <= wavefrontIntersected['flow']  # TODO: a working assumption I've made that might well not be generally true
      shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamflow, upstreamdensity, congestedflow, congesteddensity )
      if shockFrontPropagationSpeed < 0.0:
        upordownstreamText = 'travels upstream'
      elif shockFrontPropagationSpeed > 0.0:
        upordownstreamText = 'travels downstream'
      else:
        upordownstreamText = 'is stationary'
      print(f"New shockFrontPropagationSpeed is {shockFrontPropagationSpeed} km/hr: the front {upordownstreamText}")
      assert shockFrontPropagationSpeed < 0.0
      nextwavefronttopropagate['propagationspeed'] = -shockFrontPropagationSpeed
      nextwavefronttopropagate['flowaheadofwavefront'] = wavefrontIntersected['flowaheadofwavefront']
      nextwavefronttopropagate['densityaheadofwavefront'] = wavefrontIntersected['densityaheadofwavefront']
      nextwavefronttopropagate['timetilldestinationnodeonlink'] = (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed']  # TODO: probably need this line: the assertion two lines below will fail without it
      try:
        assert nextwavefronttopropagate['timetilldestinationnodeonlink'] == (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed']
      except:
        print(f"nextwavefronttopropagate['timetilldestinationnodeonlink'] - (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed'] is {nextwavefronttopropagate['timetilldestinationnodeonlink'] - (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed']}")
        raise

      print(f"As it's caught up from behind by upstream-propagating shock-front with propagation-speed {nextwavefronttopropagate['propagationspeed']}, wavefront {wavefrontIntersected} is removed from propagatingwavefronts.")
      propagatingwavefronts.remove( wavefrontIntersected )  # TODO: if nextwavefronttopropagate is trailing (has behind it, downstream) higher flow than does wavefrontIntersected, then remove nextwavefronttopropagate instead
      print(f"propagatingwavefronts is now {propagatingwavefronts}")

    elif pointReachedByNextPropagatingWaveFront == 'intersect with wavefront propagating in opposite direction (downstream) on same link':
      assert nextwavefronttopropagate['flow'] < wavefrontIntersected['flow']  # TODO: perhaps not generally true, but I expect it is: the upstream-propagating shockwave should be trailing (have behind it, downstream) a flow that's lower (because congested) than the flow upstream of the downstream-propagating wavefront it's about to intersect

      try:  # commented out for now (3rd May 2021), but this assertion fails at time-step 109, scenario cmr_1s1d1r, inflowsByNodeID = {3: 1000, 4:1500, 5:300, 6:700}
        assert nextwavefronttopropagate['flowaheadofwavefront'] == wavefrontIntersected['flowaheadofwavefront']  # FIXME: this assertion and the next one below are perhaps not generally true, but I expect they are: before the intersection, there should be no other wavefront between nextwavefronttopropagate and wavefrontIntersected
      except:
        print(f"nextwavefronttopropagate['flowaheadofwavefront'] is {nextwavefronttopropagate['flowaheadofwavefront']}, wavefrontIntersected['flowaheadofwavefront'] is {wavefrontIntersected['flowaheadofwavefront']}")
        raise
      try:  # commented out for now (3rd May 2021), but this assertion fails at time-step 109, scenario cmr_1s1d1r, inflowsByNodeID = {3: 1000, 4:1500, 5:300, 6:700}
        assert nextwavefronttopropagate['densityaheadofwavefront'] == wavefrontIntersected['densityaheadofwavefront']
      except:
        print(f"nextwavefronttopropagate['densityaheadofwavefront'] is {nextwavefronttopropagate['densityaheadofwavefront']}, wavefrontIntersected['densityaheadofwavefront'] is {wavefrontIntersected['densityaheadofwavefront']}")
        raise

#      nextwavefronttopropagate['timetilldestinationnodeonlink'] -= nexttimestepduration
      nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['propagationspeed'] * nexttimestepduration
      wavefrontIntersected['timetilldestinationnodeonlink'] -= nexttimestepduration
      wavefrontIntersected['distancereacheddownlink'] += wavefrontIntersected['propagationspeed'] * nexttimestepduration
      print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']} km, wavefrontIntersected['distancereacheddownlink'] is {wavefrontIntersected['distancereacheddownlink']} km")
#      print(f"lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink'] - wavefrontIntersected['distancereacheddownlink'] is {lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink'] - wavefrontIntersected['distancereacheddownlink']}")
      assert abs(lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink'] - wavefrontIntersected['distancereacheddownlink']) <= EPS
      upstreamflow = wavefrontIntersected['flow']
      upstreamdensity = wavefrontIntersected['density']
      congestedflow = nextwavefronttopropagate['flow']
      congesteddensity = nextwavefronttopropagate['density']
      print(f"Calculating new propagation-speed of shock-front: upstreamflow is {upstreamflow}, upstreamdensity is {upstreamdensity}, congestedflow is {congestedflow}, congesteddensity is {congesteddensity}")
#      assert nextwavefronttopropagate['flow'] < wavefrontIntersected['flow']  # TODO: a working assumption I've made that might well not be generally true; this assertion already appears as the first line of this code-block
      shockFrontPropagationSpeed = propagationSpeedOfShockFront( upstreamflow, upstreamdensity, congestedflow, congesteddensity )
      if shockFrontPropagationSpeed < 0.0:
        upordownstreamText = 'travels upstream'
      elif shockFrontPropagationSpeed > 0.0:
        upordownstreamText = 'travels downstream'
      else:
        upordownstreamText = 'is stationary'
      print(f"New shockFrontPropagationSpeed is {shockFrontPropagationSpeed} km/hr: the front {upordownstreamText}")
      assert shockFrontPropagationSpeed < 0.0
      nextwavefronttopropagate['propagationspeed'] = -shockFrontPropagationSpeed
      nextwavefronttopropagate['flowaheadofwavefront'] = wavefrontIntersected['flow']
      nextwavefronttopropagate['densityaheadofwavefront'] = wavefrontIntersected['density']
      nextwavefronttopropagate['timetilldestinationnodeonlink'] = (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed']
      try:
        assert nextwavefronttopropagate['timetilldestinationnodeonlink'] == (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed']
      except:
        print(f"nextwavefronttopropagate['timetilldestinationnodeonlink'] - (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed'] is {nextwavefronttopropagate['timetilldestinationnodeonlink'] - (lengthofcurrentlink - nextwavefronttopropagate['distancereacheddownlink']) / nextwavefronttopropagate['propagationspeed']}")
        raise
      print(f"As it's intersected by upstream-propagating shock-front with propagation-speed {nextwavefronttopropagate['propagationspeed']}, downstream-propagating wavefront {wavefrontIntersected} is removed from propagatingwavefronts.")
      propagatingwavefronts.remove( wavefrontIntersected )  # TODO: if nextwavefronttopropagate is trailing (has behind it, downstream) higher flow than does wavefrontIntersected, then remove nextwavefronttopropagate instead
      print(f"propagatingwavefronts is now {propagatingwavefronts}")

    elif pointReachedByNextPropagatingWaveFront[:60] == 'travel until the simulation clock-time reaches simulDuration':
      nextwavefronttopropagate['timetilldestinationnodeonlink'] -= nexttimestepduration
      nextwavefronttopropagate['distancereacheddownlink'] += nextwavefronttopropagate['propagationspeed'] * nexttimestepduration
      print(f"nextwavefronttopropagate['distancereacheddownlink'] is {nextwavefronttopropagate['distancereacheddownlink']} km")
      print(f"As its propagation ceases at the end of the simulation, wavefront {nextwavefronttopropagate} is removed from propagatingwavefronts.")
      propagatingwavefronts.remove( nextwavefronttopropagate )

    if nextwavefronttopropagate in propagatingwavefronts:  # if this wavefront hasn't just been removed from 'propagatingwavefronts'
      print(f"Updated propagating wavefront is {nextwavefronttopropagate}")
##    print(f"After updating downstream-propagating wavefront, propagatingwavefronts is {propagatingwavefronts}")
#    print(f"After updating downstream-propagating wavefront, len(propagatingwavefronts) is {len(propagatingwavefronts)}")

    for (i,j) in qijatnodej:
      try:
        assert qijatnodei[i,j] <= totalinflows
      except:
        print(f"qijatnodei[{i},{j}] is {qijatnodei[i,j]}")
        raise
      try:
        assert qijatnodej[i,j] <= totalinflows
      except:
        print(f"qijatnodej[{i},{j}] is {qijatnodej[i,j]}")
        raise

    if nextwavefronttopropagate['nodesonpathtoanexit'][:32] != 'upstream-propagating shock-front':
#      print(f"nextwavefronttopropagate['flow'] is {nextwavefronttopropagate['flow']}, nextwavefronttopropagate['density'] is {nextwavefronttopropagate['density']}, nextwavefronttopropagate['propagationspeed'] is {nextwavefronttopropagate['propagationspeed']}.")
#      print(f"nextwavefronttopropagate['flow'] is {nextwavefronttopropagate['flow']}, nextwavefronttopropagate['density']*nextwavefronttopropagate['propagationspeed'] is {nextwavefronttopropagate['density']*nextwavefronttopropagate['propagationspeed']}.")
#      assert nextwavefronttopropagate['flow'] == nextwavefronttopropagate['density'] * nextwavefronttopropagate['propagationspeed']
      assert abs(nextwavefronttopropagate['flow'] - nextwavefronttopropagate['density'] * nextwavefronttopropagate['propagationspeed']) <= EPS
#    else:
#      assert nextwavefronttopropagate['flow'] == nextwavefronttopropagate['density'] * nextwavefronttopropagate['propagationspeed']  # TODO: maybe replace this with an assertion that is correct for the case of an upstream-propagating wavefront
    print(f"After updating propagating wavefront, len(propagatingwavefronts) is {len(propagatingwavefronts)}")

#    otherpropagatingwavefronts = [wv for wv in propagatingwavefronts if wv['timetilldestinationnodeonlink'] > nexttimestepduration]
    otherpropagatingwavefronts = [wv for wv in propagatingwavefronts if wv['time'] == simulationclocktime]
    print(f"otherpropagatingwavefronts, i.e. with 'time' {simulationclocktime}, is {otherpropagatingwavefronts}")
#    assert len(otherpropagatingwavefronts) == len(propagatingwavefronts) - 1  # not generally true
    for wavefront in otherpropagatingwavefronts:
#      if wavefront['propagationspeed'] <= 0.0:
#        continue
      wavefront['time'] += nexttimestepduration
      wavefront['distancereacheddownlink'] += wavefront['propagationspeed'] * nexttimestepduration
      wavefront['timetilldestinationnodeonlink'] -= nexttimestepduration
      if wavefront['nodesonpathtoanexit'][:32] == 'upstream-propagating shock-front':
        direction = 'upstream'
      else:
        direction = 'downstream'
      print(f"Updated {direction}-propagating otherwavefront is {wavefront}")
  
    verifyOnEveryLinkThatPropagatingFlowMatchesqijValues(propagatingwavefronts, argsflowfunc, qijatnodei, qijatnodej)
    simulationclocktime += nexttimestepduration
    timestepcount += 1
    printSimulationClockTime(timestepcount, simulationclocktime)
#    sys.exit()

  if simulationclocktime < simulDuration:  # run final time-step, ending at the specified duration of the simulation
    assert propagatingwavefronts == []
    finaltimestepduration = simulDuration - simulationclocktime  # the final time-step's duration, by the start of which every propagating wavefront has departed the network: this implies that during the final time-step, the flow in (the entire length of) every directed link (i, j) adds to cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j]
    print(f"finaltimestepduration is {finaltimestepduration} hours")
    for (i, j) in qijatnodej:
      print(f"Current flow on entire length of link ({i}, {j}) is qijatnodej[{i}, {j}]={qijatnodej[i, j]}.")
      numVehiclesDepartingLinkDuringThisTimestep = qijatnodej[i,j] * finaltimestepduration
      print(f"numVehiclesDepartingLinkDuringThisTimestep on link ({i}, {j}) is {numVehiclesDepartingLinkDuringThisTimestep}")
      cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j] += numVehiclesDepartingLinkDuringThisTimestep
      print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation[{i}, {j}] is now {cumulNumVehiclesToHaveDepartedLinkDuringSimulation[i, j]}")
    simulationclocktime += finaltimestepduration
    timestepcount += 1
    printSimulationClockTime(timestepcount, simulationclocktime)

  return (cumulNumVehiclesToHaveDepartedLinkDuringSimulation, outflowsByNodeID, congestionAlertMessages)


def findShortestPathsSEM3or4( N, argsflowfunc, injectionnodes, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID):
  starttime = time.time()
  print("\nFinding shortest-path routes from all injection-nodes ...")
  vertices = range(N)
  predictedLinkTraversalTime = {}  # traversal-time as might be predicted by a driver
  arcs = []
  lK = argsflowfunc['lK']  # lK stores the 'free speed' of each road-link
  oneway = argsflowfunc['oneway']  # oneway stores 'True' if the road-link is one-way, 'False' otherwise
  linklength = argsflowfunc['linklength']
#  for edge in lK.keys():
  for arc in oneway.keys():
#    (a, b) = list(edge)
    (a, b) = arc
    edge = frozenset({a,b})
#    if oneway[arc]:
#      predictedLinkTraversalTime[a,b] = linklength[edge] / lK[edge]
#      print("predictedLinkTraversalTime[%d, %d] is %s" % (a, b, predictedLinkTraversalTime[a, b]))
#    else:
#      predictedLinkTraversalTime[a,b] = predictedLinkTraversalTime[b,a] = linklength[edge] / lK[edge]
#      print("predictedLinkTraversalTime[%d, %d] is %s" % (a, b, predictedLinkTraversalTime[a, b]))
#      print("predictedLinkTraversalTime[%d, %d] is %s" % (b, a, predictedLinkTraversalTime[b, a]))
    predictedLinkTraversalTime[a,b] = linklength[edge] / lK[edge]
    print("predictedLinkTraversalTime[%d, %d] is %s" % (a, b, predictedLinkTraversalTime[a, b]))
    arcs.append( (a,b) )
    if not oneway[arc]:
      predictedLinkTraversalTime[b,a] = predictedLinkTraversalTime[a,b]
      print("predictedLinkTraversalTime[%d, %d] is %s" % (b, a, predictedLinkTraversalTime[b, a]))
      arcs.append( (b,a) )
  exception = None  # no exception, such as numpy.linalg.LinAlgError, is raised by the global-optimisation solver as it's not being employed
  (pathdistance, nextNodeInShortestPath) = FloydWarshall( arcs, vertices, predictedLinkTraversalTime )
  print("pathdistance is", pathdistance)
  print("nextNodeInShortestPath is", nextNodeInShortestPath)
  endtime_runShortestPathsAlgorithm = time.time()
  timetorunShortestPathsAlgorithm = endtime_runShortestPathsAlgorithm - starttime

  parentsOnShortestPathsFromInjectionNodes = {}  # for each node, its immediate ancestors on all those shortest paths from injection-nodes that pass through that node
#  for j in range(N):
  for nodeID2 in injectionnodes+list(exitnodes):
    j = pointcoordsandindexbyID[nodeID2]['index'] 
    parentsOnShortestPathsFromInjectionNodes[j] = set()

#  shortestDistanceToAnyExit = np.full( shape=len(injectionnodes), fill_value=9e99  )  # initialise to "infinity" each shortest distance from a non-leaf to an exit-node
  shortestDistanceToAnyExit = {}  # for each non-leaf node, its shortest distance to an exit-node
#  closestExit = {}  # for each non-leaf node, its closest exit-node
  shortestPathToAnyExit = {}  # for each non-leaf node, the path to its closest exit-node
  shortestPathToExit = {}  # for each injection-node i and exit-node e, the shortest path from i to e
##  for i in range(N):  # for each node construct the shortest path to an exit, its length, and the node-ID of that closest exit
  for nodeID1 in injectionnodes:  # only for injection-nodes do we need to construct the shortest path to an exit, its length, and the node-ID of that closest exit
    i = pointcoordsandindexbyID[nodeID1]['index'] 
#    for j in range(N):
    for nodeID2 in injectionnodes+list(exitnodes):
      j = pointcoordsandindexbyID[nodeID2]['index'] 
      if j == i: continue
      (shortestpath, shortestpathlength) = constructShortestPath(i, j, nextNodeInShortestPath, predictedLinkTraversalTime)
      if isinstance(shortestpathlength, float):
        print("Shortest path from %d to %d is %s with length %.4f" % (i, j, shortestpath, shortestpathlength) )
      else:
        assert isinstance(shortestpathlength, str)
        print("Shortest path from %d to %d is %s with length %s" % (i, j, shortestpath, shortestpathlength) )
#      if j in exitnodes and shortestpathlength < shortestDistanceToAnyExit.get(i, 9e99):
      if nodeID2 in exitnodes:
        shortestPathToExit[i, j] = shortestpath
        if shortestpathlength < shortestDistanceToAnyExit.get(i, 9e99):
#          print("Shortest path from %d to exit-node %d is %s with length %.4f" % (i, j, shortestpath, shortestpathlength) )
          shortestDistanceToAnyExit[i] = shortestpathlength
#          closestExit[i] = j
          shortestPathToAnyExit[i] = shortestpath
#      print("shortestpath is", shortestpath)

  qijatnodei = {}  # store initial flows at entering-ends of all arcs, i.e. directed edges, that are in the shortest-path tree
  densityijatnodei = {}  # store initial densities at entering-ends of all arcs, i.e. directed edges, that are in the shortest-path tree
  qijatnodej = {}  # store initial flows at leaving-ends of all arcs, i.e. directed edges, that are in the shortest-path tree
  flowstateijatnodej = {}  # store flow-state, either 'freestate' or 'congestedstate', on all arcs, i.e. directed edges, that are in the shortest-path tree
##  for i in injectionnodes:  # proceed down each injection-node's shortest path to an exit, adding for each node after the first one on the path its parent on the path as one of that node's parents on all shortest paths from injection-nodes
#  directedNeighboursInShortestPathTreeOfNode = {}  # directedNeighboursInShortestPathTreeOfNode[i] is a tuple of all nodes j for which the network has an arc (i,j)
  assignedExitNodesReachedDownstreamOfArc = {}  # assignedExitNodesReachedDownstreamOfArc[i,b,c] is the set of exit-nodes assigned to injection-node i that are downstream of arc (b,c) and can be reached via the shortest paths found; if c is an assigned exit-node for i then this set includes node c
  successorNodeOnShortestPathFromNodeToExitNode = {}  # successorNodeOnShortestPathFromNodeToExitNode[b,e] is the successor-node of b on the shortest path to exit-node e
#  for i in range(N):
#    directedNeighboursInShortestPathTreeOfNode[i] = set()
##    assignedExitNodesReachedDownstreamOfArc[i] = set()
  for nodeID in injectionnodes:  # proceed down each injection-node's shortest path to an exit, adding for each node after the first one on the path its parent on the path as one of that node's parents on all shortest paths from injection-nodes
    i = pointcoordsandindexbyID[nodeID]['index'] 
    if i not in shortestPathToAnyExit:
      print(f"There is no path to an exit from node {i}.")
      continue
#    print(f"shortestPathToAnyExit[{i}] is {shortestPathToAnyExit[i]}")
#    for index, j in enumerate( shortestPathToAnyExit[i] ):
#      if index > 0:
#        parentofj = shortestPathToAnyExit[i][index-1]
#        parentsOnShortestPathsFromInjectionNodes[j].add( parentofj )
#        qijatnodei[parentofj, j] = 0.0
#        densityijatnodei[parentofj, j] = 0.0
#        qijatnodej[parentofj, j] = 0.0
#        flowstateijatnodej[parentofj, j] = 'freeflowstate'
    if nodeID not in subflowsToAssignedExitNodesByInjectionNodeID: continue
    for exitnode in subflowsToAssignedExitNodesByInjectionNodeID[nodeID]:
      e = pointcoordsandindexbyID[exitnode]['index']
#      assignedExitNodesReachedDownstreamOfArc[e].add( e )
      for index, j in enumerate( shortestPathToExit[i, e] ):
        if index > 0:
#          parentofj = shortestPathToAnyExit[i][index-1]
          parentofj = shortestPathToExit[i,e][index-1]
          parentsOnShortestPathsFromInjectionNodes[j].add( parentofj )
          qijatnodei[parentofj, j] = 0.0
          densityijatnodei[parentofj, j] = 0.0
          qijatnodej[parentofj, j] = 0.0
          flowstateijatnodej[parentofj, j] = 'freeflowstate'
#          directedNeighboursInShortestPathTreeOfNode[parentofj].add( j )
          successorNodeOnShortestPathFromNodeToExitNode[parentofj,e] = j
#          assignedExitNodesReachedDownstreamOfArc[parentofj].add( e )
          if (i,parentofj,j) in assignedExitNodesReachedDownstreamOfArc:
            assignedExitNodesReachedDownstreamOfArc[i,parentofj,j].add( e )
          else:
            assignedExitNodesReachedDownstreamOfArc[i,parentofj,j] = {e}
##  for arc in oneway.keys():
##    (a, b) = arc
##    directedNeighboursInShortestPathTreeOfNode[a].add( b )
##    if not oneway[arc]:
##      directedNeighboursInShortestPathTreeOfNode[b].add( a )
#  print(f"directedNeighboursInShortestPathTreeOfNode is {directedNeighboursInShortestPathTreeOfNode}")
  print(f"successorNodeOnShortestPathFromNodeToExitNode is {successorNodeOnShortestPathFromNodeToExitNode}")
  print(f"assignedExitNodesReachedDownstreamOfArc is {assignedExitNodesReachedDownstreamOfArc}")
#  print(f"assignedExitNodesReachedDownstreamOfArc[85] is {assignedExitNodesReachedDownstreamOfArc[85]}")
#  print(f"shortestPathToExit[85,37] is {shortestPathToExit[85,37]}")
#  print(f"shortestPathToExit[85,294] is {shortestPathToExit[85,294]}")
#  sys.exit()
  print("shortestDistanceToAnyExit is", shortestDistanceToAnyExit)
  print(f"shortestPathToExit is {shortestPathToExit}")
#  print("closestExit is", closestExit)
  print("shortestPathToAnyExit is", shortestPathToAnyExit)
  print("parentsOnShortestPathsFromInjectionNodes is", parentsOnShortestPathsFromInjectionNodes)
  print("qijatnodei is", qijatnodei)
#  print("densityijatnodei is", densityijatnodei)
  print("qijatnodej is", qijatnodej)
  print("flowstateijatnodej is", flowstateijatnodej)

  endtime_calculateshortestpaths = time.time()
  timetocalculateshortestpaths = endtime_calculateshortestpaths - endtime_runShortestPathsAlgorithm
#  return (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK, directedNeighboursInShortestPathTreeOfNode, assignedExitNodesReachedDownstreamOfArc)
  return (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK, successorNodeOnShortestPathFromNodeToExitNode, assignedExitNodesReachedDownstreamOfArc)


#def findShortestPathsSEM5(N, argsflowfunc, injectionnodes, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID):
def findShortestPathsSEM5(argsflowfunc, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID):
# Use Dijkstra's algorithm:
  starttime = time.time()
  congestionAlertMessages = []
  print(f"subflowsToAssignedExitNodesByInjectionNodeID is {subflowsToAssignedExitNodesByInjectionNodeID}")
  lengthOfShortestPath = {}  # each key is an arc (u,v)
  shortestPath = {}  # each key is an arc (u,v)
  predictedLinkTraversalTime = {}  # traversal-time as might be predicted by a driver
  arcs = []
  lK = argsflowfunc['lK']  # lK stores the 'free speed' of each road-link
  oneway = argsflowfunc['oneway']  # oneway stores 'True' if the road-link is one-way, 'False' otherwise
  linklength = argsflowfunc['linklength']
  for arc in oneway.keys():
    (a, b) = arc
    edge = frozenset({a,b})
    predictedLinkTraversalTime[a,b] = linklength[edge] / lK[edge]
    print("predictedLinkTraversalTime[%d, %d] is %s" % (a, b, predictedLinkTraversalTime[a, b]))
    arcs.append( (a,b) )
    if not oneway[arc]:
      predictedLinkTraversalTime[b,a] = predictedLinkTraversalTime[a,b]
      print("predictedLinkTraversalTime[%d, %d] is %s" % (b, a, predictedLinkTraversalTime[b, a]))
      arcs.append( (b,a) )
#  print(f"predictedLinkTraversalTime is {predictedLinkTraversalTime}")
  G = nx.DiGraph()  # shortest paths will be found in G, using a networkx method such as 'single_source_dijkstra'
  print("arcs is %s" % arcs)
  G.add_edges_from( arcs )
  print("G.nodes() is %s" % G.nodes())
#  print("G.edges() is %s" % G.edges())
#  qijatnodei = {}  # store initial flows at entering-ends of all arcs, i.e. directed edges, that are in the shortest-path tree
  qij = {}  # store initial flows in all arcs, i.e. directed edges, that are in the shortest-path tree
  for u, v in G.edges():
    G.edges[u,v]['weight'] = predictedLinkTraversalTime[u,v]
#    print(f"G.edges[{u},{v}]['weight'] is {G.edges[u,v]['weight']}")
#  for nodeID1 in injectionnodes:  # only for injection-nodes do we need to construct the shortest path to each exit-node and its length
  for nodeID1 in subflowsToAssignedExitNodesByInjectionNodeID:  # only for the injection-nodes specified as keys of 'subflowsToAssignedExitNodesByInjectionNodeID' do we need to construct the shortest path to each exit-node and its length
    i = pointcoordsandindexbyID[nodeID1]['index']
#    for nodeID2 in exitnodes:
    for nodeID2 in subflowsToAssignedExitNodesByInjectionNodeID[nodeID1]:
      j = pointcoordsandindexbyID[nodeID2]['index']
      print(f"injection-node is i={i}, exit-node is j={j}:")
#      shortestPathDijkstra = nx.dijkstra_path(G, i, j, weight = lambda u, v, predictedLinkTraversalTime: predictedLinkTraversalTime[u,v])
#      print(f"shortestPathDijkstra is {shortestPathDijkstra} km")
      length, path = nx.single_source_dijkstra(G, i, j, weight='weight')
      print(f"  path is {path}, length is {length} km.")
      lengthOfShortestPath[i,j], shortestPath[i,j] = length, path
      for k, u in enumerate(path[:-1]):
        v = path[k+1]
#        print(f"(k is {k})  u is {u}, v is {v}")
        if (u,v) in qij:  # add the subflow that is injected at i and assigned to exit-node j, on every arc in the shortest path from i to j
          qij[u,v] += subflowsToAssignedExitNodesByInjectionNodeID[nodeID1][nodeID2]
        else:
          qij[u,v] = subflowsToAssignedExitNodesByInjectionNodeID[nodeID1][nodeID2]
#      print(f"qij is now {qij}")
#      sys.exit()
      flowsonpath = [qij[u,path[k+1]] for k, u in enumerate(path[:-1])]
      print(f"flowsonpath is {flowsonpath}")
      capacitiesonpath = [argsflowfunc['linkcapacity'][frozenset({u,path[k+1]})] for k, u in enumerate(path[:-1])]
      print(f"capacitiesonpath is {capacitiesonpath}")
  print(f"qij with all subflows summed is {qij}")
  endtime_runShortestPathsAlgorithm = time.time()
  timetorunShortestPathsAlgorithm = endtime_runShortestPathsAlgorithm - starttime
  for (u,v) in qij:  # on each shortest path from an injection-node to an assigned exit-node, compare the flows on the arcs with those arcs' capacities:
    capacity_uv = argsflowfunc['linkcapacity'][frozenset({u,v})]
    if qij[u,v] > capacity_uv:
      congestionAlertMessage = f"ALERT: arc ({u},{v}) has flow-capacity of {capacity_uv} veh/hr, but sum of subflows from injection-nodes to exit-nodes in arc ({u},{v}) = {qij[u,v]} veh/hr, which exceeds capacity."
      print(congestionAlertMessage)
      congestionAlertMessages.append( congestionAlertMessage )
  print(f"A total of {len(congestionAlertMessages)} congestion-alert messages were reported.")
  congestionAlertMessages = []
  for nodeID1 in subflowsToAssignedExitNodesByInjectionNodeID:
    i = pointcoordsandindexbyID[nodeID1]['index']
#    for nodeID2 in exitnodes:
    for nodeID2 in subflowsToAssignedExitNodesByInjectionNodeID[nodeID1]:
      j = pointcoordsandindexbyID[nodeID2]['index']
      print(f"injection-node is i={i}, exit-node is j={j}")
      lengthofshortestpath, path = lengthOfShortestPath[i,j], shortestPath[i,j]
      flowsonpath = [qij[u,path[k+1]] for k, u in enumerate(path[:-1])]
      print(f"flowsonpath is {flowsonpath}")
      capacitiesonpath = [argsflowfunc['linkcapacity'][frozenset({u,path[k+1]})] for k, u in enumerate(path[:-1])]
      print(f"capacitiesonpath is {capacitiesonpath}")
  endtime_calculatecongestionpoints = time.time()
  timeToCalculateCongestionPoints = endtime_calculatecongestionpoints - endtime_runShortestPathsAlgorithm
  print(f"timetorunShortestPathsAlgorithm is {timetorunShortestPathsAlgorithm:.5f} seconds.")
  print(f"timeToCalculateCongestionPoints is {timeToCalculateCongestionPoints:.5f} seconds.")
#  sys.exit()
#  return (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK)


def findMaximumFlowSEM5(argsflowfunc, exitnodes, pointcoordsandindexbyID, pointNodeIDbyIndex, inflowsByNodeID, subflowsToAssignedExitNodesByInjectionNodeID):
# Uses Edmonds-Karp maximum-flow algorithm?
# TODO: maximum-flow method can't handle parallel arcs, so must replace each two-way edge, i.e. pair of arcs in opposite directions, by adding a dummy node; see p. 711 of 'Introduction to Algorithms', third edition.
  starttime = time.time()
  EPS = 1e-6
  FINDMAXIMUMFLOWONSHORTESTPATHARCSONLY = False  #  in calculation of maximum flow, use all arcs of the graph
#  FINDMAXIMUMFLOWONSHORTESTPATHARCSONLY = True  #  in calculation of maximum flow, use only the arcs that appear on shortest paths
  inflowNodesByIndex = [pointcoordsandindexbyID[nodeID]['index'] for nodeID in inflowsByNodeID]
  print(f"inflowNodesByIndex is {inflowNodesByIndex}")
  alertsNodesWhereDemandNotSatisfied = {}
  alertsLinksPossiblyCausingCongestion = {}
  arcs = []
  oneway = argsflowfunc['oneway']  # oneway stores 'True' if the road-link is one-way, 'False' otherwise
  for arc in oneway.keys():
    (a, b) = arc
    arcs.append( (a,b) )
    if not oneway[arc]:
      arcs.append( (b,a) )
#  print("arcs is %s" % arcs)
  Gmaximumflow = nx.DiGraph()  # maximum flow will be calculated using arcs of Gmaximumflow, using a networkx method such as 'single_source_dijkstra'
  if FINDMAXIMUMFLOWONSHORTESTPATHARCSONLY:  # use Dijkstra's algorithm:
    starttime_runShortestPathsAlgorithm = time.time()
    predictedLinkTraversalTime = {}  # traversal-time as might be predicted by a driver
    lK = argsflowfunc['lK']  # lK stores the 'free speed' of each road-link
    linklength = argsflowfunc['linklength']
    for arc in arcs:
      (a, b) = arc
      edge = frozenset({a,b})
      predictedLinkTraversalTime[a,b] = linklength[edge] / lK[edge]
#      print("predictedLinkTraversalTime[%d, %d] is %s" % (a, b, predictedLinkTraversalTime[a, b]))
    G = nx.DiGraph()  # shortest paths will be found in G, using a networkx method such as 'single_source_dijkstra'
    G.add_edges_from( arcs )
    for u, v in G.edges():
      G.edges[u,v]['weight'] = predictedLinkTraversalTime[u,v]
#      print(f"G.edges[{u},{v}]['weight'] is {G.edges[u,v]['weight']}")
    arcsInShortestPaths = []
    totalAssignedSubflowsUsingArc = {}  # for each arc (u,v) that appears in at least one shortest path, the total of all assigned subflows that use (u,v)
    for nodeID1 in subflowsToAssignedExitNodesByInjectionNodeID:  # only for the injection-nodes specified as keys of 'subflowsToAssignedExitNodesByInjectionNodeID' do we need to construct the shortest path to each exit-node and its length
      i = pointcoordsandindexbyID[nodeID1]['index']
#      for nodeID2 in exitnodes:
      for nodeID2 in subflowsToAssignedExitNodesByInjectionNodeID[nodeID1]:
        j = pointcoordsandindexbyID[nodeID2]['index']
        print(f"injection-node is i={i}, exit-node is j={j}:")
#        shortestPathDijkstra = nx.dijkstra_path(G, i, j, weight = lambda u, v, predictedLinkTraversalTime: predictedLinkTraversalTime[u,v])
#        print(f"shortestPathDijkstra is {shortestPathDijkstra} km")
        length, path = nx.single_source_dijkstra(G, i, j, weight='weight')
        print(f"  path is {path}, length is {length} km.")
        for k, u in enumerate(path[:-1]):
          v = path[k+1]
#          print(f"(k is {k})  u is {u}, v is {v}")
          if (u,v) not in arcsInShortestPaths:
            arcsInShortestPaths.append( (u,v) )
            totalAssignedSubflowsUsingArc[(u,v)] = 0.0
          totalAssignedSubflowsUsingArc[(u,v)] += subflowsToAssignedExitNodesByInjectionNodeID[nodeID1][nodeID2]
    endtime_runShortestPathsAlgorithm = time.time()
    timetorunShortestPathsAlgorithm = endtime_runShortestPathsAlgorithm - starttime_runShortestPathsAlgorithm
    print(f"timetorunShortestPathsAlgorithm is {timetorunShortestPathsAlgorithm:.5f} seconds.")
    Gmaximumflow.add_edges_from( arcsInShortestPaths )
    print(f"totalAssignedSubflowsUsingArc is {totalAssignedSubflowsUsingArc}")
    for u, v in Gmaximumflow.edges():
      Gmaximumflow.edges[u,v]['capacity'] = min( argsflowfunc['linkcapacity'][frozenset({u,v})], totalAssignedSubflowsUsingArc[(u,v)] )
#      print(f"Gmaximumflow.edges[{u},{v}]['capacity'] is {Gmaximumflow.edges[u,v]['capacity']}")

  else:  # all arcs of graph will be used in calculation of maximum flow
    Gmaximumflow.add_edges_from( arcs )
    for u, v in Gmaximumflow.edges():
#      print(f"Gmaximumflow.edges[{u},{v}]['weight'] is {Gmaximumflow.edges[u,v]['weight']}")
      Gmaximumflow.edges[u,v]['capacity'] = argsflowfunc['linkcapacity'][frozenset({u,v})]
#      print(f"Gmaximumflow.edges[{u},{v}]['capacity'] is {Gmaximumflow.edges[u,v]['capacity']}")

#  print("Gmaximumflow.nodes() is %s" % Gmaximumflow.nodes())
#  print("Gmaximumflow.edges() is %s" % Gmaximumflow.edges())
  # For each two-way link, split one arc (v1,v2) by adding a new dummy node d and replacing (v1,v2) with the pair of arcs (v1,d) and (d,v2), setting the capacity of both new arcs to that of the original arc (v1,v2):
  print("Splitting each parallel arc into a pair of new arcs with a new dummy-node in between ...")
  alreadyProcessedArcs = []
  newArcsByDummyNode = {}
#  for u in Gmaximumflow.nodes():
  for u, v in Gmaximumflow.edges():
    if (v,u) in Gmaximumflow.edges() and (v,u) not in alreadyProcessedArcs:
      dummyNode = f"dummyBetween{v}and{u}" 
      newArcsByDummyNode[ dummyNode ] = ((v,dummyNode), (dummyNode,u))
      alreadyProcessedArcs.append( (u,v) )
#  print(f"alreadyProcessedArcs is {alreadyProcessedArcs}")
#  print(f"newArcsByDummyNode is {newArcsByDummyNode}")
  for dummyNode in newArcsByDummyNode:
    ((v,dummyNode), (dummyNode,u)) = newArcsByDummyNode[ dummyNode ]
    Gmaximumflow.remove_edge(v,u)
    Gmaximumflow.add_edges_from( [(v,dummyNode), (dummyNode,u)] )
    Gmaximumflow.edges[v,dummyNode]['capacity'] = argsflowfunc['linkcapacity'][frozenset({u,v})]  # assume the capacity of each two-way link is the same in either direction
    Gmaximumflow.edges[dummyNode,u]['capacity'] = argsflowfunc['linkcapacity'][frozenset({u,v})]
#  print("After replacement of parallel arcs, Gmaximumflow.nodes() is %s" % Gmaximumflow.nodes())
#  print("Gmaximumflow.edges() is %s" % Gmaximumflow.edges())
#  for u, v in Gmaximumflow.edges():
#    print(f"Gmaximumflow.edges[{u},{v}]['capacity'] is {Gmaximumflow.edges[u,v]['capacity']}")

  totalAssignedSubflowsFromInjectionNode = {}
  totalAssignedSubflowsIntoExitNode = {}  # will store total of subflows assigned to each exit-node
  totalAssignedSubflowsFromInjectionNodes = 0.0
#  for nodeID2 in exitnodes:
#    j = pointcoordsandindexbyID[nodeID2]['index']
#    totalAssignedSubflowsIntoExitNode[j] = 0.0
  if subflowsToAssignedExitNodesByInjectionNodeID != None:
    for nodeID1 in subflowsToAssignedExitNodesByInjectionNodeID:  # only for the injection-nodes specified as keys of 'subflowsToAssignedExitNodesByInjectionNodeID' do we need to construct the shortest path to each exit-node and its length
      i = pointcoordsandindexbyID[nodeID1]['index']
      for nodeID2 in subflowsToAssignedExitNodesByInjectionNodeID[nodeID1]:
        j = pointcoordsandindexbyID[nodeID2]['index']
        subflow = subflowsToAssignedExitNodesByInjectionNodeID[nodeID1][nodeID2]
        print(f"injection-node is i={i}, exit-node is j={j}: assigned subflow is {subflow} veh/hr")
#        print(f"totalAssignedSubflowsFromInjectionNode is {totalAssignedSubflowsFromInjectionNode}")
        if i in totalAssignedSubflowsFromInjectionNode:
          totalAssignedSubflowsFromInjectionNode[i] += subflow
        else:
          totalAssignedSubflowsFromInjectionNode[i] = subflow
#        print(f"totalAssignedSubflowsFromInjectionNode[{i}] is {totalAssignedSubflowsFromInjectionNode[i]}")
        if j in totalAssignedSubflowsIntoExitNode:
          totalAssignedSubflowsIntoExitNode[j] += subflow
        else:
          totalAssignedSubflowsIntoExitNode[j] = subflow
        totalAssignedSubflowsFromInjectionNodes += subflow
      print(f"totalAssignedSubflowsFromInjectionNode[{i}] is {totalAssignedSubflowsFromInjectionNode[i]}: adding arc of that capacity from super-source (node '-1') to node {i}.")
      Gmaximumflow.add_edge(-1, i, capacity=totalAssignedSubflowsFromInjectionNode[i])  # -1 is source-node
    print(f"totalAssignedSubflowsIntoExitNode is {totalAssignedSubflowsIntoExitNode}.")
    for j in totalAssignedSubflowsIntoExitNode:
      print(f"totalAssignedSubflowsIntoExitNode[{j}] is {totalAssignedSubflowsIntoExitNode[j]}: adding arc of that capacity from node {j} to super-sink (node '-2').")
      Gmaximumflow.add_edge(j, -2, capacity=totalAssignedSubflowsIntoExitNode[j])  # -2 is sink-node
  else:  # subflowsToAssignedExitNodesByInjectionNodeID == None
    for nodeID1 in inflowsByNodeID:  # only for the injection-nodes with positive inflows do we need to  #construct the shortest path to each exit-node and its length
      i = pointcoordsandindexbyID[nodeID1]['index']
      print(f"inflowsByNodeID[{nodeID1}] is {inflowsByNodeID[nodeID1]}: adding arc of that capacity from super-source (node '-1') to node {i}.")
      Gmaximumflow.add_edge(-1, i, capacity=inflowsByNodeID[nodeID1])  # -1 is source-node
    INFINITECAPACITY = 9999999  # TODO: is this sufficiently large?
    for nodeID2 in exitnodes:
      j = pointcoordsandindexbyID[nodeID2]['index']
      print(f"Adding arc of 'infinite' capacity {INFINITECAPACITY} from node {j} to super-sink (node '-2').")
      Gmaximumflow.add_edge(j, -2, capacity=INFINITECAPACITY)  # -2 is sink-node

#  print("Gmaximumflow.nodes() is %s" % Gmaximumflow.nodes())
#  print("Gmaximumflow.edges() is %s" % Gmaximumflow.edges())
  maxPossibleFlow, maxflow_dict = nx.maximum_flow(Gmaximumflow, -1, -2)
  print(f"maxPossibleFlow is {maxPossibleFlow} veh/hr.")
#  print(f"maxflow_dict is {maxflow_dict}")
  for i in maxflow_dict[-1]:
    print(f"maxflow_dict[-1][{i}] is {maxflow_dict[-1][i]}")

  qij = {}  # store initial flows in all arcs, i.e. directed edges, that are in the shortest-path tree
  for i in [a for a in maxflow_dict if a not in (-1,-2)]:  # don't include arcs from the artificial source-node (-1) or sink-node (-2)
    for j in [b for b in maxflow_dict[i] if b not in (-1,-2)]:  # don't include arcs to the artificial source-node (-1) or sink-node (-2)
      qij[i,j] = maxflow_dict[i][j]
  for dummyNode in newArcsByDummyNode:  # remove each dummy-node, converting its new arcs back to the original parallel arc:
    ((v,dummyNode), (dummyNode,u)) = newArcsByDummyNode[ dummyNode ]
    Gmaximumflow.add_edge(v,u)
    Gmaximumflow.remove_edges_from( [(v,dummyNode), (dummyNode,u)] )
    Gmaximumflow.remove_node(dummyNode)
    Gmaximumflow.edges[v,u]['capacity'] = argsflowfunc['linkcapacity'][frozenset({u,v})]  # assume the capacity of each two-way link is the same in either direction
    assert qij[v,dummyNode] == qij[dummyNode,u]
    qij[v,u] = qij[v,dummyNode]
    del qij[v,dummyNode]
    del qij[dummyNode,u]
#  print("After restoration of every parallel arc and removal of its dummy-node, Gmaximumflow.nodes() is %s" % Gmaximumflow.nodes())
#  print("Gmaximumflow.edges() is %s" % Gmaximumflow.edges())
#  for u, v in Gmaximumflow.edges():
#    print(f"Gmaximumflow.edges[{u},{v}]['capacity'] is {Gmaximumflow.edges[u,v]['capacity']}")
#  print(f"qij is {qij}")
  outflowsByNodeID = {}
  if subflowsToAssignedExitNodesByInjectionNodeID != None:
    for j in totalAssignedSubflowsIntoExitNode:
#      outflowsByNodeID[pointNodeIDbyIndex[j]] = -totalAssignedSubflowsIntoExitNode[j]
      outflowsByNodeID[pointNodeIDbyIndex[j]] = -maxflow_dict[j][-2]
  else:
    for nodeID2 in exitnodes:
      j = pointcoordsandindexbyID[nodeID2]['index']
      outflowsByNodeID[pointNodeIDbyIndex[j]] = -maxflow_dict[j][-2]
  print(f"outflowsByNodeID is {outflowsByNodeID}")
  endtime_runMaximumFlowAlgorithm = time.time()
  timeToRunMaximumFlowAlgorithm = endtime_runMaximumFlowAlgorithm - starttime
  print(f"timeToRunMaximumFlowAlgorithm is {timeToRunMaximumFlowAlgorithm:.5f} seconds.")

  if subflowsToAssignedExitNodesByInjectionNodeID != None:
    totalSubflowsIntoAssignedExitNodes = sum(f for f in totalAssignedSubflowsIntoExitNode.values())
    assert totalAssignedSubflowsFromInjectionNodes == totalSubflowsIntoAssignedExitNodes
    if maxPossibleFlow < totalAssignedSubflowsFromInjectionNodes:
      alertCongestionInNetwork = f"ALERT: maximum possible flow through network of {maxPossibleFlow} veh/hr is exceeded by sum of subflows from injection-nodes to exit-nodes of {totalAssignedSubflowsFromInjectionNodes} veh/hr: congestion will occur."
      print(alertCongestionInNetwork)
#      alertCongestionInNetworks.append( alertCongestionInNetwork )
    else:
      assert abs(maxPossibleFlow - totalAssignedSubflowsFromInjectionNodes) <= EPS
      absenceOfCongestionMessage = f"Maximum possible flow through network of {maxPossibleFlow} veh/hr >= sum of subflows from injection-nodes to exit-nodes of {totalAssignedSubflowsFromInjectionNodes} veh/hr; N.B. congestion might still occur along subflow-paths."
      print(absenceOfCongestionMessage)
#      alertCongestionInNetworks.append( absenceOfCongestionMessage )
      alertCongestionInNetwork = None
    for i in totalAssignedSubflowsFromInjectionNode:
      if maxflow_dict[-1][i] < totalAssignedSubflowsFromInjectionNode[i]:
        alertsNodesWhereDemandNotSatisfied[i] = f"ALERT: steady-state (i.e. consistent with flows found by maximum-flow method) injection-flow {maxflow_dict[-1][i]} into injection-node {pointNodeIDbyIndex[i]} is less than total flow {totalAssignedSubflowsFromInjectionNode[i]} required there to supply its assigned exit-nodes."
  else:
    totalInflowsAtSpecifiedInjectionNodes = sum(f for f in inflowsByNodeID.values())
    if maxPossibleFlow < totalInflowsAtSpecifiedInjectionNodes:
      alertCongestionInNetwork = f"ALERT: maximum possible flow through network of {maxPossibleFlow} veh/hr is exceeded by sum of inflows at specified injection-nodes of {totalInflowsAtSpecifiedInjectionNodes} veh/hr: congestion will occur."
#      print(alertCongestionInNetwork)
#      alertCongestionInNetworks.append( alertCongestionInNetwork )
    else:
      assert abs(maxPossibleFlow - totalInflowsAtSpecifiedInjectionNodes) <= EPS
      absenceOfCongestionMessage = f"Maximum possible flow through network of {maxPossibleFlow} veh/hr >= sum of inflows at specified injection-nodes of {totalInflowsAtSpecifiedInjectionNodes} veh/hr; N.B. congestion might still occur along subflow-paths."
      print(absenceOfCongestionMessage)
#      alertCongestionInNetworks.append( absenceOfCongestionMessage )
      alertCongestionInNetwork = None
    for nodeID1 in inflowsByNodeID:
      i = pointcoordsandindexbyID[nodeID1]['index']
      if maxflow_dict[-1][i] < inflowsByNodeID[nodeID1]:
        alertsNodesWhereDemandNotSatisfied[i] = f"ALERT: steady-state (i.e. consistent with flows found by maximum-flow method) injection-flow {maxflow_dict[-1][i]} into injection-node {nodeID1} is less than specified inflow {inflowsByNodeID[nodeID1]} there."

  for (u,v) in qij:
    if abs(qij[u,v] - Gmaximumflow.edges[u,v]['capacity']) <= EPS:
#      alertsLinksPossiblyCausingCongestion[u,v] = f"ALERT: steady-state (i.e. consistent with flows found by maximum-flow method) flow {qij[u,v]} in arc ({pointNodeIDbyIndex[u]},{pointNodeIDbyIndex[v]}) is equal or close to link-capacity {Gmaximumflow.edges[u,v]['capacity']} on that link: it's possible this link's capacity being too low is a cause of congestion in the network."
      alertsLinksPossiblyCausingCongestion[u,v] = f"ALERT: steady-state (i.e. consistent with flows found by maximum-flow method) flow {qij[u,v]} in arc ({pointNodeIDbyIndex[u]},{pointNodeIDbyIndex[v]}) is equal or close to link-capacity {Gmaximumflow.edges[u,v]['capacity']} on that link."
  endtime_calculatecongestionpoints = time.time()
  timeToCalculateCongestionPoints = endtime_calculatecongestionpoints - endtime_runMaximumFlowAlgorithm
  print(f"timeToCalculateCongestionPoints is {timeToCalculateCongestionPoints:.5f} seconds.")
#  sys.exit()
###  return (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timeToRunMaximumFlowAlgorithm, timeToCalculateCongestionPoints, lK)
##  return (qij, timeToRunMaximumFlowAlgorithm, timeToCalculateCongestionPoints, lK)
#  return (qij, totalAssignedSubflowsFromInjectionNode, outflowsByNodeID, alertCongestionInNetwork, alertsNodesWhereDemandNotSatisfied, alertsLinksPossiblyCausingCongestion, timeToRunMaximumFlowAlgorithm, timeToCalculateCongestionPoints)
  return (qij, totalAssignedSubflowsFromInjectionNode, maxflow_dict, outflowsByNodeID, alertCongestionInNetwork, alertsNodesWhereDemandNotSatisfied, alertsLinksPossiblyCausingCongestion, timeToRunMaximumFlowAlgorithm, timeToCalculateCongestionPoints)


#def runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, simulDuration=14/3600):
def runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, inflowsandflowperiodsByNodeID, exitnodes, subflowsToAssignedExitNodesByInjectionNodeID, simulDuration=1.0):  # SEM4 specifies a duration (in hours) for the simulated propagation of wave- and shock-fronts
  starttime = time.time()
  print("inflowsByNodeID is", inflowsByNodeID)
  print(f"subflowsToAssignedExitNodesByInjectionNodeID is {subflowsToAssignedExitNodesByInjectionNodeID}")
  exception = None  # this will change if the global-optimisation solver is employed and it raises an exception such as numpy.linalg.LinAlgError
##  (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0, pointcoordsbyID) = setuptrafficflowproblem( JSONnetworkfilename, b0 )
#  (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, argsflowfunc, pointcoordsbyID) = setuptrafficflowproblem( JSONnetworkfilename )
#  (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, hExitNodes, argsflowfunc, pointcoordsbyID, linestringsWithCoords) = setuptrafficflowproblem( JSONnetworkfilename, exitnodes )
  (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, hExitNodes, argsflowfunc, pointcoordsandindexbyID, pointNodeIDbyIndex, linestringsWithCoords) = setuptrafficflowproblem( JSONnetworkfilename, exitnodes, SEMversion )

#  print("runSEM(): injectionnodes is", injectionnodes)
  b0 = np.zeros( len(injectionnodes) )
#  b0[3] = 1000.0  # TODO: should be determined via 'inflowsByNodeID'
#  print(f"inflowsByNodeID is {inflowsByNodeID}")
#  for i in range(N):  # node N-1 won't have any neighbours not already considered
#  for i in injectionnodes:
  for id in inflowsByNodeID:
##    b0[indexamonginjectionorexitnodes[id]] = inflowsByNodeID[id]
#    print(f"pointcoordsandindexbyID['{id}'] is {pointcoordsandindexbyID[id]}")
    b0[indexamonginjectionorexitnodes[ pointcoordsandindexbyID[id]['index'] ]] = inflowsByNodeID[ id ]
  print("b0 (inflow at injection-nodes) is", b0)
  assert len(b0) == len(injectionnodes)

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


  if SEMversion == 'SEM2':  # use Newton-Raphson or global optimiser, rather than finding shortest-path routes
    starttime = time.time()
    injectionnodeindices = [pointcoordsandindexbyID[nodeID]['index'] for nodeID in pointcoordsandindexbyID if nodeID not in exitnodes]
    print(f"injectionnodeindices is {injectionnodeindices}")
    exitnodeindices = [pointcoordsandindexbyID[nodeID]['index'] for nodeID in pointcoordsandindexbyID if nodeID in exitnodes]
    print(f"exitnodeindices is {exitnodeindices}")
#    NEWTONRAPHSON = True  # use Newton-Raphson to find h satisfying f(h) = 0
    NEWTONRAPHSON = False  # use a (local- or global-) minimisation method to find the minimum of ||f(h)||^2
    ENERGYFUNCTION = 'conservationOfMassResidual'
#    ENERGYFUNCTION = 'conservationOfMassResidualPlusTotalCongestion'
#    ENERGYFUNCTION = 'conservationOfMassResidualLessTotalCongestion'

#     The solver, whether Newton-Raphson or global minimiser, can be run more than once in order to automatically find distinct solutions h.
    NUMINITPOINTSPERDIMENSION = 1  # run the solver, either Newton-Raphson or global minimiser, exactly once

#    NUMINITPOINTSPERDIMENSION = 6  # run Newton-Raphson or global optimiser more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 11  # run Newton-Raphson or global optimiser more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 21  # run Newton-Raphson or global optimiser more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 101  # run Newton-Raphson or global optimiser more than once in order to automatically find distinct solutions h
#    NUMINITPOINTSPERDIMENSION = 1001  # run Newton-Raphson or global optimiser more than once in order to automatically find distinct solutions h
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
#      hinitList = [tuple([deltah * i for i in il]) for il in generateAllIndices( len(injectionnodes), NUMINITPOINTSPERDIMENSION ) ]
      hinitList = [tuple([deltah * i for i in il]) for il in generateAllIndices( len(injectionnodeindices), NUMINITPOINTSPERDIMENSION ) ]
      print("hinitList is %s" % hinitList)

#    else:  # NUMGLOBALOPTIMISERRUNS is 1
    else:  # NUMINITPOINTSPERDIMENSION is 1
#      hinitInjectionNodes = len(injectionnodes) * [0.0]
      hinitInjectionNodeIndices = len(injectionnodeindices) * [0.0]
      # Override the default zero-vector initial estimate:
#      hinitInjectionNodeIndices = [90]
#      hinitInjectionNodeIndices = [0, 6]
#      hinitInjectionNodeIndices = [0, 10]
#      hinitInjectionNodeIndices = [0, 12]
#      hinitInjectionNodeIndices = [0.11791782, 0.04704803]
#      hinitInjectionNodeIndices = [0.2, 0.2]
#      hinitInjectionNodeIndices = [1.6, 0.25]
#      hinitInjectionNodeIndices = [4, 2]
#      hinitInjectionNodeIndices = [6, 2]
#      hinitInjectionNodeIndices = [15.0, 20.0]
#      hinitInjectionNodeIndices = [0.0, 1.0, 2.5]
#      hinitInjectionNodeIndices = [0.0, 2.0, 2.0]
#      hinitInjectionNodeIndices = [1.0, 3.0, 3.0]

#      hinitInjectionNodes = [10*random.random(), 10*random.random()]  # choose an initial estimate randomly (when trying to find several distinct solutions)
#      hinitInjectionNodes = [3*random.random(), 3*random.random()]  # choose an initial estimate randomly (when trying to find several distinct solutions)
#      hinitInjectionNodes = [5*random.random(), 5*random.random()]  # choose an initial estimate randomly (when trying to find several distinct solutions)

#      hinitList = [ hinitInjectionNodes ]
      hinitList = [ hinitInjectionNodeIndices ]


    flowfunction = flowijfunc[flowfunc]
    derivfunction = derivijfunc[flowfunc]
    solutionsfound = []
#    NUMGLOBALOPTIMISERRUNS = NUMINITPOINTSPERDIMENSION**len(injectionnodes)
    NUMGLOBALOPTIMISERRUNS = NUMINITPOINTSPERDIMENSION**len(injectionnodeindices)
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
#      (exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, numNRiterationsused, outflowsByNodeID) = solveforhvector(NEWTONRAPHSON, N, injectionnodes, exitnodes, indexamonginjectionorexitnodes, terminatorPressure, hExitNodes, neighboursThatCanFlowIntoj, flowfunction, derivfunction, hJmin, hJmax, hinitInjectionNodes, argsflowfunc, ENERGYFUNCTION, flowijfunc, flowfunc, b0)
      (exception, hInjectionNodeIndices, qij, flowbalanced, congestionij, totalnetworkcongestion, numNRiterationsused, outflowsByNodeID) = solveforhvector(NEWTONRAPHSON, N, injectionnodeindices, exitnodeindices, indexamonginjectionorexitnodes, terminatorPressure, hExitNodes, neighboursThatCanReceiveFlowFromi, neighboursThatCanFlowIntoj, flowfunction, derivfunction, hJmin, hJmax, hinitInjectionNodes, argsflowfunc, ENERGYFUNCTION, flowijfunc, flowfunc, b0)
      endtime = time.time()
      print("congestionij is %s" % congestionij)  # these two lines check whether congestionij and totalnetworkcongestion are correctly returned by solveforhvector()
      print("Total congestion in network is %.3f" % totalnetworkcongestion)
      timetofindsoln = endtime - starttime
      sumofsolvertimes += timetofindsoln
#      print("Solver-run %d: hInjectionNodes is %s" % (w, hInjectionNodes))
      print("Solver-run %d: hInjectionNodeIndices is %s" % (w, hInjectionNodeIndices))
#      print("flowijfunc is %s" % flowijfunc)
      print("flowfunc %s, ENERGYFUNCTION %s, b0 %s" % (flowfunc, ENERGYFUNCTION, b0) )
##      print("argsflowfunc is %s" % argsflowfunc)
###      totalnetworkcongestion = totalcongestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc) )
##      (congestionij, totalnetworkcongestion) = congestioninnetwork(hInjectionNodes, injectionnodes, indexamonginjectionorexitnodes, neighboursThatCanFlowIntoj, exitnodes, hExitNodes, flowijfunc, flowfunc, b0, argsflowfunc)
##      print("congestionij is %s" % congestionij)
##      print("Total congestion in network is %.3f" % totalnetworkcongestion)
#      solutionsfound.append( (hinitInjectionNodes, exception, hInjectionNodes, qij, flowbalanced, congestionij, totalnetworkcongestion, timetofindsoln) )
      solutionsfound.append( (hinitInjectionNodes, exception, hInjectionNodeIndices, qij, flowbalanced, congestionij, totalnetworkcongestion, timetofindsoln) )
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
    nonsolutionsfound = 0  # count the number of times Newton-Raphson or global optimiser fails to find a valid, i.e. flow-balanced, solution
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


#  elif SEMversion == 'SEM3':  # don't use Newton-Raphson or global optimiser; instead find shortest-path routes
  elif SEMversion in ('SEM3', 'SEM4'):  # don't use Newton-Raphson or global optimiser; instead find shortest-path routes
#    (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK, directedNeighboursInShortestPathTreeOfNode, assignedExitNodesReachedDownstreamOfArc) = findShortestPathsSEM3or4(N, argsflowfunc, injectionnodes, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID)
    (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK, successorNodeOnShortestPathFromNodeToExitNode, assignedExitNodesReachedDownstreamOfArc) = findShortestPathsSEM3or4(N, argsflowfunc, injectionnodes, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID)

    # Every vehicle flows along a shortest path from its injection-node to an exit-node, and the traffic-flow in the entire network is the sum of the flows along these shortest paths. 

    ## (Possible future development) Could use a discretised version (see Treiber & Kesting, pp. 67-75) of the continuity equation to calculate the steady-state traffic-density and outflow at each injection-point, then at other non-leaf nodes, and then at the exit-nodes. (If a node has in-degree of at least two, then its inflows and their densities combine, and the discretised continuity equation gives the resultant outflow and its density at the far end of the adjacent down-flow link.)
    ## The link to the next node is divided into sections; at the beginning of a section the flow and speed are known, which gives the density once traffic has flown into the section, which gives, according to v = v(rho) for the link, the speed at the beginning of the next section. The flow is conserved between adjacent sections.

#    for i in injectionnodes:
#      j = shortestPathToAnyExit[i][1]
##      timeinterval = 1/3600  # in hours
##      sectionlength = v0[indexamonginjectionorexitnodes[i]] * timeinterval
#      enteringflow = b0[indexamonginjectionorexitnodes[i]]
#      enteringspeed = v0[indexamonginjectionorexitnodes[i]]
##      sectiondensity = enteringflow * timeinterval / sectionlength
#      sectiondensity = enteringflow / enteringspeed
##      print(f"After injection-flow of {enteringflow:.1f} veh/hr at node {i} has flowed towards node {j} at speed v0={v0[indexamonginjectionorexitnodes[i]]:.1f} km/hr for {timeinterval:.6f} hours, covering a distance of {sectionlength:.3f} km, the density in this section becomes {sectiondensity:.3f} veh/km.")

#    qijatnodej = {}  # store flows on all edges
    congestionij = {}  # store congestion-values on all edges
    densityijatnodej = {}  # store densities on all edges
    speedijatnodej = {}  # store densities on all edges
    for edge in lK.keys():
#      qijatnodej[edge[0], edge[1]] = 0.0
      densityijatnodej[edge] = 'linkNotUsedInAShortestPathOrLeavingEndNotReachedDuringSimul'
      speedijatnodej[edge] = 'linkNotUsedInAShortestPathOrLeavingEndNotReachedDuringSimul'
    outflowsByNodeID = {}  # store the outflows at exit-nodes, as they are calculated by the solution-method (Newton-Raphson, global optimiser, or shortest-path-based method)
    for i in exitnodes:  # define default outflow-value for each exit-node
      outflowsByNodeID[i] = f"traffic has not yet reached node {i}"

    if SEMversion == 'SEM3':
      for i in exitnodes:
#        flow, speed = flowandspeedEnteringNode( i, parentsOnShortestPathsFromInjectionNodes, qijatnodej, congestionij, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnyExit, flowijfunc, flowfunc )
        childnode = None  # an exit-node has no child-nodes on any shortest path
#        flow, speed = flowandspeedEnteringNode( i, childnode, parentsOnShortestPathsFromInjectionNodes, qijatnodej, congestionij, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, flowijfunc, flowfunc )  # get rid of injection-speed 'v0' and replace it with the link's free speed
        flow = flowEnteringNodeSEM3( i, childnode, parentsOnShortestPathsFromInjectionNodes, qijatnodej, densityijatnodej, speedijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, flowijfunc, flowfunc, injectionnodes )  # get rid of injection-speed 'v0' and replace it with the link's free speed
        outflowsByNodeID[i] = -flow
    elif SEMversion == 'SEM4':
      print(f"duration of simulation is {simulDuration} hours or {simulDuration*3600} seconds.")
#      (cumulNumVehiclesToHaveDepartedLinkDuringSimulation, outflowsByNodeID, congestionAlertMessages) = evolutionofwaveandshockfrontsSEM4( inflowsandflowperiodsByNodeID, subflowsToAssignedExitNodesByInjectionNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, densityijatnodej, speedijatnodej, flowstateijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, shortestPathToExit, flowijfunc, flowfunc, injectionnodes, outflowsByNodeID, pointcoordsandindexbyID, pointNodeIDbyIndex, nextNodeInShortestPath, directedNeighboursInShortestPathTreeOfNode, assignedExitNodesReachedDownstreamOfArc )
      (cumulNumVehiclesToHaveDepartedLinkDuringSimulation, outflowsByNodeID, congestionAlertMessages) = evolutionofwaveandshockfrontsSEM4( inflowsandflowperiodsByNodeID, subflowsToAssignedExitNodesByInjectionNodeID, simulDuration, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, densityijatnodej, speedijatnodej, flowstateijatnodej, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, shortestPathToAnyExit, shortestPathToExit, flowijfunc, flowfunc, injectionnodes, outflowsByNodeID, pointcoordsandindexbyID, pointNodeIDbyIndex, nextNodeInShortestPath, successorNodeOnShortestPathFromNodeToExitNode, assignedExitNodesReachedDownstreamOfArc )
      print(f"congestionAlertMessages is '{congestionAlertMessages}'")
      timeaveragedflowsinlinksoverSimulDuration = {key: cumulNumVehiclesToHaveDepartedLinkDuringSimulation[key] / simulDuration for key in cumulNumVehiclesToHaveDepartedLinkDuringSimulation}
      assert len(outflowsByNodeID) == len(exitnodes)
#    elif SEMversion == 'SEM5':

#    if SEMversion in ('SEM3', 'SEM4'):  # not needed
    qij = copy.deepcopy(qijatnodej)  # TODO: might need only shallow copy
    print("qij is", qij)  # check whether qij has been defined; might need to set it to qijatnodej

    print("inflowsandflowperiodsByNodeID is", inflowsandflowperiodsByNodeID)
    print("argsflowfunc['linkcapacity'] is", argsflowfunc['linkcapacity'])
    print("qijatnodei is", qijatnodei)
    print("qijatnodej is", qijatnodej)
    print("outflowsByNodeID is", outflowsByNodeID)
    if SEMversion == 'SEM3':
#      print(f"flow is {flow:.3f}, speed is {speed:.3f}.")
      print(f"flow is {flow:.3f}.")
    elif SEMversion == 'SEM4':
      print(f"cumulNumVehiclesToHaveDepartedLinkDuringSimulation is {cumulNumVehiclesToHaveDepartedLinkDuringSimulation}.")
      print(f"timeaveragedflowsinlinksoverSimulDuration is {timeaveragedflowsinlinksoverSimulDuration}.")
    endtime = time.time()
    timetorunshortestpathbasedmethod = endtime - starttime
    print(f"Shortest-path-based method took {timetorunshortestpathbasedmethod:.5f} seconds, of which {timetorunShortestPathsAlgorithm:.5f} seconds was spent running Floyd-Warshall algorithm and {timetocalculateshortestpaths:.5f} seconds was spent calculating the shortest paths.")
    print("densityijatnodej is", densityijatnodej)
#    print("speedijatnodej is", speedijatnodej)
    starttime = time.time()
#    print("Calculating congestion in edges...")
##    print(f"qijatnodej.keys() is {qijatnodej.keys()}")
#    print(f"densityijatnodej.keys() is {densityijatnodej.keys()}")
    for edge in lK.keys():
##      humpdensity = argsflowfunc['tau'][frozenset({p,i})]
##      if densityinlink > humpdensity:
##        congestiononlink = linkcapacity - abs(flowleavingparent)
      (r, s) = edge
      if (r, s) in qijatnodej:
        flowinedge = qijatnodej[r, s]
#      if flowinedge < 0.0:
      elif (s, r) in qijatnodej:
        (r, s) = (s, r)
        flowinedge = qijatnodej[r, s]
      else:
        congestionij[frozenset({r, s})] = 'linkNotUsedInAShortestPath'
        print(f"Congestion on link ({r}, {s}) is 'linkNotUsedInAShortestPath'.")
        continue
#      print(f"flowinedge qijatnodej[{r}, {s}] is {flowinedge}")
#      print(f"densityijatnodej[{edge}] is {densityijatnodej[edge]}")
      if densityijatnodej[edge] == 'linkNotUsedInAShortestPathOrLeavingEndNotReachedDuringSimul':  # qijatnodej has a value for the edge, so it must be used in a shortest path, which implies the arc's leaving-end (finish-node) wasn't reached during the duration of the simulation
        assert flowinedge == 0.0
        congestiononlink = 'leavingEndOfLinkNotReachedDuringSimul'
        densityijatnodej[edge] = 'leavingEndOfLinkNotReachedDuringSimul'
        speedijatnodej[edge] = 'leavingEndOfLinkNotReachedDuringSimul'
      else:
        congestiononlink = linkcongestion(r, s, flowinedge, densityijatnodej[edge], argsflowfunc)
      congestionij[frozenset({r, s})] = congestiononlink
      assert congestiononlink == 'leavingEndOfLinkNotReachedDuringSimul' or congestiononlink >= 0.0
      linkcapacity = argsflowfunc['linkcapacity'][edge]
#      if congestiononlink != 'leavingEndOfLinkNotReachedDuringSimul' and abs(congestiononlink) > 0.0:  # TODO: is the abs() necessary?
#        print(f"CONGESTION of {linkcapacity:.3f}-{flowinedge:.3f}={congestiononlink:.3f} occurs on link ({r}, {s}) which has capacity of {linkcapacity:.3f}.")  # TODO: this measure of congestion gives a negative value when the flow entering the parent is greater than the link-capacity: might need a better measure
#      elif congestiononlink == 'leavingEndOfLinkNotReachedDuringSimul':
#        print(f"congestion on link ({r}, {s}), which has capacity of {linkcapacity:.3f}, is {congestiononlink}.")
#      else:
#        print(f"There is NO congestion on link ({r}, {s}) which has capacity of {linkcapacity:.3f}.")
    endtime = time.time()
    timetocalculatecongestion = endtime - starttime
    print(f"Calculation of congestion took {timetocalculatecongestion:.5f} seconds.")
#    if SEMversion == 'SEM4':
##      origkeys_qij = qijatnodej.keys() 
##      for (a,b) in origkeys_qij:
#      for (a,b) in qijatnodej:
#        qijatnodej[a,b] = timeaveragedflowsinlinksoverSimulDuration[a,b]
#      print(f"qijatnodej is now defined to be the time-averaged flows over simulDuration of {simulDuration} hours ({simulDuration*3600:.1f} seconds), qijatnodej={qijatnodej}")


  elif SEMversion == 'SEM5':  # find the paths that maximise flow through the network
###    (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK) = findShortestPathsSEM3or4(N, argsflowfunc, injectionnodes, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID)  # TODO: this will be too slow: replace with Dijkstra's algorithm
##    (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK) = findShortestPathsSEM5(N, argsflowfunc, injectionnodes, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID)
#    (nextNodeInShortestPath, shortestDistanceToAnyExit, shortestPathToExit, shortestPathToAnyExit, parentsOnShortestPathsFromInjectionNodes, qijatnodei, densityijatnodei, qijatnodej, flowstateijatnodej, timetorunShortestPathsAlgorithm, timetocalculateshortestpaths, lK) = findShortestPathsSEM5(argsflowfunc, exitnodes, pointcoordsandindexbyID, subflowsToAssignedExitNodesByInjectionNodeID)

###            flowstateijatnodej  # is this important to have?
    (qij, totalAssignedSubflowsFromInjectionNode, maxflow_dict, outflowsByNodeID, alertCongestionInNetwork, alertsNodesWhereDemandNotSatisfied, alertsLinksPossiblyCausingCongestion, timeToRunMaximumFlowAlgorithm, timeToCalculateShortestPaths) = findMaximumFlowSEM5(argsflowfunc, exitnodes, pointcoordsandindexbyID, pointNodeIDbyIndex, inflowsByNodeID, subflowsToAssignedExitNodesByInjectionNodeID)
    print(f"alertCongestionInNetwork is '{alertCongestionInNetwork}'")
    print(f"alertsNodesWhereDemandNotSatisfied is {alertsNodesWhereDemandNotSatisfied}")
    print(f"alertsLinksPossiblyCausingCongestion is {alertsLinksPossiblyCausingCongestion}")
    print(f"There are {len(alertsLinksPossiblyCausingCongestion)} links with flow at capacity.")



  # Now that the congestion-finding method has run, save the output in GeoJSON format:
  outputGeoJSON = {}  # fill this with GeoJSON data:
#  BBwest = min([pointcoordsbyID[i][0] for i in range(N)])  # bounding-box limits
#  BBsouth = min([pointcoordsbyID[i][1] for i in range(N)])
#  BBeast = max([pointcoordsbyID[i][0] for i in range(N)])
#  BBnorth = max([pointcoordsbyID[i][1] for i in range(N)])
#pointcoordsandindexbyID[pointDict['matsimID']] = {'coords': xy, 'index': pointDict['index']}
  BBwest = min([pointcoordsandindexbyID[nodeID]['coords'][0] for nodeID in pointcoordsandindexbyID])  # bounding-box limits
  BBsouth = min([pointcoordsandindexbyID[nodeID]['coords'][1] for nodeID in pointcoordsandindexbyID])
  BBeast = max([pointcoordsandindexbyID[nodeID]['coords'][0] for nodeID in pointcoordsandindexbyID])
  BBnorth = max([pointcoordsandindexbyID[nodeID]['coords'][1] for nodeID in pointcoordsandindexbyID])
  print(f"BBwest {BBwest}, BBsouth {BBsouth}, BBeast {BBeast}, BBnorth {BBnorth}")
  outputGeoJSON["bbox"] = [BBwest, BBsouth, BBeast, BBnorth]
  outputGeoJSON["features"] = []
#  print("injectionnodes is", injectionnodes)
  assert len(inflowsByNodeID) <= len(injectionnodes)
  for i in [j for j in injectionnodes if j not in inflowsByNodeID]:
    inflowsByNodeID[i] = 0.0
#  print("inflowsByNodeID is", inflowsByNodeID)
  assert len(inflowsByNodeID) == len(injectionnodes)
  nonzeroinflowsByNodeID = {nodeID: inflow for (nodeID, inflow) in inflowsByNodeID.items() if inflow > 0.0}
  print("nonzeroinflowsByNodeID is", nonzeroinflowsByNodeID)
#  print("sorted( pointcoordsandindexbyID.keys() ) is", sorted( pointcoordsandindexbyID.keys() ) )
#  for nodeID in sorted( pointcoordsbyID.keys() ):
  for nodeID in sorted( pointcoordsandindexbyID.keys() ):
    pointJSON = {}
    pointJSON["geometry"] = {}
#    pointJSON["geometry"]["coordinates"] = list(pointcoordsbyID[nodeID])
    pointJSON["geometry"]["coordinates"] = list(pointcoordsandindexbyID[nodeID]['coords'])
    pointJSON["geometry"]["type"] = "Point"
    pointJSON["properties"] = {}
#    print(f"injectionnodes[:6] is {injectionnodes[:6]}")
    if nodeID in injectionnodes:
    # Previous version of output-JSON used "flow", but it's too easily confused with the "flow" property of links so I changed it to "inflow":
      if SEMversion in ('SEM1', 'SEM2', 'SEM3', 'SEM4'):
        pointJSON["properties"]["inflow"] = inflowsByNodeID[nodeID]
      elif SEMversion == 'SEM5':
        i = pointcoordsandindexbyID[nodeID]['index']
#        if i in totalAssignedSubflowsFromInjectionNode:
        if nodeID in nonzeroinflowsByNodeID:
#          print(f"maxflow_dict[-1][{i}] is {maxflow_dict[-1][i]}")
          pointJSON["properties"]["inflow"] = maxflow_dict[-1][i]
        else:
          pointJSON["properties"]["inflow"] = 0.0
      if SEMversion == 'SEM2':  # Newton-Raphson or the global optimiser was used
#        pointJSON["properties"]["head"] = hInjectionNodes[indexamonginjectionorexitnodes[nodeID]]
        pointJSON["properties"]["head"] = hInjectionNodes[indexamonginjectionorexitnodes[pointcoordsandindexbyID[nodeID]['index']]]
      if SEMversion == 'SEM5':
        i = pointcoordsandindexbyID[nodeID]['index']
        if i in alertsNodesWhereDemandNotSatisfied:
          pointJSON["properties"]["alert"] = alertsNodesWhereDemandNotSatisfied[i]
    else:  # the current node is an exit-node
      if SEMversion in ('SEM1', 'SEM2'):  # Newton-Raphson or the global optimiser was used
        pointJSON["properties"]["inflow"] = outflowsByNodeID[pointcoordsandindexbyID[nodeID]['index']]
#      elif SEMversion in ('SEM3', 'SEM4', 'SEM5'):
      else:
        pointJSON["properties"]["inflow"] = outflowsByNodeID[nodeID]
      if SEMversion == 'SEM2':  # Newton-Raphson or the global optimiser was used
#        pointJSON["properties"]["head"] = hExitNodes[indexamonginjectionorexitnodes[nodeID]]
        pointJSON["properties"]["head"] = hExitNodes[indexamonginjectionorexitnodes[pointcoordsandindexbyID[nodeID]['index']]]
#    pointJSON["properties"]["id"] = nodeID
    pointJSON["properties"]["matsim_nodeID"] = nodeID
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
###    if nodeID in injectionnodes:  # the current node is an injection-node
###    linkJSON["properties"]["K"] = linkproperties['freespeed']  # might be incorrect - is K a property of the water-flow-based monotonic-flow model?
##    print(f"linkproperties['diameter'] is {linkproperties['diameter']}")
#    linkJSON["properties"]["diameter"] = linkproperties['diameter']
##    print(f'linkJSON["properties"]["diameter"] is {linkJSON["properties"]["diameter"]}')
    (a, b) = tuple(linkproperties['endpointNodeIDs'])
    if (a,b) in qij:
      if SEMversion == 'SEM2':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = qij[a, b]
      elif SEMversion == 'SEM3':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = densityijatnodej[frozenset({a, b})] * speedijatnodej[frozenset({a, b})]
      elif SEMversion == 'SEM4':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = timeaveragedflowsinlinksoverSimulDuration[a, b]
      elif SEMversion == 'SEM5':  # found maximum possible flow in network
        linkJSON["properties"]["flow"] = qij[a, b]
        if (a,b) in alertsLinksPossiblyCausingCongestion:
          linkJSON["properties"]["alert"] = alertsLinksPossiblyCausingCongestion[a,b]
    elif (b,a) in qij:
      if SEMversion == 'SEM2':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -qij[b, a]
      elif SEMversion == 'SEM3':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -densityijatnodej[frozenset({b, a})] * speedijatnodej[frozenset({b, a})]
      elif SEMversion == 'SEM4':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -timeaveragedflowsinlinksoverSimulDuration[b, a]
      elif SEMversion == 'SEM5':  # found maximum possible flow in network
        linkJSON["properties"]["flow"] = -qij[b, a]
        if (b,a) in alertsLinksPossiblyCausingCongestion:
          linkJSON["properties"]["alert"] = alertsLinksPossiblyCausingCongestion[b,a]
    else:
#      try:
#        linkJSON["properties"]["flow"] = -densityijatnodej[frozenset({b, a})] * speedijatnodej[frozenset({b, a})]
#      except:
#        print(f"densityijatnodej[frozenset({{{b}, {a}}})] is {densityijatnodej[frozenset({b, a})]}, speedijatnodej[frozenset({{{b}, {a}}})] is {speedijatnodej[frozenset({b, a})]}")
#        raise
      linkJSON["properties"]["flow"] = 'linkNotUsedInAShortestPath'  # TODO: should the 'flow' in this case be 0.0?
#    linkJSON["properties"]["k"] =  # TODO: necessary for SEM2?
    linkJSON["properties"]["length"] = linkproperties['length']
    linkJSON["properties"]["matsim_linkID"] = linkproperties['matsim_linkID']
    linkJSON["properties"]["type"] = 1  # TODO: is the "type" property always 0 for a LineString?
    linkJSON["type"] = "Feature"
    outputGeoJSON["features"].append(linkJSON)

  outputGeoJSON["type"] = "FeatureCollection"

#  print("outputGeoJSON is", outputGeoJSON)
#  if SEMversion == 'SEM5':
#    print(f"alertCongestionInNetwork is '{alertCongestionInNetwork}'")


#  DRAWSOLUTIONGRAPH = True
  DRAWSOLUTIONGRAPH = False
  if DRAWSOLUTIONGRAPH and exception == None:
#     Draw the network as a directed graph:
    starttime = time.time()
#    networktype = 'linear'
#    networktype = 'highway'
#    networktype = 'star'
    networktype = 'nontree'
    print("networktype is", networktype)
    DRAWNODECOLOURS = True
    DRAWINTERMEDIATENODELABELS = False
    DRAWNODELABELS = False
    DRAWARCSANDARROWS = True
#    DRAWARCLABELS = True  # show flow, density, speed on edges/arcs
    DRAWARCLABELS = False  # do not show flow, density, speed on edges/arcs
#    print("pointcoordsbyID is", pointcoordsbyID)
#    if exception == None:
    debug("qij is %s" % qij)
#    print("congestionij is %s" % congestionij)
    if SEMversion in ('SEM3', 'SEM4'):  # found shortest-path routes, rather than using the global optimiser
#      print("densityijatnodej is %s" % densityijatnodej)
#      print("speedijatnodej is %s" % speedijatnodej)
      pass
    Gsolution = nx.DiGraph()  #Gsolution = nx.Graph()
#    DECIMALPLACES = 4  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 3  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 2  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 1  # number of decimal places to show in values on the graph
    DECIMALPLACES = 0  # number of decimal places to show in values on the graph

##    Gsolution.add_nodes_from(range(4))
###    Gsolution.add_edges_from( [e for e in Gsolution.edges()] )
###    Gsolution.add_edges_from( qij.keys() )
##    positiveflowarcs = [(i,j) if qij[i,j] >= 0.0 else (j,i) for (i,j) in qij]
    arcsconsideredbySEM = []
    positiveflowarcs = []
    elabels = {}
    for (i,j) in qij:
#      print("(i,j) (%d, %d)" % (i, j))
#      print("qij %f" % qij[i,j])
      if SEMversion in 'SEM4':  # time-dependent flows
        flowrounded = round(timeaveragedflowsinlinksoverSimulDuration[i,j], DECIMALPLACES)
#      elif SEMversion in ('SEM1', 'SEM2', 'SEM3', 'SEM5'):
      else:
        flowrounded = round(qij[i,j], DECIMALPLACES)
#      print(f"flowrounded is {flowrounded} veh/hr, congestionrounded is {congestionrounded} veh/hr")
      if SEMversion in ('SEM2', 'SEM3', 'SEM4'):
#        print(f"congestionij[frozenset(set(({i},{j})))] is {congestionij[frozenset(set((i,j)))]}")
        if congestionij[frozenset(set((i,j)))] == 'leavingEndOfLinkNotReachedDuringSimul':
          congestionrounded = "ENR"
        else:
          congestionrounded = round(congestionij[frozenset(set((i,j)))], DECIMALPLACES)
      if SEMversion == 'SEM2':  # used Newton-Raphson or global optimiser, rather than finding shortest-path routes
##        hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
##        hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
#        if i in exitnodes:
        if i in exitnodeindices:
          hi = terminatorPressure
        else:
          hi = hInjectionNodes[indexamonginjectionorexitnodes[i]]
#        if j in exitnodes:
        if j in exitnodeindices:
          hj = terminatorPressure
        else:
          hj = hInjectionNodes[indexamonginjectionorexitnodes[j]]
##        print("hi %f, hj %f" % (hi, hj))
        densityonlink = abs(hi - hj)
        if densityonlink == 'leavingEndOfLinkNotReachedDuringSimul':
          densityrounded = "ENR"
        else:
          densityrounded = round(densityonlink, DECIMALPLACES)
      elif SEMversion in ('SEM3', 'SEM4'):
        densityonlink = densityijatnodej[frozenset({i, j})]
#        print("densityonlink is %s" % densityonlink)
        if densityonlink == 'leavingEndOfLinkNotReachedDuringSimul':
          densityrounded = "ENR"
        else:
          densityrounded = round(densityonlink, DECIMALPLACES)
#        print(f"densityrounded is {densityrounded}")
      if SEMversion == 'SEM2':  # used Newton-Raphson or global optimiser, rather than finding shortest-path routes
        speedonlink = qij[i,j] / densityonlink
#      else:
      elif SEMversion in ('SEM3', 'SEM4'):
        speedonlink = speedijatnodej[frozenset({i, j})]
#      print("speedonlink is %s" % speedonlink)
      if SEMversion in ('SEM2', 'SEM3', 'SEM4'):
        if speedonlink == 'leavingEndOfLinkNotReachedDuringSimul':
          speedrounded = "ENR"
        else:
          speedrounded = round(speedonlink, DECIMALPLACES)
#        print("speedrounded is %s" % speedrounded)
      if qij[i,j] >= 0.0:
        arcsconsideredbySEM.append( (i,j) )
        flowvaluetodisplay = flowrounded
#        elabels[i,j] = "q %.*f, \u03b2 %.*f, v %.*f" % (DECIMALPLACES, flowrounded, DECIMALPLACES, congestionrounded, DECIMALPLACES, -speedrounded)  # display q, \beta_{ij} (congestion), and speed on each edge {i,j}
      else:
#      elif SEMversion in ('SEM1', 'SEM2'):
        arcsconsideredbySEM.append( (j,i) )
        flowvaluetodisplay = -flowrounded
#        elabels[i,j] = "q %.*f, \u03b2 %.*f, v %.*f" % (DECIMALPLACES, -flowrounded, DECIMALPLACES, congestionrounded, DECIMALPLACES, -speedrounded)  # display q, \beta_{ij} (congestion), and speed on each edge {i,j}
      if qij[i,j] > 0.0:
        positiveflowarcs.append( (i,j) )
      if SEMversion in ('SEM2', 'SEM3', 'SEM4'):
        if densityrounded == 'ENR':
          assert densityrounded == speedrounded == congestionrounded == 'ENR'
##          elabels[i,j] = "q %.*f, \u03c1 %s, v %s, \u03b2 %s" % (DECIMALPLACES, flowvaluetodisplay, densityrounded, speedrounded, congestionrounded)  # display q, rho (density), speed, and beta (congestion) on each edge
#          elabels[i,j] = "q %.*f, \u03c1=v=\u03b2=%s" % (DECIMALPLACES, flowvaluetodisplay, densityrounded)  # display q, and rho (density) on each edge
          elabels[i,j] = "q %.*f, \u03c1=v=%s" % (DECIMALPLACES, flowvaluetodisplay, densityrounded)  # display q, and rho (density) on each edge
        else:
          if abs(densityrounded) < 1e4:
#            elabels[i,j] = "q %.*f, \u03c1 %.*f, v %.*f, \u03b2 %.*f" % (DECIMALPLACES, flowvaluetodisplay, DECIMALPLACES, densityrounded, DECIMALPLACES, speedrounded, DECIMALPLACES, congestionrounded)  # display q, rho (density), speed, and beta (congestion) on each edge
            elabels[i,j] = "q %.*f, \u03c1 %.*f, v %.*f" % (DECIMALPLACES, flowvaluetodisplay, DECIMALPLACES, densityrounded, DECIMALPLACES, speedrounded)  # display q, rho (density), and speed on each edge
          else:
#            elabels[i,j] = "q %.*f, \u03c1 %.1e, v %.*f, \u03b2 %.*f" % (DECIMALPLACES, flowvaluetodisplay, densityrounded, DECIMALPLACES, speedrounded, DECIMALPLACES, congestionrounded)  # display q, rho (density), speed, and beta (congestion) on each edge
            elabels[i,j] = "q %.*f, \u03c1 %.1e, v %.*f" % (DECIMALPLACES, flowvaluetodisplay, densityrounded, DECIMALPLACES, speedrounded)  # display q, rho (density), and speed on each edge
      elif SEMversion == 'SEM5':
        linkcapacity = argsflowfunc['linkcapacity'][frozenset({i,j})]
#        elabels[i,j] = "q %.*f [%.*f]" % (DECIMALPLACES, flowvaluetodisplay, DECIMALPLACES, linkcapacity)  # display q and capacity on each arc
        if flowvaluetodisplay > 0:
          elabels[i,j] = "%.*f  [%.*f]" % (DECIMALPLACES, flowvaluetodisplay, DECIMALPLACES, linkcapacity)  # display q and capacity on each positive-flow arc
#          elabels[i,j] = "%.*f" % (DECIMALPLACES, flowvaluetodisplay)  # display q on each arc
        else:
#          elabels[i,j] = "[%.*f]" % (DECIMALPLACES, linkcapacity)  # display only the capacity on each zero-flow arc
          elabels[i,j] = ""
    Gsolution.add_edges_from( arcsconsideredbySEM )
    debug("Gsolution.nodes() is %s" % Gsolution.nodes())
    for additionalNode in [node for node in range(N) if node not in Gsolution.nodes()]:
      Gsolution.add_node( additionalNode )  # in case not all nodes were originally included in Gsolution.nodes() when 'arcsconsideredbySEM' was added to the edges of Gsolution, as can happen if, for example, an injection-node has no path to an exit-node
    debug("Gsolution.nodes() is now %s" % Gsolution.nodes())
#    allnodes = range(N)
#    print("allnodes is %s" % list(allnodes))
    debug("Gsolution.edges() is %s" % Gsolution.edges())
    debug("arcsconsideredbySEM is %s" % arcsconsideredbySEM)
    debug("positiveflowarcs is %s" % positiveflowarcs)
###    print("elabels is %s" % elabels)
####    colourvalmap = {0: 1.0, 2: 0.5714285714285714, 3: 0.0}
###    colourvalmap = {}
###    for t in exitnodes:
###      colourvalmap[t] = 0.55  # colour terminal nodes red
###    nodecolourvalues = [colourvalmap.get(node, 0.45) for node in Gsolution.nodes()]  # colour non-terminal nodes light blue
##    nodecolourvalues = ['lightblue' if node in exitnodes else 'yellow' for node in Gsolution.nodes()]  # colour non-terminal nodes light blue
#    nodecolourvalues = ['lightblue' if pointNodeIDbyIndex[node] in exitnodes else 'yellow' for node in Gsolution.nodes()]  # colour exit-nodes light blue
    nodecolourvalues = []  # colour exit-nodes light blue, and injection-nodes with positive inflows yellow; make other nodes transparent
    for node in Gsolution.nodes():
      if pointNodeIDbyIndex[node] in exitnodes:
        nodecolourvalues.append('lightblue')
##      elif pointNodeIDbyIndex[node] in subflowsToAssignedExitNodesByInjectionNodeID:
#      elif pointNodeIDbyIndex[node] in inflowsByNodeID:
      elif pointNodeIDbyIndex[node] in nonzeroinflowsByNodeID:
#        nodecolourvalues.append('yellow')
        nodecolourvalues.append('red')
      else:
#        if DRAWINTERMEDIATENODELABELS:
#          nodecolourvalues.append('yellow')
#        else:
#          nodecolourvalues.append('#00000000')
        nodecolourvalues.append('#00000000')
#    nodecolourvalues = ['lightblue' if pointNodeIDbyIndex[node] in exitnodes else 'yellow' for node in allnodes]  # colour exit-nodes light blue
    debug("nodecolourvalues is %s" % nodecolourvalues)

#    red_edges = [(0, 1), (1, 2)]
#    black_edges = [edge for edge in Gsolution.edges() if edge not in red_edges]
    inflowlabels = {}
    hvaluelabels = {}
#    print("exitnodes is", exitnodes)
#    print("injectionnodes is", injectionnodes)
#    print("b0 is", b0)

    for i in range(N):  # node N-1 won't have any neighbours not already considered
#      if i in exitnodes:
      if pointNodeIDbyIndex[i] in exitnodes:
#        inflowlabels[i] = terminatorPressure  # not used at exit-nodes
#        hvaluelabels[i] = "h = %.*f" % (DECIMALPLACES, terminatorPressure)
        hvaluelabels[i] = "h = %.0f" % terminatorPressure
      else:
        if b0[indexamonginjectionorexitnodes[i]] > 0.0:  # don't display an inflow if it's zero
          inflowlabels[i] = b0[indexamonginjectionorexitnodes[i]]
          print(f"inflowlabels[{i}] is now {inflowlabels[i]}")
        if SEMversion == 'SEM2':  # used Newton-Raphson or global optimiser, rather than finding shortest-path routes
#          hvaluerounded = round(hInjectionNodes[indexamonginjectionorexitnodes[i]], DECIMALPLACES)
          hvaluelabels[i] = "h = %.*f" % (DECIMALPLACES, hInjectionNodes[indexamonginjectionorexitnodes[i]])
#          hvaluelabels[i] = "h = %e" % hInjectionNodes[indexamonginjectionorexitnodes[i]]

#    print("inflowlabels is %s" % inflowlabels)
    if SEMversion == 'SEM2':  # used Newton-Raphson or global optimiser, rather than finding shortest-path routes
      debug("hvaluelabels is %s" % hvaluelabels)

#     Need to create a layout when making separate calls to draw nodes and edges.
    if networktype == 'linear':
      pos = nx.bipartite_layout(Gsolution, Gsolution.nodes(), align='horizontal')
#      pos = nx.spectral_layout(Gsolution)
#      pos = nx.planar_layout(Gsolution)
#      pos = nx.kamada_kawai_layout(Gsolution)
#      pos = nx.shell_layout(Gsolution)
#      pos = nx.spiral_layout(Gsolution)
#      pos = nx.circular_layout(Gsolution)
    elif networktype == 'highway':  # for plotting as bipartite set, put highway nodes in first subset, remaining nodes in second subset
##      print("Gsolution.nodes()[:-1] is %s" % repr(Gsolution.nodes()[:-1]) )  # incorrect: type(Gsolution.nodes) is <class 'networkx.classes.reportviews.NodeView'>
#      print("type(Gsolution.nodes.items()) is %s" % type(Gsolution.nodes.items()))
#      print("type(Gsolution.nodes()) is %s" % type(Gsolution.nodes()))
#      print("list(Gsolution.nodes(data=True)) is %s" % list(Gsolution.nodes(data=True)))
      nodesinfirstset = [n for n in Gsolution.nodes() if n < len(Gsolution.nodes())-1]
      debug("nodesinfirstset is %s" % nodesinfirstset)
#      pos = nx.bipartite_layout(Gsolution, Gsolution.nodes()[:-1], align='horizontal')
      pos = nx.bipartite_layout(Gsolution, nodesinfirstset, align='horizontal')
      pos[4][0] += .75  # ad-hoc line for five-node ``highway'' network
    elif networktype == 'star':
      pos = nx.spring_layout(Gsolution)
    elif networktype == 'nontree':
##      pos = nx.spring_layout(Gsolution)
#      pos = nx.spring_layout(Gsolution, seed=0)  # set 'seed' for deterministic node-layout
      pos = nx.spring_layout(Gsolution, seed=0, iterations=1)  # set 'seed' for deterministic node-layout
      print("type(pos) is {type(pos)}")
    try:
#      assert len(pos) == len(pointcoordsbyID)
      assert len(pos) == len(pointcoordsandindexbyID)
    except:
      print(f"WARNING: len(pos) is {len(pos)}, len(pointcoordsandindexbyID) is {len(pointcoordsandindexbyID)}, which are not equal: a potential problem or not?")
#      raise

##    assert set(pos.keys()) == set(pointcoordsbyID.keys())
#    print(f"set(pos.keys()) is {set(pos.keys())}")
#    print(f"pointcoordsandindexbyID.keys() is {pointcoordsandindexbyID.keys()}")
##    assert set(pos.keys()) == set(pointcoordsandindexbyID.keys())
#    assert len(pos.keys()) == len(pointcoordsandindexbyID.keys())

#    for nodeID in pointcoordsbyID:  # set the nodes' positions
#      pos[nodeID] = np.array( pointcoordsbyID[nodeID] )
    for nodeID in pointcoordsandindexbyID:  # set the nodes' positions
      pos[ pointcoordsandindexbyID[nodeID]['index'] ] = np.array( pointcoordsandindexbyID[nodeID]['coords'] )
###    plt.figure(1, figsize=(9,5))
##    plt.figure(figsize=(9,5))
#    fig, ax = plt.subplots(figsize=(9,5))  # for report
    fig, ax = plt.subplots(figsize=(11,5))  # for presentation
    xmin, xmax, ymin, ymax = plt.axis()
#    print("xmin %.2f, xmax %.2f, ymin %.2f, ymax  %.2f" % (xmin, xmax, ymin, ymax))
    if DRAWNODECOLOURS:
#      nx.draw_networkx_nodes(Gsolution, pos, cmap=plt.get_cmap('jet'), node_color = nodecolourvalues, node_size = 500)
      nx.draw_networkx_nodes(Gsolution, pos, node_color = nodecolourvalues, node_size = 500)
##    nodelabels = {v: 'v'+str(v) for v in Gsolution.nodes()}
#    nodelabels = {v: str(v) if pointNodeIDbyIndex[v] in exitnodes else '' for v in Gsolution.nodes()}
    nodelabels = {}
    for v in Gsolution.nodes():
      if pointNodeIDbyIndex[v] in exitnodes:
        nodelabels[v] = str(v)
#      elif pointNodeIDbyIndex[v] in subflowsToAssignedExitNodesByInjectionNodeID:
      elif pointNodeIDbyIndex[v] in inflowsByNodeID:
        nodelabels[v] = str(v)
      else:
#        nodelabels[v] = ''
        if DRAWINTERMEDIATENODELABELS:
          nodelabels[v] = str(v)
        else:
          nodelabels[v] = ''
    debug("nodelabels is %s" % nodelabels)
    if DRAWNODELABELS:
      nx.draw_networkx_labels(Gsolution, pos, nodelabels)

#    if networktype == 'linear':
#      pos_hvalues = {v: 1*(coords+np.array([0.0, -.005])) for v, coords in pos.items()}
#    elif networktype == 'highway':
#      pos_hvalues = {v: 1*(coords+np.array([0.0, -.105])) for v, coords in pos.items()}
#    elif networktype == 'star':
#      pos_hvalues = {v: coords+np.array([0.155, -.105]) for v, coords in pos.items()}
#    elif networktype == 'nontree':
#      pos_hvalues = {v: coords+np.array([0.002, -.005]) for v, coords in pos.items()}
#    print("pos_hvalues is %s" % pos_hvalues)
    BBwest = min([pos[i][0] for i in range(N)])  # bounding-box limits
    BBeast = max([pos[i][0] for i in range(N)])
    BBnorth = max([pos[i][1] for i in range(N)])
    BBsouth = min([pos[i][1] for i in range(N)])
#    print(f"BBwest {BBwest}, BBeast {BBeast}, BBnorth {BBnorth}, BBsouth {BBsouth}")
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
##      print("np.array([x, y] is", np.array([x, y]))
#    print("pos_hvalues is %s" % pos_hvalues)
    if SEMversion == 'SEM2':  # used Newton-Raphson or global optimiser, rather than finding shortest-path routes
      nx.draw_networkx_labels(Gsolution, pos_hvalues, hvaluelabels, font_color='k')

    if DRAWARCLABELS:
#      nx.draw_networkx_edge_labels(Gsolution, pos, edge_labels=elabels)
      nx.draw_networkx_edge_labels(Gsolution, pos, edge_labels=elabels, font_color='g')
    if DRAWARCSANDARROWS:
##      nx.draw_networkx_edges(Gsolution, pos, edgelist=red_edges, edge_color='r', arrows=True)
#      nx.draw_networkx_edges(Gsolution, pos, edgelist=arcsconsideredbySEM, edge_color='b', arrows=True, arrowsize=25)  # solid (non-transparent) blue
      nx.draw_networkx_edges(Gsolution, pos, edgelist=arcsconsideredbySEM, edge_color='#0000FF22', arrows=True, arrowsize=25)  # transparent blue
#      nx.draw_networkx_edges(Gsolution, pos, edgelist=black_edges, arrows=False)
      nx.draw_networkx_edges(Gsolution, pos, edgelist=positiveflowarcs, edge_color='b', arrows=True, arrowsize=25)  # solid (non-transparent) blue

    for i in range(N):
      if i not in exitnodes:
#        if networktype == 'linear':
#          x = (pos[i][0]+1.58)/3.12  # for ``linear'' network
#          y = pos[i][1]+.57
#        elif networktype == 'highway':
#          x = (pos[i][0]+.86)/2.03  # for five-node ``highway'' network
#          y = pos[i][1]+.57
#        elif networktype == 'star':
#          x = (pos[i][0]+.86)/2.03
#          y = pos[i][1]+.57
#        elif networktype == 'nontree':
#          x = (pos[i][0]+1.06)/2.03 * .8 + .10
#          y = (pos[i][1]+.85) / 1.85 * .8 + .12
        if BBeast-BBwest > 0:
          x = (pos[i][0]-BBwest)/(BBeast-BBwest)*(.74-.2)+.2
        else:
          x = .5
        if BBnorth-BBsouth > 0:
          y = (pos[i][1]-BBsouth)/(BBnorth-BBsouth)*(.8-.2)+.24
        else:
#          y = .5
          y = .55
#        print("%d: x %.5f, y %.5f, inflowlabels[i] %.4f" % (i, x, y, inflowlabels[i]))
        if i in inflowlabels:  # don't display an inflow if it's zero
          fig.text(x, y, str(inflowlabels[i]), bbox=dict(facecolor='red', alpha=.5))  # commented out for now so as to be able to see links when zoom in

###    fig.text(.3, .8, 'Testing', size=10, transform=ax.transAxes)
##    fig.text(.3, .8, 'Testing', size=10)
    endtime = time.time()
    timetodrawgraph = endtime - starttime
    print(f"Drawing of graph took {timetodrawgraph:.5f} seconds.")
    plt.show()

  return outputGeoJSON


#############################################################################
if __name__ == '__main__':  # Needed only if another file shall import this one as a module
#  DEBUG = False  # flag whether to print debugging-messages

  SEMversion = 'SEM2'  # use the basin-hopping global-optimisation solver
#  SEMversion = 'SEM3'  # find shortest-path routes: this is (see thesis by Gao, 2009) a static (link-volumes, travel-costs, and origin-destination demands are constant over time), deterministic (each vehicle takes with probability one its shortest path, rather than each vehicle having an estimate for each path's cost that differs from the actual cost by a random error), equilibrium (the travel-costs on all paths used from any origin to any destination are equal, while all unused paths have equal or greater costs) method.
  print("__name__ == '__main__': SEMversion is", SEMversion)

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

  outputGeoJSON = runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, exitnodes)
#  outputGeoJSON = runSEM( SEMversion, JSONnetworkfilename, inflowsByNodeID, exitnodes, simulDuration=2.0)
