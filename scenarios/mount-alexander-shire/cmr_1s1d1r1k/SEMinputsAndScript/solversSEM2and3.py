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
#EPS = 1e-3  # this shouldn't be a global variable, so its definition has been moved inside each of the functions 'adjustSpeedandDensityToLink()' and 'flowandspeedEnteringNode()'


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


def geoJSONtoAdjacencyMatrix( JSONnetworkfilename ):
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
  exitnodes = []
  for l in linestringsWithCoords:  # add any missing Points at the ends of LineStrings
    for xy in l:
      if xy not in pointsWithCoords:
        print(f"Endpoint {xy} of LineString {l} is not in pointsWithCoords: adding it")
        pointsWithCoords[ xy ] = currentPointID
        exitnodes.append( currentPointID )
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
def setuptrafficflowproblem( JSONnetworkfilename ):
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

  (adjacency, exitnodes, pointcoordsbyID, linestringsWithCoords, lK, numlanes, linklength, linkcapacity, matsimlinkID, oneway) = geoJSONtoAdjacencyMatrix( JSONnetworkfilename )
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
    maximumdensityperlane = 120.0  # in vehicles per kilometre per lane; a constant over all edges - this value taken from p. 93, 'Traffic Flow Dynamics' (2013) by Treiber and Kesting
    print(f"maximumdensityperlane is {maximumdensityperlane:.4f}")
#    T = 1.4 / 3600  # in hours; a constant over all edges - this value taken from p. 93 of Treiber and Kesting
#    print(f"T is {T:.9f}")
    T = {}
    for edge in lK.keys():
      T[edge] = numlanes[edge] / linkcapacity[edge] - 1/(lK[edge]*maximumdensityperlane)  # in hours; from eq. (8.15), p. 93 of Treiber and Kesting
#      rhoC[edge] = 1 / (lK[edge] * T + 1/maximumdensityperlane)  # density per lane at capacity; from p. 94 of Treiber and Kesting
      tau[edge] = numlanes[edge] / (lK[edge] * T[edge] + 1/maximumdensityperlane)  # density over entire link-width at capacity; from p. 94 of Treiber and Kesting
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
#    print("flow is %.4f, congestiononlink is %.4f" % (flow, congestiononlink))
  else:
    congestiononlink = 0.0
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
        print("=========== WARNING: THE SOLUTION FAILS at iteration %d: np.linalg.LinAlgError exception raised. ========================" % iteration)
        return (e, None, None, False, iteration)  # e, hInjectionNodes, qij, flowbalanced
      print("dh is %s" % dh)
      h += dh  # if h falls outside the bounds if any, e.g. |h_i| <= 4 for flowfunc = 'concaveQuadratic', should h be adjusted to fall within the bounds? If so, should h be placed at the boundary or should it be "reflected" back within the bounds? (Probably the former, I imagine, as I don't know a way to justify "reflection," which might anyway cause the adjusted h to fall outside the bounds in the opposite direction.) But I believe that h will fall outside the bounds only if there is no root within the bounds for it to find, so adjusting h would merely delay the information that no root has been found.
      try:
#        dh_dot = scipy.linalg.norm( dh )  # incorrect?
        dh_dot = scipy.linalg.norm( dh )**2
      except ValueError as e:  # e.g., dh contains 'inf' or 'NaN'
        print("=========== WARNING: THE SOLUTION FAILS at iteration %d: ValueError exception raised, perhaps by 'inf' or 'NaN'. ========================" % iteration)
        return (e, None, None, False, iteration)  # e, hInjectionNodes, qij, flowbalanced
      except OverflowError as e:  # e.g., (34, 'Numerical result out of range')
        print("=========== WARNING: THE SOLUTION FAILS at iteration %d: OverflowError exception raised, perhaps because numerical result is out of range. ========================" % iteration)
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
      print("=========== WARNING: Newton-Raphson method has not converged: has reached maximum number (%d) of iterations without ||dh||^2 (%.12f) falling to %.4e * ||dh0||^2 (||dh0||^2 = %.5f); ||dh||^2/||dh0||^2 = %.6e" % (MAXITERS, dh_dot, eps_NRconvergence, dh0_dot, normsratio))
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
      print("=========== WARNING: THE SOLUTION FAILS at iteration %d: flow does not balance over network. ========================" % iteration)
    else:
      print("=========== WARNING: THE SOLUTION FAILS: flow does not balance over network. ========================")
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
  if densityinlink >= humpdensity:
    print(f"WARNING: density in link {{{p},{j}}} of {densityinlink:.3f} is over the hump {humpdensity:.3f}.")
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


def flowandspeedEnteringNode( j, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnExit, flowijfunc, flowfunc ):
  # Traffic injected at a node flows towards the next node on its shortest path towards an exit. At the beginning of the first link the flow and speed are known, which gives the density; for traffic flowing into the link this density gives, according to v = v(rho) for the link, the speed in the link. The flow is conserved between each injection-node and the end of the first link on its shortest path to an exit; if a link is divided for numerical reasons into sub-sections, then the flow is conserved between adjacent sections.
  EPS = 1e-3
#  print(f"Entering flowandspeedEnteringNode( {j}, parentsOnShortestPathsFromInjectionNodes )")
  print(f"Entering flowandspeedEnteringNode( {j}, ... ):")
  if j in exitnodes:
    print(f"Node {j} is an exit-node, so has zero injection-flow and -density, and injection-speed is undefined")
#    print(f"{j} is an exit-node, so has zero injection-flow and injection-speed is undefined.")
    injectionflow = injectiondensity = 0.0
    injectionspeed = 'undefined'
    print(f"At node {j} injectionflow is {injectionflow:.3f}, injectionspeed is {injectionspeed}, and injectiondensity is {injectiondensity:.3f}.")
  else:
    injectionflow = b0[indexamonginjectionorexitnodes[j]]
    print(f"injectionflow is {injectionflow:.5f}")
    injectionspeed = v0[indexamonginjectionorexitnodes[j]]
#    print(f"injectionspeed is {injectionspeed:.5f}")
    print(f"injectionspeed is {injectionspeed}")
    injectiondensity = injectionflow / injectionspeed
#    print(f"injectiondensity is {injectiondensity:.5f}")
    print(f"injectiondensity is {injectiondensity}")
    print(f"At node {j} injectionflow is {injectionflow:.3f} and injectionspeed is {injectionspeed:.3f}, so injectiondensity is {injectiondensity:.3f}.")

  if parentsOnShortestPathsFromInjectionNodes[j] == set():
    print(f"Node {j} has no parent-nodes on shortest paths from injection-nodes.")  # N.B. node j might be an exit-node with no flow entering it
    if j in exitnodes:
      print(f"Node {j} is an exit-node, so should have zero injection-flow and injection-speed")
      assert injectionflow == 0.0
      assert injectionspeed == 0.0
      print(f"Node {j} is an exit-node; returning injection-flow {injectionflow:.3f} and injection-speed {injectionspeed:.3f}")
      return injectionflow, injectionspeed
    print(f"Node {j} is an injection-node; returning injection-flow {injectionflow:.3f} and injection-speed {injectionspeed:.3f}")
    return injectionflow, injectionspeed

  # Node j has at least one parent-node on shortest paths from injection-nodes:
  totalflow = injectionflow  # total flow so far: flows from parent-nodes will be added to this
  totaldensity = injectiondensity  # vehicles per longitudinal metre of road, no matter how many lanes it has laterally
  speedsinlinks = []  # store speeds in all links flowing into j
  print(f"Node {j}'s parent-nodes on shortest paths from injection-nodes are {parentsOnShortestPathsFromInjectionNodes[j]}")
  for p in parentsOnShortestPathsFromInjectionNodes[j]:
    flowenteringparent, speedenteringparent = flowandspeedEnteringNode( p, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnExit, flowijfunc, flowfunc )
    densityenteringparent = flowenteringparent / speedenteringparent
#    print(f"Node {p} has closest exit-node {closestExit[p]}, at distance {shortestDistanceToAnExit[p]:.4f} with shortest path {shortestPathToAnExit[p]}")
#    print(f"p is {p}, j is {j}")
    assert j == shortestPathToAnExit[p][1]
    linkcapacity = argsflowfunc['linkcapacity'][frozenset({p,j})]
    if flowenteringparent > linkcapacity:
      # TODO: is it physically realistic to have traffic-flow in excess of the link-capacity use a different route to an exit, perhaps a second-shortest path? Note there might be several different paths all with the shortest length.
      print(f"ALERT: entering-flow of {flowenteringparent:.1f} is too high for link ({p}, {j}) which has a capacity of {linkcapacity:.3f}; high congestion and density, and low speed, will result.")
    print(f"The total inflow at node {p} of {flowenteringparent:.3f} veh/hr starts to flow towards node {j} with density {densityenteringparent:.3f} veh/km and speed v0={speedenteringparent:.3f} km/hr.")
    flowfunction = flowijfunc[flowfunc]
    (speedinlink, densityinlink, numiterations) = adjustSpeedandDensityToLink( flowfunction, p, j, densityenteringparent, argsflowfunc, flowenteringparent, speedenteringparent )
    print(f"After {numiterations} iterations, the speed adjusts to {speedinlink:.3f} km/hr and density to {densityinlink:.3f}.")
    qij[p, j] = flowenteringparent
    densityij[frozenset({p, j})] = densityinlink
    speedij[frozenset({p, j})] = speedinlink
    if flowenteringparent > linkcapacity:
#      assert abs(speedinlink) < 1e-1  # all the precision that's required
      assert abs(speedinlink) < 1e-30  # true in tests so far
#    humpdensity = argsflowfunc['tau'][frozenset({p,j})]
#    if densityinlink > humpdensity:
#      congestiononlink = linkcapacity - abs(flowenteringparent)
    congestiononlink = linkcongestion(p, j, flowenteringparent, densityinlink, argsflowfunc)
    congestionij[frozenset({p, j})] = congestiononlink
    if abs(congestiononlink) > 0.0:  # TODO: is the abs() necessary?
      print(f"CONGESTION of {linkcapacity:.3f}-{flowenteringparent:.1f}={congestiononlink:.3f} occurs on link ({p}, {j}) which has capacity of {linkcapacity:.3f}.")  # TODO: this measure of congestion gives a negative value when the flow entering the parent is greater than the link-capacity: might need a better measure
    totalflow += flowenteringparent
    totaldensity += densityinlink
    speedsinlinks.append( speedinlink )
    print(f"Adding {flowenteringparent:.3f} to totalflow and {densityinlink:.3f} to totaldensity that enter node {j}.")
#  speed = totalflow / totaldensity
  if abs(totaldensity) >= EPS:
    speed = totalflow / totaldensity  # v = q/rho
  elif abs(totalflow) < EPS:
    speed = np.mean( speedsinlinks )  # both total flow and total density at node j are close to zero, so we set speed to the mean value of the speeds in inbound links (FIXME: a bit of a hack, although might well be fine in practice)
  else:  # flow in link is far from zero, density in link is close to zero
    speed = totalflow / totaldensity  # v = q/rho
  print(f"At node {j}, totalflow is {totalflow:.3f} and totaldensity is {totaldensity:.3f}, so speed is {speed:.3f}.")
  return totalflow, speed


def runSEM2orSEM3( SEMversion, JSONnetworkfilename, inflowsByNodeID):
  starttime = time.time()
  print("inflowsByNodeID is", inflowsByNodeID)
##  (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, totalinflow, argsflowfunc, lK, linklength, v0, pointcoordsbyID) = setuptrafficflowproblem( JSONnetworkfilename, b0 )
#  (terminatorPressure, eps, eps_flowbalance, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, argsflowfunc, pointcoordsbyID) = setuptrafficflowproblem( JSONnetworkfilename )
  (terminatorPressure, eps, flowfunc, flowijfunc, derivijfunc, adjacency, N, exitnodes, injectionnodes, indexamonginjectionorexitnodes, neighboursofj, hExitNodes, argsflowfunc, pointcoordsbyID, linestringsWithCoords) = setuptrafficflowproblem( JSONnetworkfilename )

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
  SPEEDOFINJECTEDTRAFFIC = 20  # assumed speed of injected traffic at every injection-node (speed of traffic as it enters an injection-node), in kilometres per hour
#  SPEEDOFINJECTEDTRAFFIC = 30  # assumed speed of injected traffic at every injection-node (speed of traffic as it enters an injection-node), in kilometres per hour

  v0 = np.full( shape=len(b0), fill_value=SPEEDOFINJECTEDTRAFFIC)  # speed of traffic as it enters each injection-node, in kilometres per hour
  print("v0 (speed of traffic as it enters an injection-node) is", v0)
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


#  else:  # don't use the solver; instead find shortest-path routes
  elif SEMversion == 'SEM3':  # don't use the solver; instead find shortest-path routes
    starttime = time.time()
    print("\nFinding shortest-path routes from all injection-nodes ...")
    vertices = range(N)
    predictedLinkTraversalTime = {}  # traversal-time as might be predicted by a driver
    edges = []  # lK stores the 'free speed' of each road-link
    lK = argsflowfunc['lK']
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

    for i in injectionnodes:  # proceed down each injection-node's shortest path to an exit, adding for each node after the first one on the path its parent on the path as one of that node's parents on all shortest paths from injection-nodes
      for index, j in enumerate( shortestPathToAnExit[i] ):
        if index > 0:
          parentsOnShortestPathsFromInjectionNodes[j].add( shortestPathToAnExit[i][index-1] )
    print(f"shortestDistanceToAnExit is", shortestDistanceToAnExit)
    print("closestExit is", closestExit)
    print("shortestPathToAnExit is", shortestPathToAnExit)
    print("parentsOnShortestPathsFromInjectionNodes is", parentsOnShortestPathsFromInjectionNodes)
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

    qij = {}  # store flows on all edges
    congestionij = {}  # store congestion-values on all edges
    densityij = {}  # store densities on all edges
    speedij = {}  # store densities on all edges
    outflowsByNodeID = {}  # store the out-flows at exit-nodes, as they are calculated by the solution-method (Newton-Raphson, global optimiser, or shortest-path-based method)
    for j in exitnodes:
      flow, speed = flowandspeedEnteringNode( j, parentsOnShortestPathsFromInjectionNodes, qij, congestionij, densityij, speedij, argsflowfunc, exitnodes, b0, indexamonginjectionorexitnodes, v0, shortestPathToAnExit, flowijfunc, flowfunc )
      outflowsByNodeID[j] = -flow
    print(f"flow is {flow:.3f}, speed is {speed:.3f}.")
    endtime = time.time()
    print("densityij is", densityij)
    print("speedij is", speedij)
    timetorunshortestpathbasedmethod = endtime - starttime
    print(f"Shortest-path-based method took {timetorunshortestpathbasedmethod:.5f} seconds.")



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
    else:
      if SEMversion == 'SEM2':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -qij[b, a]
      elif SEMversion == 'SEM3':  # found shortest-path routes, rather than using the global optimiser
        linkJSON["properties"]["flow"] = -densityij[frozenset({b, a})] * speedij[frozenset({b, a})]
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
    if SEMversion == 'SEM3':  # found shortest-path routes, rather than using the global optimiser
      print("densityij is %s" % densityij)
      print("speedij is %s" % speedij)
    G = nx.DiGraph()  #G = nx.Graph()
#    DECIMALPLACES = 4  # number of decimal places to show in values on the graph
#    DECIMALPLACES = 3  # number of decimal places to show in values on the graph
    DECIMALPLACES = 1  # number of decimal places to show in values on the graph
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
      elif SEMversion == 'SEM3':
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
    print("BBeast-BBwest is", BBeast-BBwest)
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
      print("np.array([x, y] is", np.array([x, y]))
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
  networktype = 'linear'

  print(f"JSONnetworkfilename is {JSONnetworkfilename}")

#  outputGeoJSON = runSEM2orSEM3( SEMversion, inputfilename, inflowsByNodeID)
  outputGeoJSON = runSEM2orSEM3( SEMversion, JSONnetworkfilename, inflowsByNodeID)
