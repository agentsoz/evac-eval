# Risks posed by four fires (Mount Alexander Shire)

![Diagram: fires ffdi100{a,b,c,d}, Mount Alexander Shire](fourOverlappingFiresMtAlex.png)

### Population archetypes superimposed:
![Diagram: with population archetypes superimposed](populationArchetypesSuperimposedOnFires-fourfires.png)



# The maximum-flow method run on fire ffdi100a (the south-east fire)

### Inside the fire-perimeter (road-links coloured by capacity, e.g. Calder Freeway in yellow, Midland Highway in red):
![Diagram: inside perimeter of fire ffdi100a](insideFirePerimeter-ffdi100a-1000.png)

### Non-zero population inside the fire:
![Diagram: non-zero population inside perimeter of fire ffdi100a](populationInsideFirePerimeter-ffdi100a-1000.png)

## Each population raster-cell with non-zero population and inside the fire is assigned a population node (injection-node), i.e. the network node having largest maximum out-capacity

### Population nodes coloured by largest maximum out-capacity (yellow high, red low):
![Diagram: network vector-layer with population nodes coloured by largest maximum out-capacity](linksColouredByCapacity_populationNodesColouredByLargestMaxOutCapacity-1000-fourfires.png)

### "Safe" buffer-area surrounding fire, between one and two kilometres from fire-perimeter:
![Diagram: "safe" buffer-area surrounding fire](safeAreaSurroundingFire-1000-fourfires.png)

### Define exit-nodes as network-nodes lying in the buffer-area (thickness of a road-link indicates capacity):
![Diagram: exit-nodes in "safe" buffer-area](exitNodesInSafeBuffer-fourfires.png)


### Maximum-flow solution (thickness of a green link indicates flow; population nodes coloured by evacuation-time, defined to be population/inflow):
Problem: each white node has zero inflow in the maximum-flow solution, so its evacuation-time is infinite. This arises from there being effectively unlimited supply of flow at each injection-node (because even a tiny population could enter the network with a flow of 4000 vehicles per hour), so that the inflow at one node can "flood" a road's (even freeway's) capacity and prevent other nodes' populations from flowing onto that road; this follows from the max-flow model's assumption that each inflow persists for all time. The problem is reduced or goes away altogether, depending on population sizes, if we set the inflow at each node to a small multiple of its population, but that might be too unrealistic.

One approach is to estimate the time it would take for the "flooding" traffic (i.e. from nodes with positive inflows) to exit the network, and then run the max-flow method again with those nodes' inflows now set to zero; this would be repeated until no population-node had an inflow of zero in the max-flow solution. The sum of estimated exit-times would then be an upper bound on the time for all nodes' populations to have reached the exit-nodes.

The literature has variations on the maximum-flow problem, e.g. 'Maximum-Flow Evacution-Planning Problem with Non-conservation Flow Constaint' (Bhandari and Khadka,2020), some of which might be useful.
![Diagram: maximum-flow solution](maxFlowSoln_linksColouredByFlow_populationNodesColouredByEvacuationTime-fourfires.png)



# The three types of evacuation risk, measured across all four fires

## Community risk per population node (mean of the node's evacuation-times, i.e. population/inflow, over all fires that contain the node within their perimeters) - a white node has zero inflow in the maximum-flow solution for at least one fire, so its mean evacuation-time is infinite:
![Diagram: community risk per population node](communityRiskPerPopulationNode-1000-fourfires.png)


## Ignition risk per fire (number of population nodes plus number of critical links within fire-perimeter) - where fires overlap, the ignition risks are added (should we take the maximum instead?):
![Diagram: Ignition risk per fire](ignitionRiskPerFire-1000-fourfires.png)


## Point-of-impact risk per link (number of evacuation-flows the link is critical to, over all fires - black 1, red 2, yellow 3):
![Diagram: Point-of-impact risk per link](pointOfImpactRiskPerLink-1000-fourfires.png)
