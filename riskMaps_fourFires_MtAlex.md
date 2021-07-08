# Risks posed by four fires (Mount Alexander Shire)

![Diagram: fires ffdi100{a,b,c,d}, Mount Alexander Shire](fourOverlappingFiresMtAlex.png)

![Diagram: raster-layer showing how many fire-perimeters each cell lies inside](rasterLayer-numFirePerimetersInside-fourfires.png)

### Population archetypes superimposed, as a coarser-grained raster:
![Diagram: with population archetypes as coarse-grained raster](populationArchetypesSuperimposed-coarseGrainedRaster-1000-fourfires.png)



# Maximum-flow method run on fire ffdi100a

### Inside the fire-perimeter (road-links coloured by capacity):
![Diagram: inside perimeter of fire ffdi100a](insideFirePerimeter-ffdi100a-1000.png)

### Non-zero population inside the fire:
![Diagram: non-zero population inside perimeter of fire ffdi100a](populationInsideFirePerimeter-ffdi100a-1000.png)

## Each population raster-cell with non-zero population and inside the fire is assigned a "population node" (injection-node) - the network node having largest maximum out-capacity

### Population nodes coloured by largest maximum out-capacity:
![Diagram: network vector-layer with population nodes coloured by largest maximum out-capacity](linksColouredByCapacity_populationNodesColouredByLargestMaxOutCapacity-1000-fourfires.png)

### "Safe" area surrounding fire:
![Diagram: "safe" area surrounding fire](safeAreaSurroundingFire-1000-fourfires.png)

### Define exit-nodes as network-nodes lying in the "safe" buffer-area:
![Diagram: exit-nodes in "safe" area](exitNodesInSafeBuffer-fourfires.png)


### Maximum-flow solution:
![Diagram: maximum-flow solution](maxFlowSoln_linksColouredByFlow_populationNodesColouredByEvacuationTime-fourfires.png)

### Number of evacuation-flows each road-link is critical to:


## Community risk per population node (mean time taken to evacuate node, over all fires that contain it within their perimeters):
![Diagram: community risk per population node](communityRiskPerPopulationNode-1000-fourfires.png)


## Ignition risk per fire (number of population nodes plus number of critical links within fire-perimeter):
![Diagram: Ignition risk per fire](ignitionRiskPerFire-1000-fourfires.png)


## Point-of-impact risk per link (number of evacuation-flows the link is critical to, over all fires):
![Diagram: Point-of-impact risk per link](pointOfImpactRiskPerLink-1000-fourfires.png)
