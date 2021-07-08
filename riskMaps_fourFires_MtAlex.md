# Risks posed by four fires (Mount Alexander Shire)

![Diagram: fires ffdi100{a,b,c,d}, Mount Alexander Shire](fourOverlappingFiresMtAlex.png)

### Population archetypes superimposed:
![Diagram: with population archetypes superimposed](populationArchetypesSuperimposedOnFires-fourfires.png)



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

### Define exit-nodes as network-nodes lying in the "safe" buffer-area (thickness of a road-link indicates capacity):
![Diagram: exit-nodes in "safe" area](exitNodesInSafeBuffer-fourfires.png)


### Maximum-flow solution (thickness of a green link indicates flow; population nodes coloured by evaucation-time):
![Diagram: maximum-flow solution](maxFlowSoln_linksColouredByFlow_populationNodesColouredByEvacuationTime-fourfires.png)



# The three types of evacuation risk, measured across all four fires

## Community risk per population node (mean time taken to evacuate node, over all fires that contain it within their perimeters):
![Diagram: community risk per population node](communityRiskPerPopulationNode-1000-fourfires.png)


## Ignition risk per fire (number of population nodes plus number of critical links within fire-perimeter):
![Diagram: Ignition risk per fire](ignitionRiskPerFire-1000-fourfires.png)


## Point-of-impact risk per link (number of evacuation-flows the link is critical to, over all fires):
![Diagram: Point-of-impact risk per link](pointOfImpactRiskPerLink-1000-fourfires.png)
