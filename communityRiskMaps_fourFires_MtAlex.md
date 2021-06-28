# Community risk per population node, as posed by four fires (Mount Alexander Shire)

![Diagram: fires ffdi100{a,b,c,d}, Mount Alexander Shire](fourOverlappingFiresMtAlex.png)

![Diagram: raster-layer showing how many fire-perimeters each cell lies inside](rasterLayer-numFirePerimetersInside-fourfires.png)

### Population archetypes superimposed:
![Diagram: with population archetypes superimposed, defined by 'home' activity](populationArchetypesSuperimposed-fourfires.png)

### As a coarser-grained raster:
![Diagram: with population archetypes as coarse-grained raster](populationArchetypesSuperimposed-coarseGrainedRaster-1000-fourfires.png)

### Number of fires each population-cell sits inside (community risk)
![Diagram: Number of fires each population-cell sits inside (community risk)](numFiresPopulationIsInside-1000-fourfires.png)



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

### Define exit-nodes as network-nodes that are both in the "safe" buffer-area and close to population?
