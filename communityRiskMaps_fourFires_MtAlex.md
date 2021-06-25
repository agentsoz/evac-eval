# Community risk per population node, as posed by four fires (Mount Alexander Shire)

![Diagram: fires ffdi100{a,b,c,d}, Mount Alexander Shire](fourOverlappingFiresMtAlex.png)

![Diagram: raster-layer showing how many fire-perimeters each cell lies inside](rasterLayer-numFirePerimetersInside-fourfires.png)

### Population archetypes superimposed:
![Diagram: with population archetypes superimposed, defined by 'home' activity](populationArchetypesSuperimposed-fourfires.png)

### Population as a coarser-grained raster:
![Diagram: with population archetypes as coarse-grained raster](populationArchetypesSuperimposed-coarseGrainedRaster-fourfires.png)

### Number of fires each population-cell sits inside (community risk)
![Diagram: Number of fires each population-cell sits inside (community risk)](numFiresPopulationIsInside-fourfires.png)



# Maximum-flow method run on fire ffdi100a

## Each population raster-cell with non-zero population and inside the fire is assigned a "population node" (injection-node) - the network node having largest maximum out-capacity

### Links coloured by capacity and population nodes coloured by largest maximum out-capacity:
![Diagram: network vector-layer with links coloured by capacity and population nodes coloured by largest maximum out-capacity](linksColouredByCapacity_populationNodesColouredByLargestMaxOutCapacity-fourfires.png)

### Each population node that's inside at least one fire, its raster-cell coloured by the number of fires it sits inside:
![Diagram: network vector-layer with each population node that's inside at least one fire, its raster-cell coloured by the number of fires it sits inside](populationNodesAndRasterOfNumFiresPopulationIsInside-fourfires.png)

![Diagram: maximum-flow applied to fire ffdi100a, Mount Alexander Shire](maxFlowFireA.png)

![Diagram: maximum-flow applied to fire ffdi100a, in close-up](maxFlowFireA_closeUp.png)

![Diagram: maximum-flow applied to fire ffdi100a, links coloured by flow](maxFlowFireA_linksColouredByFlow.png)
