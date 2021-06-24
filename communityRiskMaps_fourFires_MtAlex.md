# Community risk per population-node, as posed by four fires (Mount Alexander Shire)

![Diagram: fires ffdi100{a,b,c,d}, Mount Alexander Shire](fourOverlappingFiresMtAlex.png)

![Diagram: raster-layer showing how many fire-perimeters each cell lies inside](rasterLayer-numFirePerimetersInside-fourfires.png)

![Diagram: with population archetypes superimposed, defined by 'home' activity](populationArchetypesSuperimposed-fourfires.png)

![Diagram: with population archetypes as coarse-grained raster](populationArchetypesSuperimposed-coarseGrainedRaster-fourfires.png)

![Diagram: network vector-layer with each population-node that's inside at least one fire, its raster-cell coloured by the number of fires it sits inside](populationNodesAndRasterOfNumFiresPopulationIsInside-fourfires.png)



# Each population raster-cell inside at least one fire is assigned an injection-node - the network node having largest maximum out-capacity

![Diagram: network vector-layer with links coloured by capacity and population-node coloured by largest maximum out-capacity](linksColouredByCapacity_populationNodesColouredByLargestMaxOutCapacity-fourfires.png)


# Maximum-flow method run on fire ffdi100a

![Diagram: maximum-flow applied to fire ffdi100a, Mount Alexander Shire](maxFlowFireA.png)

![Diagram: maximum-flow applied to fire ffdi100a, in close-up](maxFlowFireA_closeUp.png)

![Diagram: maximum-flow applied to fire ffdi100a, links coloured by flow](maxFlowFireA_linksColouredByFlow.png)
