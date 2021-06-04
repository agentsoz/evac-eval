# Overview of static-evacuation models (SEM)

## Inputs for a scenario
* A network of links (edges) and nodes.
* Each link has a 
  * flow-capacity, 
  * free speed, 
  * length, 
  * number of road-lanes, and 
  * flag to indicate whether the link is one-way or both directions can be used.
* A subset of *injection-nodes,* at each of which there is a specified in-flow of traffic from nearby homes.
* A subset of *exit-nodes* where traffic can flow out of the network.
* Every node is either an injection-node or an exit-node, and not both.
* The in-flow at an injection-node is sub-divided into specified subflows, one for each exit-node that has been assigned to that injection-node.

## Comparison of versions of the model

Model | Why currently fails | Advantages | Drawbacks | Improvements (if any) needed to make model work
------ |------------------- |---------- |--------- |------------------------
SEM1  | Each link has infinite flow-capacity, meaning congestion can't arise | Simple to implement; Newton-Raphson method quickly solves model | Congestion can't arise in model, owing to no upper limit on flow through a link |
SEM2  | Traffic can flow out at any exit-node: can't send a subflow to its assigned exit-node. On at least one network having a solution, model can find no solution. | Quite simple to implement | Global optimiser is slow, and must often be run several times from different starting-vectors | Instead of one variable h per node, could try one variable (density) per arc, but without additional constraints the assignment to each arc of free or congested flow would be arbitrary
SEM3  | Implementation is incomplete | Identifies some congested links; sends subflows to their assigned exit-nodes | Upstream propagation of congestion is perhaps too simplified; model and implementation are complex | Extension to correctly propagate shockfronts through nodes, with resulting downstream effects
SEM4  | Implementation is incomplete | Identifies some congested links; sends subflows to their assigned exit-nodes; reports time- and link-dependent extent of upstream-propagating congestion | Model and implementation are complex | Extension to correctly propagate shockfronts through nodes, with resulting downstream effects
SEM5  | If capacities allow, maximum-flow method correctly sends total assigned flow from each injection-node and total assigned flow to each exit-node; but doesn't satisfy individual subflows between origin/destination pairs | Maximum-flow is implemented in libraries, and code is fast | Doesn't consider individual subflows between injection/exit pairs, so fails to detect all congestion along subflows (reports false negatives) | Apply maximum-flow method to the sub-network containing only links used by assigned subflows

## Versions of the model in more detail

### SEM1 (base version) – hydraulic-based flow through pipes
* Flow-conservation constraint applies at each non-exit node: sum(inflows) = sum(outflows).
* Flow is a monotonic function of head-loss, meaning multivariate Newton-Raphson method can be used to find pressure-head values that imply flows obeying the constraint.
* Flow-capacity in each link is infinite: can obtain an arbitrarily high flow by making the head-loss sufficiently large. This implies congestion never arises.

### SEM2 – hydraulic-based flow through pipes, with varying fluid-density
* Based on SEM1, but a link's head-loss is interpreted as its density: implies flow is density-dependent.
* Flow is a non-monotonic function of density, with a "hump": densities below the hump-density correspond to freely-flowing traffic, densities above it to congested traffic.
* Non-monotonicity of flow means that Newton-Raphson fails: instead must use a global optimiser to find head-values that obey the flow-conservation constraint.
* Each link now has a flow-capacity (the upper bound on its flow-function).

* The five-node network illustrated below, with all links having equal free speed and number of lanes, and assigned subflows of 
  - 100 from node 1 to node 0
  - 500 from node 2 to node 0
  - 500 from node 2 to node 4

has the unique solution shown (flow-capacities in square brackets). But with one h-variable per node, SEM2 cannot find this solution, even if we allow arbitrary h-values at exit-nodes: two links with flows of 500 implies h2-h0 = h2-h3, but h2-h0 = (h2-h3)+(h3-h0) > h2-h3, a contradiction.
![Diagram: five-node network with no SEM2 solution](fivenode-nopossibleSEM2soln-triangularflowfn-100-1000.png)
* Could seek to improve model by replacing h-variables at nodes with one density-variable per arc, but then the model wouldn't be able to choose on each link between low (free) and high (congested) densities, so there would be at least 2^|E| distinct solutions and the model wouldn't reliably report congestion: would need additional constraints to enforce relationships between adjacent links' densities.

### SEM3 - traffic-flow along road-links: congestion results from insufficient capacity downstream
* Models traffic flowing downstream and encountering a link of insufficient capacity, so that links upstream of the bottleneck can be identified as congested.
* Traffic entering at an injection-node takes a shortest-path route to each of its assigned exit-nodes.
* Theory taken from Chapter 8 of Treiber & Ketting, *Traffic Flow Dynamics* (2013).
* Piece-wise linear 'triangular' flow-function gives speed and density within each link: has two pieces, corresponding to free and congested flow.
![Diagram: triangular fundamental diagram](triangularFundamentalDiagram.png)
* Traffic-flow is divided into distinct ''zones'', each with constant values for flow, density, and speed; the zones' boundaries are wavefronts that propagate with time. A link can contain several zones, or just one.
* Free flow propagates downstream as a wavefront, with propagation-speed equal to the free speed of its current link.
* If a free-flow wavefront reaches a link with insufficient capacity to carry its flow, then congested flow is assumed to back up, i.e. propagate upstream, all the way to the injection-nodes.
![Diagram: congestion propagates upstream](congestionPropagatesUpstream.png)
* Once all downstream-propagating wavefronts have departed the network, the resulting flows on links give a static solution.

### SEM4 - traffic-flow along road-links, with time-dependent upstream-propagating congestion
* Essentially a dynamic model, as the wavefronts' positions are updated through time. The simulation has a specified duration (in hours), which if set sufficiently long gives a static solution.
* Congested flow backs up through time as a ''shockfront'' (upstream-propagating wavefront).
* A shockfront can catch up with a slower-propagating shockfront, or meet a downstream-propagating wavefront, on the same link.
* The next time-step is defined as the shortest elapsed time until a wavefront reaches the end of its link or intersects with another wavefront on that link.
* Can indicate differing extents of congestion in links by reporting each link's time-averaged flow, i.e. the cumulative number of vehicles (in general a floating-point value) to have departed the link during the simulation, divided by the simulation's duration.

### SEM5 - detect congestion by finding the network's maximum flow
* Reports the maximum-possible flow through network, using supplies and demands derived from the assigned subflows from injection- to exit-nodes; reports links that attain their flow-capacity at maximum flow.
* If the maximum flow is less than total of assigned subflows, then congestion on at least some links is unavoidable.
* Maximum-flow method doesn't consider individual subflows between origin/destination pairs, so congestion along subflow-paths is possible even when the maximum flow equals total of assigned subflows; consider the network below (flow-capacities in square brackets), with assigned subflows of 
  - 800 from node 0 to node 2
  - 200 from node 0 to node 3.
When the subflows use only links on shortest paths, which exclude arc (2,3), there is congestion on arc (0,1):
![Diagram: four-node network with congestion along assigned subflows despite sufficiently-high maximum flow found by maximum-flow method](fournode-congestionOnAssignedSubflowDespiteSufficientMaxFlow-800-200.png)
Maximum-flow method finds maximum flow (1000) equal to the total of assigned subflows, but this flow is achieved only by using arc (2,3) - a false negative:
![Diagram: four-node network on which maximum-flow method fails to find congestion along assigned subflows](fournode-maxFlowMethodFailsToFindCongestionOnAssignedSubflow-800-200.png)
* This problem might be alleviated by calculating maximum flow on the sub-network that contains only links carrying assigned subflows.


## Next steps for development
* Complete implementation of SEM3/4 to correctly propagate congested flow upstream.
* Perhaps optimise SEM3/4 to obtain an improved time-complexity, which is currently as high as O(E N^4).


## Authors of software-implementations
* **James Hilton** - *SEM1*
* **Stephen Taylor** - *Generalisation of SEM1 to SEM2, implementation of SEM3/4 and SEM5*
