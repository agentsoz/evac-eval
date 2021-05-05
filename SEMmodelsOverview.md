# Overview of static-evacuation models (SEM)

## Inputs for a scenario
* A network of nodes and edges
* Each edge has a flow-capacity, free speed, and number of lanes
* A subset of ''injection-nodes,'' at each of which there is a specified in-flow of traffic
* A subset of ''exit-nodes'' where the traffic can flow out of the network.

## Versions of the model

### SEM1 (base version) – hydraulic-based flow through pipes
* Properties of SEM (base version) – model of hydraulic flow through pipes
* Constraint: flow is conserved at every node: sum(inflows) = sum(outflows).
* Flow is a monotonic function of head-loss, meaning multivariate Newton-Raphson method can be used to find head-values that obey the constraint.
* Flow is proportional to pipe-diameter, but pipe-diameter doesn't constrain capacity: can obtain an arbitrarily high flow as a result of sufficiently large head-loss. So no link has an upper limit on fluid-speed or -flow.

### SEM2 – hydraulic-based flow through pipes, varying fluid-density
* Based on SEM1, but a link's head-loss is interpreted as its density: implies that flow is density-dependent.
* Flow is a non-monotonic function of density, with a "hump": densities below the hump-density correspond to freely-flowing traffic, densities above it to congested traffic.
* Non-monotonicity of flow means that Newton-Raphson fails: instead must use a global optimiser to find head-values that obey the flow-conservation constraint.
* Each link now has a flow-capacity (the upper bound on its flow-function).
* Link-lengths play no role, so traffic can take unlikely routes to the exit-nodes: this makes false negatives possible, when the model doesn't find a highly likely route that would give rise to congestion.
- The simple three-node network with one injection-node can, besides a false positive as previously established, also provide a false negative: if the righthand link is very long relative to the lefthand link, then no traffic in reality would choose to travel down it, yet the model predicts that densities on the two links are equal. Hence if the lefthand link is not sufficiently capacious to take all inflowing traffic without congestion, yet both links together can take all traffic without congestion, then the model will report a solution with no congestion even though all traffic travels to the left and becomes congested.
![title](threenodefalsenegative-triangularflowfn-1000.png)

## SEM3 - traffic-flow along road-links, with flow dependent on density
* Use piece-wise linear ''triangular'' flow-function: has two pieces, corresponding to free and congested flow, respectively (and is non-monotonic in density).
* Traffic entering at any injection-node takes its shortest-path route to the nearest exit-node (yet to be generalised to the case in which some traffic-flow heads for one exit-node while other flow heads for another).
* Traffic-flow is divided into distinct ''zones'', each with its own values for flow, density, and speed; the zones' boundaries are wavefronts that propagate with time. A link can contain several zones, or just one.
* Free flow propagates downstream as a wavefront, with propagation-speed equal to the free speed of its current link.
* If a free-flow wavefront reaches a link with insufficient capacity to carry its entire flow, then congested flow propagates upstream as a ''shockwave'' (upsteram-propagating wavefront). 
* Every congestion-zone propagates upstream.
* Piece-wise linear triangular flow-function gives speed and density within each link.


The script `solveFlow-cmr_1s1d1r1k.py` runs the static-evacuation model (SEM) on the simple scenario cmr_1s1d1r1k ("Castlemaine region, one source, one destination, one road, one thousand vehicles"), for purposes of comparison with the agent-based model (ABM). It specifies traffic in-flows at certain nodes (the ''injection-nodes'') of a linear sequence of road-links, then runs the SEM and writes the results (flows on the links) to the file `resultsSEM-cmr_1s1d1r1k.geojson`.


# Unused notes
SEM1:
*hydraulic-based flow model*
SEM1 assumes constant fluid-density?
· Two-way pipes for water vs. one-way links for traffic
. Pipe diameter vs. road capacity: 
. Flow-density relationship (water vs cars) - the model implies that density is constant
* Flow is a monotonic function of head-loss (implying constant fluid-density?)
SEM2:
* Hazen-Williams flow-equation is used.
* Non-monotonic flow captures density-dependence of traffic-speed.
* I thought the flow-conservation constraints might allow all flows to be free (not congested), but that's not true.
SEM3
* Use Lighthill-Whitham-Richards model to calculate flows in links.
* Means a low-capacity link constrains distant upstream links: unrealistic, or implicit in any static-flow model?
* A shockwave propagates (upstream) with lower speed than a (downstream-propagating) wavefront - but is this always the case?

# Scenarios
* QGIS map of SEM3 link-flow for one hour
SEM3 flow for scenario CMR_1s1d1r1k
* Test scenario – 1 source, 1 destination, 1 route (50 links), 1000 vehicles



### Commands (a possible sequence in Linux Bash)
```
conda config --set auto_update_conda False
conda deactivate  # deactivate the current environment ('gsenv') if desired - not strictly necessary
```


### Next steps for development
* Allow the user to specify the input-file and the output-file via options on the command-line.
* 

## Authors
* **James Hilton** - *SEM1*
* **Stephen Taylor** - *Generalisation of SEM1 to SEM2, implementation of SEM3 and SEM4* 

## License
The script relies on the Geostack library and Python bindings to it, released under the Mozilla Public Licence, version 2.0 (https://choosealicense.com/licenses/mpl-2.0/).
