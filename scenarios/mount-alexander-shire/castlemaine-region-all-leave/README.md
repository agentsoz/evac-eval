# Evaluating bushfire evacuation models

The scenarios contained in this directory represent the population (all-leave) that strictly follows the message to leave and evacuate.
This population is specified in the Matsim plan file "population-archetypes.xml.gz".
This scenario is available in [scenarios/mount-alexander-shire/castlemaine-region](scenarios/mount-alexander-shire/castlemaine-region).
In this scenario, which is based on the EES [castlemaine-region-archetypes](https://github.com/agentsoz/ees/tree/master/scenarios/mount-alexander-shire/castlemaine-region-archetypes) scenario, the resident population of Castlemaine region (technically the 2016 Australian census population in statistical areas *Castlemaine SA2* and *Castlemaine Region SA2*) is evacuated in response to a Catastrophic (Code Red) fire to the northwest of Maldon. Details of the original scenario including a video of the simulation run are available [here](https://github.com/agentsoz/ees/blob/master/scenarios/mount-alexander-shire/castlemaine-region-archetypes/README.md#commit-495a7f8--9-sep-2019).

The evacuation message is given in the file "scenario_messages.json.

The network used is given in the file "loddon_mallee_northern_cluster_shires_network.xml.gz".

The fire hazard  is given in the file "20181109_mountalex_evac_ffdi100d_grid.json".

 
## Scenarios
There are two scenarios described in this directory.

### 1. Castlemaine region all-leave evacuation

This scenario simulates the simple evacuation of the (all-leave) population according to the evacuation message.
The EES config file is "archetypes.xml".
The Matsim config file is "archetypes_matsim_main.xml".

### 2. Castlemaine region all-leave evacuation  with road closure

This scenario simulates the simple evacuation of the (all-leave) population according to the evacuation message (as in Scenario 1).
However, a single vital road is closed 4 hours into the evacuation.
The road closure is described in the change events file "networkChangeEvents.xml"
The EES config file for this scenario is "archetypes_rc.xml".
The Matsim config file is "archetypes_matsim_main_rc.xml".
