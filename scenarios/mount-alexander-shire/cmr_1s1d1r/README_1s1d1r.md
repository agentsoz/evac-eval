# Scenario "Castlemaine Region - 1 source/ 1 destination / 1 road"

The label "cmr_1s1d1r"  refers to the scenario "Castlemaine Region - 1 source/ 1 destination / 1 road" which has been created to provide comparison with similar outputs from the SEM model.

CMR_1S1D1R is a version of the "castlemaine-region-all-leave" scenario with the following modifications.:

### 1. The population of 10867 inthe plans XML (cmr_1s1d_plans.xml) consists of identical persons.
### 2. All persons leave from the same source (252696.7593, 5909967.483).
### 3. All persons travel to the same destination (262886.1854, 5890668.996).
### 4. The road network consists of one long road given in cmr_1s1d1r_network.xml. 

The main files in the scenario include:
* cmr_1s1d1r_ees.xml - this is the EES configuration file to be included in the java call to execute. 
* cmr_1s1d1r_config.xml - this is the Matsim configuration file. 
* cmr_1s1d1r_network.xml - this is the Matsim network configuration file. 
* cmr_1s1d_plans.xml - this is the Matsim plan configuration file. 
* cmr_1s1d_plans.csv - this is the CSV version of the Matsim plan configuration file. 

