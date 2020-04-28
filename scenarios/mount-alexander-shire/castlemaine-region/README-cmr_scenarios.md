## Running the test Castlemaine Region (CMR) scenarios
There are four test scenarios.

1. CMR0 - This scenario is based on the original scenario (population-archetypes.xml).

2. CMR1 -  This scenario uses the same population attributes as in CMR0 
        -  The same plan activities (type="home") as in CMR0
        -  But provides 3 evacuation destinations (instead of 1)
        -  AND issues 1 EVACUATE_NOW message at 1300HRS (instead of multiple messages) 

3. CMR2 -  This scenario uses the same population attributes as in CMR0 
        -  But provides New "other"(evacuation-type) activities at 1500HRS
        -  But provides 3 evacuation destinations (instead of 1)
        -  AND issues 1 EVACUATE_NOW message at 1300HRS (instead of multiple messages) 

4. CMR3 -  This scenario uses the "dumbed down" population attributes as in CMR0 
        -  But provides New "other"(evacuation-type) activities at 1500HRS
        -  But provides 3 evacuation destinations (instead of 1)
        -  AND issues 1 EVACUATE_NOW message at 1300HRS (instead of multiple messages) 

IN Windows:
1. To run EES Jar for scenario CMR0 use:
>>> java -Xms2g -Xmx2g -cp libs\*;ees-2.1.1-SNAPSHOT.jar io.github.agentsoz.ees.Run --config scenarios\mount-alexander-shire\castlemaine-region\cmr0_ees.xml

2. To run EES Jar for scenario CMR1 use:
>>> java -Xms2g -Xmx2g -cp libs\*;ees-2.1.1-SNAPSHOT.jar io.github.agentsoz.ees.Run --config scenarios\mount-alexander-shire\castlemaine-region\cmr1_ees.xml

3. To run EES Jar for scenario CMR2 use:
>>> java -Xms2g -Xmx2g -cp libs\*;ees-2.1.1-SNAPSHOT.jar io.github.agentsoz.ees.Run --config scenarios\mount-alexander-shire\castlemaine-region\cmr2_ees.xml

4. To run EES Jar for scenario CMR3 use:
>>> java -Xms2g -Xmx2g -cp libs\*;ees-2.1.1-SNAPSHOT.jar io.github.agentsoz.ees.Run --config scenarios\mount-alexander-shire\castlemaine-region\cmr3_ees.xml

