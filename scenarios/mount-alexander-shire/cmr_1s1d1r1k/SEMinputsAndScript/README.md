# Input-files and Python script to run the SEM on the simple scenario cmr_1s1d1r1k (outputs to be stored on CloudStor)

The script `solveFlow-cmr_1s1d1r1k.py` runs the static-evacuation model (SEM) on the simple scenario cmr_1s1d1r1k ("Castlemaine region, one source, one destination, one road, one thousand vehicles"), for purposes of comparison with the agent-based model (ABM). The script is based on `Vivian.py` (created in 2020?). It specifies a traffic in-flow at the northernmost node of a linear sequence of road-links, then runs the SEM and writes the results (flows on the links) to the file `resultsSEM-cmr_1s1d1r1k.geojson`.

## Prerequisites
It might be best and easiest to create a virtual environment, and use conda to install into it the Python bindings to Geostack and any other required packages.
Follow the instructions on https://gitlab.com/geostack/library/-/wikis/Installing-Geostack-for-Python-using-conda to install miniconda for Python 3.x.

### Commands (a possible sequence in Linux Bash)
```
conda config --set auto_activate_base False
conda config --set auto_update_conda False
conda create --name gsenv python=3.7  # create an environment called 'gsenv' containing Python 3.7
conda activate gsenv  # activate the environment 'gsenv'
conda config --add channels conda-forge
conda config --append channels geostack  # need channel 'geostack' in order to install bindings ('geostack')
# Before the following command, might need to install an OpenCL implementation with 'conda install pocl' or similar
conda install --channel geostack geostack 
git clone https://gitlab.com/geostack-applications/applications.git  # install networkFlow.py (solver)
# [A sequence of the user's commands] 
conda deactivate  # deactivate the current environment ('gsenv') if desired - not strictly necessary
```

### Installing
Install and run the script by entering a sequence of commands such as
```
git clone https://github.com/agentsoz/evac-eval
cd evac-eval/scenarios/mount-alexander-shire/cmr_1s1d1r1k/SEMinputsAndScript/
python solveFlow-cmr_1s1d1r1k.py
```

### Next steps for development
* Allow the user to specify the input-file and the output-file via options on the command-line.
* 

## Authors
* **Vivian Dabre** - *Original script*
* **Stephen Taylor** - *Modification and generalisation* 

## License
The script relies on the Geostack library and Python bindings to it, released under the Mozilla Public Licence, version 2.0 (https://choosealicense.com/licenses/mpl-2.0/).

## Acknowledgments
* Dhirendra Singh, James Hilton, Leorey Marquez, and others provided guidance.
