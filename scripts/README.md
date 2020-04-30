
## How to convert MATSim network XML to GeoJSON

**Important**: the following script requires `Java11`, `Python3.6`, and `GDAL 2.4`.

The script downloads all required files, including MATSim itself, as well as the regional MATSim XML network files, and converts them all to `GeoJSON` format. To run, do:
```
./createExampleMATSimJsonNetworks.sh
```

If all went well, you should have the following files created in your current directory:
```
./loddon_mallee_northern_cluster_shires_network/*json
./surf_coast_shire_network/*.json
./mount_alexander_shire_network/*.json

```

**NOTE**: Currently the output files are generated in current directory but these files needs to be moved into `/scenarios/mount-alexander-shire/castlemaine-region/gssem`. I am working on this.

## How to convert Phoenix fire model Shapefile to GeoTiff

**Important**: requires `GDAL2.4`, `Python3.6`.

The script that converts the Phoenix fire model shapefile into Geotiff is `shp_to_geotiff_ffdi.py`. This script will take `20181109_mountalex_evac_ffdi100d_grid.shp` from `/scenarios/mount-alexander-shire/castlemaine-region/` as input shape file, convert into GeoTiff and save the output GeoTiff file into `/scenarios/mount-alexander-shire/castlemaine-region/gssem` as `20181109_mountalex_evac_ffdi100d_grid.tif`.

Example Phoenix fire shapefiles for the regions are available in the `ees-data` repository:
* Mount Alexander Shire: https://github.com/agentsoz/ees-data/tree/master/mount-alexander-shire/phoenix-shapefiles/20181109
* Surf Coast Shire: https://github.com/agentsoz/ees-data/tree/master/surf-coast-shire/phoenix

## How to convert ABS 2016 Synthetic Population CSV to GeoTiff

**Important**: requires `GDAL2.4`, `Python3.6`, `R3.6`.

The full synthetic population for Greater Melbourne for the 2016 census is available here: https://github.com/agentsoz/synthetic-population/tree/master/data.

The script used to convert the `population-archetypes.csv.gz` into GeoJSON is `Population_to_geojson.R`. To run this script from commandline go to the /script directory and type the following command `Rscript Population_to_geojson.R` 

Once the GeoJSON file is converted then following command is used to convert the `castlemaine-area_pop.json` GeoJSON file to shape file.

```
ogr2ogr -nlt POINT -skipfailures castlemaine-area_pop.shp castlemaine-area_pop.json OGRGeoJSON  

```

Once the shape file is generated then use the script `shp_to_tiff_pop.py` to convert the shape file into GeoTiff


**Note**: Full automation of the script is in progress.

## How to create MATSim Population XML from ABS 2016 Synthetic Population

The following script takes the archetypes population CSV-that is based on the ABS 2016 Synthetic Population-and converts it to the MATSim XML file. To run, do
```
./createMATSimPopnFromArchetypesCsv.sh
```

**TODO**: To create a *dumber* population file for the EES model we need to change the relevant columns in the CSV and generate simpler versions of the population using this script. For instance, to disable the behaviour of going to check on dependants we can just set the column `HasDependents` to `false`. If we set all attributes for every person to the same value (other than their XY locations), then they will all display the same behaviour.
