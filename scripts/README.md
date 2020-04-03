
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

**TODO**: The produced GeoJSON files then need be be converted to GeoTIFF format as required by GSSEM.

## How to convert Phoenix fire model Shapefile to GeoJSON

**Important**: requires `GDAL2.4`

Example command to convert the fire Shapefile `ffdi100d_grid.shp` to GeoJSON format:
```
ogr2ogr -f "GeoJson" -t_srs EPSG:28355 ffdi100d_grid.json ffdi100d_grid.shp

```

Example Phoenix fire shapefiles for the regions are available in the `ees-data` repository:
* Mount Alexander Shire: https://github.com/agentsoz/ees-data/tree/master/mount-alexander-shire/phoenix-shapefiles/20181109
* Surf Coast Shire: https://github.com/agentsoz/ees-data/tree/master/surf-coast-shire/phoenix

**TODO**: The produced GeoJSON files then need be be converted to GeoTIFF format as required by GSSEM.


## How to convert ABS 2016 Synthetic Population CSV to GeoTiff

The full synthetic population for Greater Melbourne for the 2016 census is available here: https://github.com/agentsoz/synthetic-population/tree/master/data.

The household files for a given region can be converted to the GeoTIFF format required by GSSEM using the script here: https://bitbucket.csiro.au/users/for321/repos/emv2/browse/data-munging/synthetic-population-csv-to-vector.R

**TODO**: The script should be copied here and adjusted as needed.

## How to create MATSim Population XML from ABS 2016 Synthetic Population

The following script takes the archetypes population CSV-that is based on the ABS 2016 Synthetic Population-and converts it to the MATSim XML file. To run, do
```
./createMATSimPopnFromArchetypesCsv.sh
```

**TODO**: To create a *dumber* population file for the EES model we need to change the relevant columns in the CSV and generate simpler versions of the population using this script. For instance, to disable the behaviour of going to check on dependants we can just set the column `HasDependents` to `false`. If we set all attributes for every person to the same value (other than their XY locations), then they will all display the same behaviour.
