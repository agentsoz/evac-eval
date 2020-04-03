
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


## How to convert Phoenix fire model Shapefile to GeoJSON

**Important**: requires `GDAL2.4`

Example command to convert the fire Shapefile `ffdi100d_grid.shp` to GeoJSON format:
```
ogr2ogr -f "GeoJson" -t_srs EPSG:28355 ffdi100d_grid.json ffdi100d_grid.shp

```

Example Phoenix fire shapefiles for the regions are available in the `ees-data` repository:
* Mount Alexander Shire: https://github.com/agentsoz/ees-data/tree/master/mount-alexander-shire/phoenix-shapefiles/20181109
* Surf Coast Shire: https://github.com/agentsoz/ees-data/tree/master/surf-coast-shire/phoenix
