
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
