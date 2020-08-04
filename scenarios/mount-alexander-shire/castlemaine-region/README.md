## Castlemaine Region scenario

This scenario is based on the [original scenario here](https://github.com/agentsoz/ees/tree/master/scenarios/mount-alexander-shire/castlemaine-region-archetypes) (see [video](https://cloudstor.aarnet.edu.au/plus/s/phuCz9Q3suj9slb)). The colour scheme in the video is as follows:

Cars | Meaning |
--- | --- |
Blue | Worried Waverers |
Yellow | Considered Evacuators |
Green | Community Guided |
Pink | Responsibility Deniers |
Orange | Threat Deniers |
Red | Experienced Independents|

Dots | Meaning |
--- | --- |
Blue | People at home |
Pink | People deliberating/preparing (lag time before driving)
Purple | People at Invac/Evac place |
White | People at Dependents place |

The [final messaging regime](./scenario_messages.json ) is as follows:

Message | Time sent | Zones sent to |
--- | --- | --- |
`ADVICE`  |  `1110hrs` | [western SA1s](MountAlexander_SA1s/west.geojson)
`WATCH_AND_ACT` | `1200hrs` | [western SA1s](MountAlexander_SA1s/west.geojson)
`EMERGENCY_WARNING` | `1230hrs` | [western SA1s](MountAlexander_SA1s/west.geojson)
`EVACUATE_NOW` | `1300hrs` | [western SA1s](MountAlexander_SA1s/west.geojson)
`ADVICE` | `1330hrs` | [inner SA1s](MountAlexander_SA1s/inner.geojson)
`EMERGENCY_WARNING` | `1500hrs` | [inner west SA1s](MountAlexander_SA1s/inner-west.geojson)
`WATCH_AND_ACT` | `1600hrs` | [inner SA1s](MountAlexander_SA1s/inner.geojson)

*Note that the EES model (actually the underlying BDI-ABM action/percept model) does not support sending two messages (percepts) of the same type at the same time to an agent. Therefore the timing of the different messages should be unique.*


## Phoenix fire GeoJSON

File `20181109_mountalex_evac_ffdi100d_grid.shp` is the original shape file used as input file to convert the shp file to GeoTIFF file.

File `20181109_mountalex_evac_ffdi100d_grid.json` was generated using:
```
ogr2ogr -f "GeoJson" -t_srs EPSG:28355 \
  20181109_mountalex_evac_ffdi100d_grid.json \
  ../../../../ees-data/mount-alexander-shire/phoenix-shapefiles/20181109/Evac_Phoenix_runs/20181109_mountalex_evac_ffdi100d/20181109_mountalex_evac_ffdi100d_grid.shp

```
