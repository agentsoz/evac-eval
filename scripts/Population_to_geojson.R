library(rgdal)

#ca.pop <- read.csv(data, header = TRUE, sep = ",", quote = "\"")
ca.pop <- read.csv(gzfile("population-archetypes.csv.gz"), as.is = TRUE)

ca.pop$Geographical.Coordinate <- as.character(ca.pop$Geographical.Coordinate)

# split coordinates into separate columns
coords <- strsplit(substr(ca.pop$Geographical.Coordinate, 2, nchar(ca.pop$Geographical.Coordinate) - 1), ", ")

ca.pop$Geographical.Coordinate.trimmed <- substr(ca.pop$Geographical.Coordinate, 2, nchar(ca.pop$Geographical.Coordinate) - 1)

ca.pop.coords <- data.frame(do.call("rbind", strsplit(ca.pop$Geographical.Coordinate.trimmed, ", ", fixed = TRUE)))

ca.pop$Geographical.Coordinate <- NULL
ca.pop$Geographical.Coordinate.trimmed <- NULL

str(ca.pop)

ca.pop$Long <- as.numeric(as.character(ca.pop.coords$X1))
ca.pop$Lat <- as.numeric(as.character(ca.pop.coords$X2))

#ca.pop$Long <- as.numeric(ca.pop.coords$X1)
#ca.pop$Lat <- as.numeric(ca.pop.coords$X2)

# convert to Spatial*DataFrame
ca.pop.SP <- SpatialPointsDataFrame(ca.pop[, c(37, 38)], ca.pop[, -c(37, 38)], proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Write as geojson
writeOGR(ca.pop.SP, "castlemaine-area_pop.json", "castlemaine-area-population", driver = "GeoJSON", overwrite_layer = TRUE)

# write as shapefile
#writeOGR(ca.pop.SP, "castlemaine_area_pop.shp", "castlemaine-area-population", driver = "ESRI Shapefile", overwrite_layer = TRUE)
