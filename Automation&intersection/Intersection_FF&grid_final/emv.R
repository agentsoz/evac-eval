#loading libraries
library(geojsonsf)
library(sf)
library(dplyr)
library(tidyverse)
library(lubridate)
library(dbscan)
library(caret)
library(sp)
library(units)
library(rio)
library(ndjson)
library(jsonlite)
library(magrittr)
library(stringr)
library(plotly)
library(geojsonio)


# +
#Function for EPSG:4326 CRS
crs_epsg4326 <- function() {
    crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
    return(crs)}
    
crs_epsg32754 <- function() {
    crs="+proj=utm +zone=54 +south +datum=WGS84 +units=m +no_defs"
    return(crs)}

# +
#Loading the fire file and the grid file
gridfile <- geojson_sf("grid.geojson")

firefile <- geojson_sf("Final_fireFile.json")

# -

#Transforming the grid file
gridfile %>% 
    select(gridID,geometry) %>%
    mutate(gridID=as.character(gridID)) %>%
    st_transform(crs=crs_epsg32754())->grid


#Checking the structure of the data.
grid %>% as.data.frame() %>% head(5)

# Replace hour burnt NA to zeros
firefile %>% select(geometry, HL_PROB)->ffile


#Set crs
ffile = st_set_crs(ffile, 32754)

#Structure of the data
ffile %>% as.data.frame() %>% head(5)

# set as Sf and change crs
ffile %>% 
    #st_as_sfc(crs = st_crs(grid))%>%
    #st_sf(ffile = c('x','y'), geoms = ., stringsAsFactors = FALSE)
    st_as_sf(coords = c("Y_COORD","X_COORD"),crs=crs_epsg32754())%>%
    st_transform(crs=crs_epsg32754())->ffile_final

ffile_final_area <- mutate(ffile_final, ffile_final_area = st_area(ffile_final))

ffile_final_area %>% as.data.frame()
#Intersection
intersection = as.tibble(st_intersection(grid,ffile_final))

plot(grid$geometry, axes = TRUE)
plot(ffile_final$geometry, add = TRUE, col = "white")
plot(intersection$geometry, add = TRUE, col = 'red')

#Calculate the area
intersection$area <- st_area(intersection$geometry)

#Groupby
intersection %>% 
  as_tibble() %>% 
  group_by(HOUR_BURNT, gridID) %>% 
  summarize(intersect_area = sum(area)) -> final_area_int # werite a function to calculate a cumulitive sum 

final_area_int %>% as.data.frame() %>% head(15)
final_area_int %>% group_by(gridID) %>%  mutate(cum_area = cumsum(intersect_area)) %>% ungroup() %>% 
  select(HOUR_BURNT, gridID, cum_area) ->filter_grid

#filter_grid%>% as.data.frame()
#final_area_int %>% group_by(gridID)%>% summarize(a = sum(intersect_area)) -> final_a

filter_grid %>% as.data.frame()

#Calculating the grid cell area
grid_area <- mutate(grid, grid_area = st_area(grid))
grid_area %>% as.data.frame() %>% head()

#Merging two tables
grid_final <- merge(grid_area, filter_grid, by = "gridID", all.y = TRUE)

# Calculate the covergae area %
grid_final <- grid_final %>% 
   mutate(area_coverage = as.numeric(cum_area/grid_area))

grid_final %>% as.data.frame() %>% head()

#writing the file.
inter <- sf_geojson(grid_final)
geojson_write(inter, file = "grid_intersection_cumulative.geojson")

# # END


