library(dplyr)
library(sp)
library(cism)
library(rgeos)
library(geosphere)
library(readxl)
if('2016-12-07_HOUSEHOLD.RData' %in% dir('data')){
  load('data/2016-12-07_HOUSEHOLD.RData')
} else {
  HOUSEHOLD <- get_data(dbname = 'MALTEM',
                tab = 'HOUSEHOLD')
  save(HOUSEHOLD,
       file = 'data/2016-12-07_HOUSEHOLD.RData')
}

# Read in the bairros names
dict <- read_excel('public_data/lista total de agregados por bairro.xlsx',
                   skip = 1)

# Get a bairro number
dict$bairro_number <- 
  as.numeric(unlist(lapply(dict$`Cod do Bairro`,
                           function(x){
                             substr(x, 3, 4)
                           })))

# Clean up names and narrow
dict <- dict %>%
  dplyr::select(bairro_number,
                `posto Administrativo`,
                Bairros) %>%
  rename(posto_administrativo = `posto Administrativo`,
         bairro_name = Bairros)

# Get shorter name
df <- HOUSEHOLD

# Specify coordinates
df$lng <- df$HOUSEHOLD_HEAD_GPS_LNG
df$lat <- df$HOUSEHOLD_HEAD_GPS_LAT

# Get neighborhood number
df$bairro_number <-
  unlist(lapply(df$HOUSEHOLD_HEAD_AGREG_NUM,
                function(x){
                  substr(x, 3,4)
                }))
df$bairro_number <- as.numeric(as.character(df$bairro_number))

# Narrow data
df <- df %>%
  dplyr::select(bairro_number,
                lng,
                lat) %>%
  mutate(x = lng,
         y = lat) %>%
  filter(!is.na(lng),
         !is.na(lat))

# Make spatial
coordinates(df) <- ~x+y
proj4string(df) <- proj4string(mag2)

df@data$id <- 1:nrow(df)
df@data$suspect <- FALSE

# Get unique bairros
bairros <- data_frame(bairro_number = sort(unique(df$bairro_number)))

# Remove form our household data, those points which 
# are suspiciously far from others
for (i in 1:nrow(bairros)){
  # Get info just for this bairro
  this_bairro <- bairros$bairro_number[i]
  sub_df <- df[df@data$bairro_number == this_bairro,]
  # Get distances between points
  distances <- geosphere::distm(sub_df, fun = distVincentySphere)
  median_distances <- apply(distances, 1, function(x){quantile(x, 0.75, na.rm = TRUE)})
  # Keep only the closest 50 percent
  close_cluster <- sub_df[which(median_distances <= median(median_distances, na.rm = TRUE)),]
  # Flag those which are suspicious
  far_cluster <- sub_df[which(median_distances > median(median_distances, na.rm = TRUE)),]
  df$suspect[df$id %in% far_cluster$id] <- TRUE
}

# Remove those which are suspicious
df <- df[!df@data$suspect,]

# Get centroids
bairros$x <- bairros$y <- NA
for (i in 1:nrow(bairros)){
  # Get info just for this bairro
  this_bairro <- bairros$bairro_number[i]
  sub_df <- df[df@data$bairro_number == this_bairro,]
  bairros$x[i] <- mean(sub_df$lng, na.rm = TRUE)
  bairros$y[i] <- mean(sub_df$lat, na.rm = TRUE)
}

# Create convex hulls
ch <- list()
for (i in 1:nrow(bairros)){
  this_bairro <- bairros$bairro_number[i]
  sub_df <- df[df@data$bairro_number == this_bairro,]
  x <- rgeos::gConvexHull(sub_df)
  ch[[i]] <- x
}

# Create delaunay triangulation / voronoi tiles for entire surface
voronoi <- function(shp = df){
  
  shp@data <- data.frame(shp@data)
  
  # Fix row names
  row.names(shp) <- 1:nrow(shp)
  
  # Remove any identical ones
  shp <- shp[!duplicated(shp$lng,
                                                     shp$lat),]
  
  # Helper function to create coronoi polygons (tesselation, not delaunay triangles)
  # http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/
  voronoipolygons = function(layer) {
    require(deldir)
    crds = layer@coords
    z = deldir(crds[,1], crds[,2])
    w = tile.list(z)
    polys = vector(mode='list', length=length(w))
    require(sp)
    for (i in seq(along=polys)) {
      pcrds = cbind(w[[i]]$x, w[[i]]$y)
      pcrds = rbind(pcrds, pcrds[1,])
      polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
    }
    SP = SpatialPolygons(polys)
    voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1], 
                                                           y=crds[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                        function(x) slot(x, 'ID'))))
  }
  # http://gis.stackexchange.com/questions/180682/merge-a-list-of-spatial-polygon-objects-in-r
  appendSpatialPolygons <- function(x) {
    ## loop over list of polygons
    for (i in 2:length(x)) {
      # create initial output polygon
      if (i == 2) {
        out <- maptools::spRbind(x[[i-1]], x[[i]])
        # append all following polygons to output polygon  
      } else {
        out <- maptools::spRbind(out, x[[i]])
      }
    }
    return(out)
  }
  
  tile_polys <- voronoipolygons(shp)
  # Add the bairro numbers
  tile_polys@data$bairro_number <- the_bairros <- shp$bairro_number
  cols <- rainbow(as.numeric(factor(tile_polys@data$bairro_number)))
  
  # Disolve borders
  x = gUnaryUnion(tile_polys, id = tile_polys$bairro_number)
  
  jdata = SpatialPolygonsDataFrame(Sr=x, 
                                   data=data.frame(bairro_number = as.numeric(as.character(names(x)))),FALSE)
  
  return(jdata)
}

# Get voronoi tesselations
dfv <- voronoi(shp = df)

# Plot in current state
plot(dfv)
plot(mag2, add = T, col = adjustcolor('blue', alpha.f = 0.2))

# Narrow down so as to only keep those areas which are IN Magude
proj4string(dfv) <- proj4string(mag2)
out <- gIntersection(dfv, mag2, byid=TRUE)

# Join with data
row.names(out) <- as.character(1:length(out))
out <- SpatialPolygonsDataFrame(out, data.frame(bairros), match.ID = TRUE)

# Plot again
plot(mag2, col = adjustcolor('red', alpha.f = 0.2))
plot(out, add = T, lwd = 0.2)

# Join in info on bairros
dict$full_name <- paste0(dict$posto_administrativo,
                         ': ',
                         dict$bairro_name)
out@data <- left_join(out@data,
                      dict,
                      by = 'bairro_number')
bairros <- left_join(bairros,
                     dict, 
                     by = 'bairro_number')


# Plot with ggplot2 framework and labels
library(ggplot2)
library(maptools)
library(ggrepel)
row.names(out) <- as.character(1:nrow(out@data))

outgg <- fortify(out, region = 'bairro_number')

ggplot() +
  geom_polygon(data = outgg,
               aes(x = long,
                   y = lat,
                   group = group),
               color = 'black',
               fill = 'white') +
  coord_map() +
  geom_label_repel(data = bairros,
             aes(x = x,
                 y = y,
                 label = full_name),
             size = 1,
             alpha = 0.5,
             label.padding = unit(0.05, "lines"),
             box.padding = unit(0.05, "lines"))

# Plot leaflet
library(leaflet)
leaflet() %>%
  addProviderTiles(provider = 'Stamen.Toner') %>%
  addPolygons(data = out,
              color = 'green',
              popup = out@data$bairro_name) 
