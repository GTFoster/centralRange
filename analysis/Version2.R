## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------------------------------------------
library(geometry)
library(sp)
library(igraph)
library(dplyr)
library(ggplot2)
library(gridExtra)
#library(rgbif)
library(curl)
#library(BIEN)
library(maps)
library(sf)
sf_use_s2(FALSE)
library(terra)
library(maps)
library(data.table)
library(geodata)
#library(rangeBuilder) #For creating alpha hulls #Edit: This is a good approach, but takes way too long. 
#library(alphahull)
#Erase later-only for debugging
library(tictoc)
#library(exactextractr)
library(geosphere)
#library(hypervolume)
library(raster)

## ---------------------------------------------------------------------------------------------------------------------------------------
load("data/MetaWebCentrality.Rda") #Load in Global Metaweb
metaCent$biennm <- gsub(pattern=" ", replacement="_",   metaCent$node) #Make name convention same a bien
Clima <- geodata::worldclim_global("worldclim", var="bio", res=5) #Load in global climate data (just as a mask)
wadmin <- geodata::world(resolution = 5, level =0, path="data/climate") #Download world administrative boundaries
ccs <- geodata::country_codes() #get country codes-includes continent key
wadminrast <- terra::rasterize(wadmin, Clima, field="GID_0") #rasterize so it's faster to join

fls <- list.files(path="data/GBIF/occs/") #Grab all the filenames in our occurence data dir
nms <- gsub("_", " ", fls) %>% gsub(".csv", "", .) #Refit to match data
nms <- paste(toupper(substr(nms, 1,1)), substr(nms,2,100), sep="")
#missing <- metaCent[metaCent$node %in% nms==FALSE,] #Extract those that aren't in our dat

dat <- metaCent[metaCent$node %in% nms==TRUE,] #Nodes which occur in our dataset
dat$flnm <- paste(gsub(" ", "_", tolower(dat$node)),".csv", sep="") #add the filenames of their occ points

temp <- dat[595,]
filtemp <- paste("data/GBIF/occs/", temp$flnm, sep="") #set up temp filename
tempdat <- read.csv(file=filtemp) #read in file

tempdat$cells <- terra::cellFromXY(Clima, xy=tempdat[,5:4]) #Grab the cell ids
indcells <- tempdat %>% dplyr::group_by(., cells) %>% dplyr::slice_sample(n = 1) #remove redundant cells (randomly choose 1 point)
indcells$ISO3 <- wadminrast[indcells$cells]$GID_0 #grab country ID
indcells <- ccs %>% dplyr::select(., ISO3, continent) %>% left_join(indcells, ., by=c("ISO3")) #Grab continent

points <- dplyr::filter(indcells, continent=="South America")


## ---------------------------------------------------------------------------------------------------------------------------------------
#' @points A dataframe with lat-long points in the 5th and 4th column, respectively
#' @threshold quantile value for the amount of data points to remove, starting with those farthest from the centroid first. 
#' 
#' @output a spatial polygon of MCH
#' 
# Function to calculate the convex hull area
cut_alphahull <- function(points, threshold=.95) {
  hull.pol<-chull(points[,5:4]) #Calculate alpha-hull for full data
  coords <- points[c(hull.pol, hull.pol[1]), ]
  PolOur <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords[,5:4])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
  center <- sf::st_centroid(PolOur)
  ptsob <- sf::st_as_sf(points[,5:4], coords=c("decimalLongitude", "decimalLatitude"), crs=crs(PolOur))
  dists <- sf::st_distance(center, ptsob)
  under <- dists<quantile(dists,threshold)
  cut <- points[under[1,]==TRUE,]
  
  hull.pol2<-chull(cut[,5:4]) #Calculate alpha-hull for full data
  coords2 <- cut[c(hull.pol2, hull.pol2[1]), ]
  PolOur2 <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords2[,5:4])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
  return(PolOur2)
}

maps::map('world', fill = TRUE, col = "white")
plot(cut_alphahull(points=points, threshold=0.99), col="red", add=T, alpha=0.5)
plot(cut_alphahull(points=points, threshold=0.7), col = "blue", add=T)

output <- cut_alphahull(points=points, threshold=0.95)
area <- st_area(output)
units(area)$numerator <- c("km", "km") # Get Niche Area


## ---------------------------------------------------------------------------------------------------------------------------------------
load(file="data/climate/env.RData")


## ---------------------------------------------------------------------------------------------------------------------------------------
points$envCell <- terra::cellFromXY(env, xy=as.matrix(points[,5:4]))


indcells_env <- points %>% dplyr::group_by(., envCell) %>% dplyr::slice_sample(n = 1) #remove redundant cells (randomly choose 1 point)

vals <- terra::extract(env, as.matrix(indcells_env[,5:4]))
indcells_env <- cbind(indcells_env, vals)


## ---------------------------------------------------------------------------------------------------------------------------------------
cut_alphahull_env <- function(input, threshold=.95) {
  input <- input[is.finite(input$fpcat.PC1)==T,]
  hull.pol<-chull(input[,10:11]) #Calculate alpha-hull for full data
  coords <- input[c(hull.pol, hull.pol[1]), ]
  PolOur <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords[,10:11])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
  center <- sf::st_centroid(PolOur)
  ptsob <- sf::st_as_sf(input[,10:11], coords=c("fpcat.PC1", "fpcat.PC2"), crs=crs(PolOur))
  dists <- sf::st_distance(center, ptsob)
  under <- dists<quantile(dists,threshold)
  cut <- input[under[1,]==TRUE,]
  
  hull.pol2<-chull(cut[,10:11]) #Calculate alpha-hull for the cut down data
  coords2 <- cut[c(hull.pol2, hull.pol2[1]), ]
  PolOur2 <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords2[,10:11])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
  return(PolOur2)
}
samples <- sample(1:nrow(env.df), 10000, replace=F)

plot(x=env.df$fpcat.PC1[samples], y=env.df$fpcat.PC2[samples], col="grey", xlab="PC1", ylab="PC2")
#plot(cut_alphahull_env(input), col="blue", add=T)
#points(x=input$fpcat.PC1, y=input$fpcat.PC2, col="red")

#EnvOutput <- cut_alphahull_env(input=input, threshold=0.95)
#area <- st_area(EnvOutput)
#units(area)$numerator <- c("km", "km") # Get Niche Area


## ---------------------------------------------------------------------------------------------------------------------------------------
makefuncs <- function(){
  library(geometry)
  library(sp)
  library(igraph)
  library(dplyr)
  library(ggplot2)
 # library(gridExtra)
 # library(rgbif)
  library(curl)
 # library(BIEN)
  library(maps)
  library(sf)
  sf_use_s2(FALSE)
  library(terra)
  library(maps)
  library(data.table)
  library(geodata)
  #library(rangeBuilder) #For creating alpha hulls #Edit: This is a good approach, but takes way too long. 
  #library(alphahull)
  #Erase later-only for debugging
  library(tictoc)
 # library(exactextractr)
  library(geosphere)
 # library(hypervolume)
  library(raster)  

  Clima <- unwrap(packed_Clima)
  
  #' @points A dataframe with lat-long points in the 5th and 4th column, respectively
  #' @threshold quantile value for the amount of data points to remove, starting with those farthest from the centroid first. 
  #' 
  #' @output a spatial polygon of MCH
  #' 
  # Function to calculate the convex hull area
  cut_alphahull <- function(points, threshold=.95) {
    hull.pol<-chull(points[,5:4]) #Calculate alpha-hull for full data
    coords <- points[c(hull.pol, hull.pol[1]), ]
    PolOur <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords[,5:4])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
    center <- sf::st_centroid(PolOur)
    ptsob <- sf::st_as_sf(points[,5:4], coords=c("decimalLongitude", "decimalLatitude"), crs=crs(PolOur))
    dists <- sf::st_distance(center, ptsob)
    under <- dists<quantile(dists,threshold)
    cut <- points[under[1,]==TRUE,]
    
    hull.pol2<-chull(cut[,5:4]) #Calculate alpha-hull for full data
    coords2 <- cut[c(hull.pol2, hull.pol2[1]), ]
    PolOur2 <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords2[,5:4])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
    return(PolOur2)
  }
  
  #' @points A dataframe with lat-long points in the 10 and 11th columns, respectively
  #' @threshold quantile value for the amount of data points to remove, starting with those farthest from the centroid first. 
  #' 
  #' @output a spatial polygon of MCH
  #' 
  # Function to calculate the convex hull area 
    cut_alphahull_env <- function(input, threshold=.95) {
    input <- input[is.finite(input$fpcat.PC1)==T,]
    hull.pol<-chull(input[,10:11]) #Calculate alpha-hull for full data
    coords <- input[c(hull.pol, hull.pol[1]), ]
    PolOur <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords[,10:11])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
    center <- sf::st_centroid(PolOur)
    ptsob <- sf::st_as_sf(input[,10:11], coords=c("fpcat.PC1", "fpcat.PC2"), crs=crs(PolOur))
    dists <- sf::st_distance(center, ptsob)
    under <- dists<quantile(dists,threshold)
    cut <- input[under[1,]==TRUE,]
    
    hull.pol2<-chull(cut[,10:11]) #Calculate alpha-hull for the cut down data
    coords2 <- cut[c(hull.pol2, hull.pol2[1]), ]
    PolOur2 <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords2[,10:11])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs')))
    return(PolOur2)
  }
}


## ---------------------------------------------------------------------------------------------------------------------------------------
tictoc::tic()
library(parallel)


packed_Clima <- terra::wrap(Clima)
packed_wadminrast <- terra::wrap(wadminrast)


cl <- makeCluster(4, outfile="")

clusterExport(cl, c("dat", "packed_Clima", "packed_wadminrast", "ccs", "env", "makefuncs", "cut_alphahull", "cut_alphahull_env"))

clusterCall(cl, makefuncs)

# Run the simulations in parallel

list_results <- clusterApplyLB(cl, dat$flnm, function(x) {
#list_results <- lapply(dat$flnm, function(x) {
  Clima <- unwrap(packed_Clima)
  wadminrast <- unwrap(packed_wadminrast)
  filtemp <- paste("data/GBIF/occs/",x , sep="") #set up temp filename
  tempdat <- read.csv(filtemp)
  tempdat <- na.omit(tempdat)
  tempdat$cells <- terra::cellFromXY(Clima, xy=tempdat[,5:4]) #Grab the cell ids
  tempdat$envCell <- terra::cellFromXY(env, xy=as.matrix(tempdat[,5:4])) #Grab cell ids from env data

  indcells <- tempdat %>% dplyr::group_by(., cells, envCell) %>% dplyr::slice_sample(n = 1) #remove redundant cells (randomly choose 1 point)
  indcells$ISO3 <- wadminrast[indcells$cells]$GID_0 #grab country ID
  indcells <- ccs %>% dplyr::select(., ISO3, continent) %>% left_join(indcells, ., by=c("ISO3")) #Grab continent
  Envals <- terra::extract(env, as.matrix(indcells[,5:4]))
  indcells <- cbind(indcells, Envals)

  summaryout <- NULL
    for(cont in unique(ccs$continent)){
      temp <- dplyr::filter(indcells, continent==cont)
      if(nrow(temp)==0){ #If continent is absent, skip
        next
    }
      if(nrow(temp)<10){ #Assign NAs if less than 10 occurance points in a continent
      hull_area <- NA
      envArea <- NA
      }
    if(nrow(temp)>9){
      hull <- try(cut_alphahull(temp, threshold=.95)) #Threshold area
      tempMask <-  mask(Clima, hull)[[1]] #Mask around oceans
      #poly <- sf::st_polygon(cbind(list(cbind(c(xvals, xvals[1]), c(yvals, yvals[1]))))) #Create a sf polygon from the mch coords
      #polyvec <- terra::vect(poly) #convert to terra vect iobject for terra::mask()
      #crs(polyvec) <- crs(Clima) #assign the right crs
      
      pol <- terra::as.polygons(tempMask) %>% aggregate()
      crs(pol) <- crs(hull)
      hull_area <- st_area(st_as_sf(pol))
      units(hull_area)$numerator <- c("km", "km") # Get Niche Area
    }
    summary_temp <- data.frame(continent=cont, spArea=hull_area)
    summaryout <- rbind(summaryout, summary_temp)
  }
  envhull <- try(cut_alphahull_env(indcells, threshold = .95))
  envArea <- NA #set a possible default
  if(class(envhull)[1]!="try-error"){
    envArea <- sf::st_area(envhull) }    
  summaryout$envArea <- envArea #if you have a good hull, calculate its area
  return(summaryout)
})

stopCluster(cl=cl)
time <- tictoc::toc()

time <- time$toc-time$tic

save(list_results, file="data/ranges/Rangesizes_Real.Rda")
save(time, file="timing.Rda")

