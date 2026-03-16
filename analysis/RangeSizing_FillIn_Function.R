library(igraph)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(rgbif)
library(curl)
library(BIEN)
library(maps)
library(sf)
library(sp)
sf_use_s2(FALSE)
library(terra)
library(maps)
library(data.table)
library(geodata)
library(geosphere)


#' @points A dataframe with lat-long points in the 5th and 4th column, respectively
#' @threshold quantile value for the amount of data points to remove, starting with those farthest from the centroid first. 
#' @env
#' 
#' @output a spatial polygon object of MCH
#' 
# Function to calculate the convex hull area
cut_mch <- function(points, env, threshold=.95) {
  hull.pol <-chull(points[,5:4]) #Calculate alpha-hull for full data
  coords <- points[c(hull.pol, hull.pol[1]), ] #grab the data rows that make the edges. 
  PolOur <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords[,5:4])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'))) #Convert as a sf object
  center <- sf::st_centroid(PolOur) #Grab our centroid
  ptsob <- sf::st_as_sf(points[,5:4], coords=c("decimallongitude", "decimallatitude"), crs=crs(PolOur)) #make sf point object
  dists <- sf::st_distance(center, ptsob) #Calculate distance from the centroid
  under <- dists < quantile(dists,threshold) #grab the points below a certain distance from the centroid 
  cut <- points[under[1,]==TRUE,]#Subset
  
  hull.pol2<-chull(cut[,5:4]) #Calculate mc-hull for our newly trimmed dataset
  coords2 <- cut[c(hull.pol2, hull.pol2[1]), ] #grab the points definind our newly defined mch
  cutPolygon <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(coords2[,5:4])), ID=1)), proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'))) #Create a spatial polygon object of our cut mch
  
  #Now doing environmental space
  cut$envCell <- terra::cellFromXY(env, xy=as.matrix(cut[,5:4])) #extract env values based on x-y space
  indcells_env <- cut %>% dplyr::group_by(., envCell) %>% dplyr::slice_sample(n = 1) #remove redundant cells (randomly choose 1 point)
  envvals <- terra::extract(env, as.matrix(indcells_env[,5:4])) #Extract environmental coordinates from our PCA
  indcells_env <- cbind(indcells_env, envvals) #add these back into ourdf
  indcells_env <- dplyr::filter(indcells_env, is.na(fpcat.PC1)==F) #Remove NA values
  
  env_hullID <-chull(indcells_env[,c("fpcat.PC1", "fpcat.PC2")]) #Create environmental alpha-hull
  env_coords <- indcells_env[c(env_hullID, env_hullID[1]), ] #Grab the cells defining the environmental edge

  output <- list(geob=cutPolygon, envdf=env_coords)
  return(output)
}

#### Load in our necessary data
Clima <- geodata::worldclim_global("worldclim", var="bio", res=5) #Load in global climate data (just as a mask)
wadmin <- geodata::world(resolution = 5, level =0, path="data/climate") #Download world administrative boundaries
ccs <- geodata::country_codes() #get country codes-includes continent key
wadminrast <- terra::rasterize(wadmin, Clima, field="GID_0") #rasterize so it's faster to join
load(file="data/climate/env.RData") #load in our climate PCA


flnames <- list.files("data/GBIF/occs")
matches <- read.csv(file="backbone_matches.csv")
envsummaryout <- NULL
geosummaryout <- NULL
queryprob <- NULL

envexist <- read.csv(file="data/envsummaryout.csv")
geoexist <- read.csv(file="data/geosummaryout.csv")



for(i in 1:length(flnames)){
  if(flnames[i] %in% c(envexist$file, geoexist$file)){
    next #skip this species if we already have estimates for it
  }
  print(paste("Attempting", i, "of ~7k", sep=" ")) 
  tempdat <- read.csv(file=paste("data/GBIF/occs/",flnames[i], sep=""))
  if(nrow(tempdat)<1){
    print("no occurences")
next
  }
  if(is.null(tempdat$taxonkey)==TRUE){ #Skip species missing the taxonkey column (probably queried with older code)
    queryprob <- c(queryprob, i)
    next
  }
  
  tempdat$cells <- terra::cellFromXY(wadminrast, xy=tempdat[,5:4]) #Grab the cell ids
  indcells <- tempdat %>% dplyr::group_by(., cells) %>% dplyr::slice_sample(n = 1) #remove redundant cells (randomly choose 1 point)
  indcells$ISO3 <- wadminrast[indcells$cells]$GID_0 #grab country ID
  indcells <- ccs %>% dplyr::select(., ISO3, continent) %>% left_join(indcells, ., by=c("ISO3")) #Grab continent
  
  nameinfo <- dplyr::filter(matches, speciesKey==as.numeric(names(sort(-table(tempdat$taxonkey)))[1])) #looks hacky, but just grabing most common element
  if(nrow(nameinfo)>1){ #If multiple accepted names, choose the exact match
    nameinfo <- dplyr::filter(nameinfo, status=="ACCEPTED") #remove synonymns if needed
    nameinfo <- dplyr::filter(nameinfo, matchType=="EXACT")
    if("Bombus terrestris/lucorum" %in% nameinfo$verbatim_name){ #Deal with one problematic split
      nameinfo <- dplyr::filter(nameinfo, verbatim_name=="Bombus terrestris")
    }
    nameinfo <- nameinfo[1,] #If all else fails, just take the first element
  }
  
  
  if(nrow(nameinfo)==0){
    print("no name match")
    next
  }
  
  #### Run geo on each continent
  geosummary <- NULL
  envcords <- NULL
  for(cont in unique(ccs$continent)){
    temp <- dplyr::filter(indcells, continent==cont)
   colnames(temp) <- tolower(colnames(temp))
    if(nrow(temp)==0){ #If continent is absent, skip
      next
    }
    if(nrow(temp)<10){ #Assign NAs if less than 10 occurance points in a continent
      hull_area <- NA
      envArea <- NA
    }
    if(nrow(temp)>9){
      hull <- cut_mch(temp, env=env, threshold=.95) #Threshold area
      tempMask <-  mask(Clima[[1]], hull[[1]]) #Mask around oceans
      contgeo <- terra::expanse(tempMask, unit="km") #Grab extent size
      tempgeo <- data.frame(continent=cont, contgeokm=contgeo) #make a temp dataframe with continent name and range size
      
      print(cont)
      geosummary <- rbind(geosummary, tempgeo) #add to geo summary output for the species
      envcords <- rbind(envcords, hull[[2]]) #Add our environmental range edges for the continent to the total output
    }
  }
  if(is.null(geosummary)==T){ #If we don't have any hulls, skip this species
     next
}
  geosummary$species <- temp$species[1] #Add species ID
  geosummary <- cbind(geosummary, nameinfo)
  geosummary$file <- flnames[i]
  
  ### Run full environmental area
  envcords <- dplyr::filter(envcords, is.na(fpcat.PC1)==F) #If we're missing envirionmental data for any of our points, remove. Should be an exceptionally small number of points. 
  env_hullID <-chull(envcords[,c("fpcat.PC1", "fpcat.PC2")]) #Make mch from the boundaries of our continent-edges
  env_coords <- envcords[c(env_hullID, env_hullID[1]), ] #Grab the coordinants
  cutPolygon <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(env_coords[,c("fpcat.PC1", "fpcat.PC2")])), ID=1)))) #Create a spatial polygon object of our cut mch
  envarea <- sf::st_area(cutPolygon) #Check its area
  
  envsummary <- data.frame(sp=temp$species[1], envarea=envarea, file=flnames[i]) #Save the summary output of the environment
  envsummary <- cbind(envsummary, nameinfo) #add in our naming
  
  envsummaryout <- rbind(envsummaryout, envsummary) #add to our existing stack
  write.csv(envsummaryout, file="data/envsummaryout_fill.csv")
  geosummaryout <- rbind(geosummaryout, geosummary)
  write.csv(geosummaryout, file="data/geosummaryout_fill.csv")
}


#species with missing taxonkeys (probably occurences queried using old code)
#save(queryprob, file="queryprobIDs.RDA")

#Species with issues
#for(q in queryprob[-1]){
# file.rename(from=paste("data/GBIF/occs/",flnames[q], sep=""), to=gsub(".csv", replacement="_old.csv", paste("data/GBIF/occs/",flnames[q], sep="")))
#}
#gsub(".csv", replacement="_old.csv", flnames[q])
#flnames[q]
#paste("old_",flnames[q], sep="")
#paste("data/GBIF/occs/",flnames[q], sep="")
