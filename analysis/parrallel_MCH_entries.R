library(igraph)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(rgbif)
library(curl)
library(BIEN)
library(maps)
library(sf)
sf_use_s2(FALSE)
library(terra)
library(maps)
library(data.table)
library(geodata)
library(parallel)

load("data/MetaWebCentrality.Rda") #Load in Global Metaweb


#'
#' @param x latitude
#' @param y longitude
#' @param msk environmental mask to use, so that we don't count ocean cells
#' @param method method for calculating hulls (case insentitive). Current acceptable values are "mch" and "alphahull"
#' @param alphaval Alpha value used for calculating alpha hull. Defaults to 20.
#'
#' @returns geographic range size in terms of # of cells covered in environmental mask

getArea <- function(x,y, equalArea=TRUE, msk, method, alphaval=20, returnHull=TRUE){
  if(tolower(method) %in% c("mch", "alphahull")==FALSE){
    stop("No usable method provided")
  }

  rem <- which(is.na(x) | is.na(y))
  if(any(rem)){
    x <- x[-rem]
    y <- y[-rem]
  }
  if(tolower(method)=="mch"){ #Use grDevices to create convex hull
    tmp <- grDevices::chull(x,y)
    x <- x[tmp]
    y <- y[tmp]
    poly <- sf::st_polygon(cbind(list(cbind(c(x, x[1]), c(y, y[1]))))) #Create a sf polygon from the mch coords
    polyvec <- terra::vect(poly) #convert to terra vect iobject for terra::mask()
  }
  if(tolower(method)=="alphahull"){ #Use alphaHull to creat alpha hulls
    tmp <-ashape(x = x, y=y, alpha = alphaval)
    tmp <- data.frame(tmp$edges)[,c( 'x1', 'y1', 'x2', 'y2')]
    l <- st_linestring(matrix(as.numeric(tmp[1,]), ncol=2, byrow = T))
    for(i in 2:nrow(tmp)){
      l <- c(l, st_linestring(matrix(as.numeric(tmp[i,]), ncol=2, byrow = T)))
    }
    polyvec <- vect(st_sf(geom = st_sfc(l), crs = crs(msk)) %>% st_polygonize() %>% st_collection_extract())
  }
  crs(polyvec) <- crs(msk)
  tempMask <-  mask(msk, polyvec) #mask our temperature data by the polygon (this maintains 0 for oceans and whatnot)
  df_dat <- terra::as.data.frame(tempMask[[1]], na.rm=T) #Extract the non NA values within the mask
  cellcount <- nrow(df_dat)

  if(returnHull==FALSE){
    return(cellcount)
  }
  if(returnHull==TRUE){
    ret <- list(count=cellcount, pgon=polyvec)
    return(ret)
  }
}




#Example code from stuff.
tictoc::tic()
cl <- makeCluster(detectCores()-2, outfile="")

clusterExport(cl, c("nsites", "runSimulation", "runDispersalSim","simModularity"))

# Run the simulations in parallel
list_results <- clusterApplyLB(cl, 1:num_iterations, function(i) {
  result <- runSimulation()
  name <- paste("../Data/Mem_usage/Neutime", i, ".txt", sep="")
  system(paste('free -g | cat > ', name, sep = ""))
  return(result)
})


stopCluster(cl)
time3 <- tictoc::toc()
