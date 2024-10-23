library(tidyverse)
library(gbifdb)
library(rgbif) #Only used to find taxonkey

load(file="toQuery_centralRange.Rda") #Load in our total species list
gbif <- gbif_local(dir='../../../GBIF_local/occurrence/2024-04-01') #Tell R where to find our local GBIF copy

nmlist <- data.frame(Name=NA, ID=1) #Start an empty df to keep track of query progress

for(name in unique(unlist(toQuery))){ #For each species
spdat <- gbif |> 
  filter(species == name) |> 
  collect() #Grab all occurence records for the species of interest and call is spdat

filname <- name %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
  paste("data/GBIF/occs/", ., ".csv", sep="") #Create a filename

spdat <- spdat %>% #Do some data cleaning
    dplyr::filter(occurrenceStatus  == "PRESENT") %>%
    dplyr::filter(!basisOfRecord %in% c("FOSSIL_SPECIMEN")) %>% #Remove fossils
    CoordinateCleaner::cc_cen( #Remove records within 1k buffer of country centroids
    lon = "decimalLongitude", 
    lat = "decimalLatitude", 
    buffer = 1000, # radius of circle around centroid to look for centroids
    value = "clean",
    test="both")  %>% 
    cc_sea( #Remove oceanic records
  lon = "decimalLongitude",
  lat = "decimalLatitude"
  )# %>%
  #dplyr::filter(., establishmentMeans != "introduced") #Remove known introduced  spp

occs <- dplyr::select(spdat, genus, species,scientificName, decimalLatitude, decimalLongitude, elevation, day, month, year, taxonkey) %>% unique()
write.csv(occs, file=filname, row.names = FALSE)

tempcheck <- data.frame(Name=name, ID=max(nmlist$ID)+1)
nmlist <- rbind(nmlist, tempcheck)

write.csv(nmlist, file="Progresscheck.csv")
}