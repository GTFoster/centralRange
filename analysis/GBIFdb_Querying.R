library(tidyverse)
library(gbifdb)
library(rgbif) #Only used to find taxonkey
library(CoordinateCleaner)

load(file="toQuery_centralRange.Rda") #Load in our total species list
gbif <- gbif_local(dir='/home/shared/occurrence/2024-04-01') #Tell R where to find our local GBIF copy
#I apologize to the reproducibility gods for absolute paths

nmlist <- data.frame(Name=NA, ID=1, NROW=0) #Start an empty df to keep track of query progress

for(name in unique(unlist(toQuery))[886:length(toQuery)]){ #For each species, starting at last break
spdat <- gbif |> 
  filter(species == name) |> 
  collect() #Grab all occurence records for the species of interest and call is spdat

filname <- name %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
  paste("data/GBIF/occs/", ., ".csv", sep="") #Create a filename

spdat <- spdat %>% #Do some data cleaning
    dplyr::filter(occurrencestatus  == "PRESENT") %>%
    dplyr::filter(!basisofrecord %in% c("FOSSIL_SPECIMEN")) %>% #Remove fossils
    CoordinateCleaner::cc_cen( #Remove records within 1k buffer of country centroids
    lon = "decimallongitude", 
    lat = "decimallatitude", 
    buffer = 1000, # radius of circle around centroid to look for centroids
    value = "clean",
    test="both")  %>% 
    cc_sea( #Remove oceanic records
  lon = "decimallongitude",
  lat = "decimallatitude"
  )# %>%
  #dplyr::filter(., establishmentMeans != "introduced") #Remove known introduced  spp

occs <- dplyr::select(spdat, genus, species,scientificname, decimallatitude, decimallongitude, elevation, day, month, year, taxonkey) %>% unique()
write.csv(occs, file=filname, row.names = FALSE)

tempcheck <- data.frame(Name=name, ID=max(nmlist$ID)+1, NROW=try(nrow(occs)))
nmlist <- rbind(nmlist, tempcheck)

write.csv(nmlist, file="Progresscheck.csv")
}
