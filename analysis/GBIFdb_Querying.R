library(tidyverse)
library(gbifdb)
library(rgbif) #Only used to find taxonkey
library(CoordinateCleaner)


load(file="toQuery_centralRange.Rda") #Load in our total species list
gbif <- gbif_local(dir='/home/shared/occurrence/2024-04-01') #Tell R where to find our local GBIF copy
#I apologize to the reproducibility gods for absolute paths

matches <- read.csv(file="backbone_matches.csv")

nmlist <- data.frame(Name=NA, ID=0, NROW=0, time=Sys.time()) #Start an empty df to keep track of query progress

keyquery <- NULL
for(name in matches$verbatim_name){
  filname <- name %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
      paste("data/GBIF/occs/", ., ".csv", sep="") #Create a filename
    if(file.exists(filname)==TRUE){ #If the occ file already exists locally and is filled, skip and move to the next species
      if(file.info(filname)$size > 150){
      next
      }
    }
  key <- matches$speciesKey[matches$verbatim_name==name] #Save the key to search by
  keyquery <- rbind(keyquery, data.frame(name=name, key=key))
}

keyquery$name[keyquery$name=="Bombus terrestris/lucorum"] <- "Bombus terrestris-lucorum"

#Each query seems to take ~8min, pretty much reguardless of size. So to speed things up, we're going to call our gbifdb query in chunks and then split the outputs up afterwards
idex <- c(0, seq(from=1, to=nrow(keyquery), by=30)[-1], nrow(keyquery))

#for(j in idex[-length(idex)]){
for(j in c(241,271,301,331,361,391,421,451,481,511,511)){
  keystemp <- keyquery$key[(j+1):(j+30)]
  Sys.time()
    spdat <- gbif |> 
      #filter(taxonkey %in% key) |> 
      filter(taxonkey %in% keystemp) |> 
      collect() #Grab all occurence records for the species of interest and call is spdat, using our backbone key to search
    Sys.time()
    outrun <- spdat
    spdat <- dplyr::filter(spdat, is.na(decimallongitude)==F) #remove recs without latlong
    
    spdat <- spdat %>% #Do some data cleaning
        dplyr::filter(occurrencestatus  == "PRESENT") %>%
        dplyr::filter(!basisofrecord %in% c("FOSSIL_SPECIMEN")) %>% #Remove fossils
        CoordinateCleaner::cc_cen( #Remove records within 1km buffer of country centroids
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
    #subset across species and save into seperate files
    
    for(i in (j+1):(j+30)){
      filname <- keyquery$name[i] %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
        paste("data/GBIF/occs/", ., ".csv", sep="") #Create a filename
      if(filname=="data/GBIF/occs/bombus_terrestris/lucorum.csv"){
        next
      }
      tempocc <- dplyr::filter(occs, taxonkey==keyquery$key[i])
      write.csv(tempocc, file=filname, row.names = FALSE) 
      IDkey <- matches$verbatim_index[matches$verbatim_name==keyquery$name[i]]
      if(length(IDkey)==0){
        IDkey <- NA
      }
      tempcheck <- data.frame(Name=keyquery$name[i], ID=IDkey, NROW=try(nrow(tempocc)), time=Sys.time())
    }
    
    nmlist <- rbind(nmlist, tempcheck)
    rm(occs, tempocc)
    gc()
    write.csv(nmlist, file="Progresscheck_fill2.csv")
}
