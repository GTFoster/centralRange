library(parallel)
library(igraph)
library(dplyr)

#Most of the parallel version of this code was written by ChatGPT using the following prompt:Hey ChatGPT, use the parallel package in R to parallelize this code across an arbitrary number of cores. 

load(file="Parallelwebs.RDA")


#Creating the site index
Birdweb <- Birdweb %>%
  group_by(latitude, longitude) %>%
  mutate(site_index = cur_group_id()) %>%
  ungroup()

Insectweb <- Insectweb %>%
  group_by(latitude, longitude) %>%
  mutate(site_index = cur_group_id()) %>%
  ungroup()

# Define number of cores and simulations
ncores <- detectCores() - 2 #leave an extra core to appease the daemons
cl <- makeCluster(ncores)

clusterEvalQ(cl, { #Make sure workers have necessary libraries
  library(igraph)
  library(dplyr)
  library(vegan)
})

# Export necessary data and functions to the cluster
clusterExport(cl, varlist = c("Birdweb", "Insectweb"))

# Function to process Birdweb simulations
bird_sim_fun <- function(i) {
  Vset <- NULL
  Eset <- NULL
  #lay <- layout_as_bipartite(rg)
  for (j in unique(Birdweb$site_index)) {
    temp <- Birdweb[Birdweb$site_index == j, ]
    if (nrow(temp) < 3) next
    temp <- dplyr::select(temp, species1, species2) %>% unique()
    nullsetup <- vegan::nullmodel(table(temp), method = "quasiswap")
    out <- simulate(nullsetup, nsim = 1)
    if(j!=20) {rg <- igraph::graph_from_incidence_matrix(out[, , 1])}
    if(j==20) {rg <- igraph::graph_from_incidence_matrix(table(temp$species1, temp$species2))}
    
    #plot(rg, layout=lay)
    Vtemp <- igraph::as_data_frame(rg, "vertices") %>% unique()
    Vset <- rbind(Vset, Vtemp) #Vertex list across all sites
    Etemp <- igraph::as_data_frame(rg)
    Eset <- rbind(Eset, Etemp) #Edge list across all sites
    print(j)
  }
  if (is.null(Eset) || is.null(Vset)) return(NULL)
  Vset <- unique(Vset)
  Eset <- unique(Eset)
  metatemp <- graph_from_data_frame(Eset, directed = FALSE, vertices = Vset$name)
  
  bmetaCent_temp <- data.frame(
    between = igraph::centralization.betweenness(metatemp)$res,
    close   = igraph::centralization.closeness(metatemp)$res,
    deg     = igraph::centralization.degree(metatemp, mode = "all")$res,
    eig     = igraph::centralization.evcent(metatemp)$vector,
    node    = V(metatemp)$name,
    indexID = i
  )
  return(bmetaCent_temp)
}

# # Function to process Insectweb simulations
insect_sim_fun <- function(i) {
  Vset <- NULL
  Eset <- NULL

  for (j in unique(Insectweb$site_index)) {
    print(paste("starting", j))
    temp <- Insectweb[Insectweb$site_index == j, ]
    if (nrow(temp) < 3) next
    temp <- dplyr::select(temp, species1, species2) %>% unique()
    nullsetup <- vegan::nullmodel(table(temp), method = "greedyqswap")
    if (nullsetup$ncol == 1 || nullsetup$nrow == 1) next
    out <- simulate(nullsetup, nsim = 1, thin=1)
    if(j!="space") {rg <- igraph::graph_from_incidence_matrix(out[, , 1])}
    #if(j==47) {rg <- igraph::graph_from_incidence_matrix(table(temp$species1, temp$species2))}
    #rg <- igraph::graph_from_incidence_matrix(out[, , 1])
    Vtemp <- igraph::as_data_frame(rg, "vertices") %>% unique()
    Vset <- rbind(Vset, Vtemp)
    Etemp <- igraph::as_data_frame(rg)
    Eset <- rbind(Eset, Etemp)
  }
  if (is.null(Eset) || is.null(Vset)) return(NULL)
  Vset <- unique(Vset)
  Eset <- unique(Eset)
  metatemp <- graph_from_data_frame(Eset, directed = FALSE, vertices = Vset$name)

  imetaCent_temp <- data.frame(
    between = igraph::centralization.betweenness(metatemp)$res,
    close   = igraph::centralization.closeness(metatemp)$res,
    deg     = igraph::centralization.degree(metatemp, mode = "all")$res,
    eig     = igraph::centralization.evcent(metatemp)$vector,
    node    = V(metatemp)$name,
    indexID = i
  )
  return(imetaCent_temp)
}

# Run simulations in parallel
nsim <- 2000
bmetaCent_list <- parLapply(cl, 1:nsim, bird_sim_fun)
tictoc::tic()
imetaCent_list <- parLapply(cl, 1:nsim, insect_sim_fun)
tictoc::toc()


# Combine and post-process
bmetaCent_null <- do.call(rbind, bmetaCent_list)
imetaCent_null <- do.call(rbind, imetaCent_list)



# Save results
save(bmetaCent_null, file = "BmetawebSwaps_raw2.Rda")
save(imetaCent_null, file = "ImetawebSwaps_raw2.Rda")

# Stop the cluster
stopCluster(cl)
