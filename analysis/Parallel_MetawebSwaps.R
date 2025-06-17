library(parallel)
library(igraph)
library(dplyr)

#Most of the parallel version of this code was written by ChatGPT using the following prompt:Hey ChatGPT, use the parallel package in R to parallelize this code across an arbitrary number of cores. 

load(file="data/Parallelwebs.RDA")


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
nsim <- 1000
bmetaCent_list <- parLapply(cl, 1:nsim, bird_sim_fun)
tictoc::tic()
imetaCent_list <- parLapply(cl, 1:nsim, insect_sim_fun)
tictoc::toc()


# Combine and post-process
bmetaCent_null <- do.call(rbind, bmetaCent_list)
imetaCent_null <- do.call(rbind, imetaCent_list)



# Save results
save(bmetaCent_null, file = "BmetawebSwaps_raw.Rda")
save(imetaCent_null, file = "ImetawebSwaps_raw.Rda")

# Stop the cluster
stopCluster(cl)

###########################################Testing/messing around

#bmetaCent <- left_join(bmetaCent, dist_bmetaCent_null, by="node")
#bmetaCent$close.z.up <- (bmetaCent$close-bmetaCent$meanclose.y)/bmetaCent$sdclose.y


imetaCent <- data.frame(between=igraph::centralization.betweenness(imeta)$res,
                        close=igraph::centralization.closeness(imeta)$res,
                        deg=igraph::centralization.degree(imeta, mode="all")$res,
                        eig=igraph::centralization.evcent(imeta)$vector,
                        node=V(imeta)$name)


library(dplyr)

# Make sure your input data frames look like:
# emp_df: node | between
# sim_df: node | between | indexID

empirical_pvals <- imetaCent %>% 
  rename(obs_close = close) %>%
  left_join(dplyr::select(imetaCent_null, node, close), by = "node") %>%
  group_by(node) %>%
  summarize(
    obs = first(obs_close),
    p_lessCentral = (sum(close >= obs) + 1) / (n() + 1))

sums <- dplyr::group_by(imetaCent_null, node) %>%
  dplyr::summarise(., meancl=mean(close), sdcl=sd(close))

bd <- left_join(imetaCent,sums, by="node")

bd$close.z <- (bd$close-bd$meancl)/bd$sdcl


hist(bd$close.z, breaks=500)


empirival_pval <- imetaCent %>%
  left_join(imetaCent_null, ., by="node")



imetaCent2 <- left_join(imetaCent, dist_imetaCent_null_2, by="node")
imetaCent2$close.z.up <- (imetaCent2$close-imetaCent2$meanclose)/imetaCent2$sdclose


plot(x=imetaCent2$close.z.up, y=imetaCent2_100$close.z.up)
abline(a=0, b=1, col="red")

plot(x=imetaCent2$sdclose, y=imetaCent2_100$sd)
abline(a=0, b=1, col="red")

plot(x=imetaCent2$meanclose, y=imetaCent2_100$meanclose)
abline(a=0, b=1, col="red")

summary(lm(imetaCent2_100$sd~imetaCent2$sdclose))


summary(lm(imetaCent2_100$close.z.up~imetaCent2$close.z.up))

summary(lm(imetaCent2_100$meanclose~imetaCent2$meanclose))



imetaCent2_100 <- imetaCent2

bmetaCent %>% 
  ggplot(., aes(x=scale(log(envarea)), y=close.z, col=net))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm")+
  facet_wrap(~type)+
  theme_classic()+
  xlab("Scaled Log Environmental Area")+
  ylab("Closeness Centrality (z-transformed)")+guides(col="none")+
  scale_color_manual(labels = c("Avian-Plant", "Insect-Plant"), values = c("#0F882C", "dodgerblue"))

for (j in unique(Insectweb$site_index)) {
  temp <- Insectweb[Insectweb$site_index == j, ]
  print(paste("site ",j, " has ", nrow(temp), " rows", sep=""))
}

