# load some packages
library(readr)
library(dplyr)
library(tidyr)
library(sfnetworks)
library(sf)
library(tidygraph)
library(ggplot2)

# Some ideas here:
#    https://r-spatial.org//r/2019/09/26/spatial-networks.html
# But this code is adapted from here:
#    https://luukvdmeer.github.io/sfnetworks/articles/sfn04_routing.html

# load data
streams <- st_read("data/major_clipped_streams.shp")
gp <- st_read("data/snapped_coordinates_full_siblings.shp")
gp_half <- st_read("data/snapped_half_sib_coordinates.shp")
gp_pairs <- read_csv("data/full_sib_pairs.csv")
gp_half_pairs <- read_csv("data/half_sib_pairs.csv")

# plot this to check it out
streams %>%
  ggplot() +
  geom_sf() +
  geom_sf(data = gp)

# turn the stream network into a network object that works
#   with the sfnetworks package
net <- streams %>%
  st_cast("LINESTRING") %>%
  as_sfnetwork(directed = FALSE) %>%
  st_transform(3111) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

# match the CRS of the network
gp <- gp %>%
  st_transform(st_crs(net))
gp_half <- gp_half %>%
  st_transform(st_crs(net))

# drop the GP points into the network and cut at each, so the
#   network has nodes close to each sample
net_gp <- st_network_blend(net, gp)
net_gp_half <- st_network_blend(net, gp_half)

# grab the edge list from the network
net_gp_edge <- net_gp %>%
  activate("edges")
net_gp_half_edge <- net_gp_half %>%
  activate("edges")

# function to calculate the length of a path from a list of edges
path_distance <- function(idx, edges) {
  edges %>%
    slice(idx) %>%
    mutate(edge_length = st_length(geometry)) %>%
    pull(edge_length) %>%
    sum()
}

# use st_network_paths to find the shortest path between each pair of
#   of samples and then calculate summed edge lengths (weights) over
#   each path
gp_dist <- rep(NA, nrow(gp_pairs))
for (i in seq_len(nrow(gp_pairs))) {
  
  # calculate the shortest path from sample i to all other samples
  paths <- net_gp %>%
    st_network_paths(
      from = gp %>% filter(offspring2 == gp_pairs$offspring1_id[i]), 
      to = gp %>% filter(offspring2 == gp_pairs$offspring2_id[i]), 
      weights = "weight"
    )

  # grab the list of edges connecting each pair of points
  edge_list <- paths %>% pull(edge_paths)
  
  # and then calculate distances from this
  gp_dist[i] <- sapply(edge_list, path_distance, edges = net_gp_edge)

}

# repeat for half sibs
gp_half_dist <- rep(NA, nrow(gp_half_pairs))
for (i in seq_len(nrow(gp_half_pairs))) {
  
  # calculate the shortest path from sample i to all other samples
  paths <- net_gp_half %>%
    st_network_paths(
      from = gp_half %>% filter(offspring_ == gp_half_pairs$offspring1_id[i]), 
      to = gp_half %>% filter(offspring_ == gp_half_pairs$offspring2_id[i]), 
      weights = "weight"
    )
  
  # grab the list of edges connecting each pair of points
  edge_list <- paths %>% pull(edge_paths)
  
  # and then calculate distances from this
  gp_half_dist[i] <- sapply(edge_list, path_distance, edges = net_gp_half_edge)
  
}

# add to pairs files
gp_pairs <- gp_pairs %>%
  mutate(distance_m = gp_dist)
gp_half_pairs <- gp_half_pairs %>%
  mutate(distance_m = gp_half_dist)

# and save to file
write_csv(gp_pairs, file = "outputs/gp-full-sib-distance.csv")
write_csv(gp_half_pairs, file = "outputs/gp-half-sib-distance.csv")


# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

#subset of stocked and wild origin pairs
first_order_pairs <- read.csv("./data/fulls_sib_pairs_plus_Rdistances_origins_subset_to_plot_WildStocked.csv", 
                              header = TRUE)


second_order_pairs <- read.csv("./data/half_sib_pairs_plus_Rdistances_origins_subset_to_plot_wild_stocked_pairs.csv",
                               header=TRUE)


# Represent it
  p <- first_order_pairs  %>%
  ggplot( aes(x=(distance_m/1000), fill=pair_origin)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#bdbfbf", "#404080")) +
  theme_ipsum() +
  labs(fill="")
  
  p2 <- second_order_pairs  %>%
    ggplot( aes(x=(distance_m/1000), fill=pair_origin)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#bdbfbf", "#404080")) +
    theme_ipsum() +
    labs(fill="")
  
  
  p3 <- second_order_pairs  %>%
    ggplot( aes(x=(distance_m/1000), fill=pair_origin_one_stocked_counts_as_stocked)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#bdbfbf", "#404080")) +
    theme_ipsum() +
    labs(fill="")
