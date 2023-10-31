## Script to estimate calibration locations using multilateration-round 1

## Housekeeping ######
library(tidyverse)
source("/Users/tyson/Documents/git/wes-tracking-test/R/localizing/ipft_functions.R")

## Hard coded values
dist <- 100*1.5 ## Distance filter
mod_alpha = 4.169569
mod_R0 = -44
max_RSSI = -41
min_RSSI = -115

## Process detection data for multilateration ######
dets <- readr::read_csv("./manuscript_files/data/fingerprint_detections.csv") %>% 
  rename(node = receiver)
nodes_df <- readr::read_csv("./manuscript_files/data/receiver_coordinates.csv") %>% 
  rename(node = receiver)


## Summarize
dets_sum <- dets %>% 
  dplyr::group_by(tag,
                  cp,
                  node) %>% 
  dplyr::summarise(mean_RSSI = round(mean(rssi,
                                          na.rm = T), 0),
                   dets_node = n()) %>% 
  dplyr::ungroup()

## Identify node with strongest signal at each calibration point
cp_max_node <- dets_sum %>% 
  
  ## Get max RSSI values within each window
  dplyr::group_by(cp) %>%
  dplyr::mutate(max_RSSI = max(mean_RSSI)) %>%
  dplyr::filter(mean_RSSI == max_RSSI) %>% 
  dplyr::distinct(cp,
                  .keep_all = T)

## Empty data frame
dets_all <- data.frame()

## Set progress bar
pb2 <- txtProgressBar(min = 0, max = length(unique(cp_max_node$cp)), style = 3)

## For each window
for(cp_f in unique(cp_max_node$cp)){
  
  ## Progress bar
  Sys.sleep(0.1)
  setTxtProgressBar(pb2, which(unique(cp_max_node$cp) == cp_f))
  
  # cp_f = unique(cp_max_node$cp)[1]
  # cp_f = "Z_13"
  
  ## Get grid_points with detections
  cp_gps <- dets_sum %>% 
    dplyr::filter(cp == cp_f) %>%
    dplyr::pull(node) %>% 
    unique() 
  
  ## If possible:
  if( length(cp_gps) > 0 ){
    
    ## Filter node pts based on those with detections
    node_pts_df <- nodes_df %>% 
      dplyr::filter(node %in% cp_gps)
    
    ## Get distance between nodes with detections during the point
    n_dist <- raster::pointDistance(node_pts_df[,c("x", "y")], 
                                    node_pts_df[,c("x", "y")], 
                                    lonlat = F,
                                    allpairs = T)
    
    # Make matrix into a dataframe with a row for NodeId
    n_dist_df <- data.frame(n_dist)
    colnames(n_dist_df) <- node_pts_df$node
    n_dist_df$node <- colnames(n_dist_df)
  }
  
  ## Keep nodes within specified distance filter
  nodes_dist_f <- n_dist_df %>%
    dplyr::filter(node == cp_max_node[cp_max_node$cp == cp_f,]$node) %>% 
    tidyr::gather(key = "node", 
                  value = "distance") %>%
    dplyr::filter(distance <= dist) 
  
  ## Filter calibration detections based on nodes within distance filter
  dt_r_dets_f <- dets_sum %>% 
    na.omit() %>% 
    dplyr::filter(node %in% nodes_dist_f$node,
                  cp == cp_f) 
  
  ## Bind to other cp dets
  dets_all <- dplyr::bind_rows(dt_r_dets_f,
                               dets_all)
  
}

## End progress bar
close(pb2)

## Read in calibration locations
cal_locs <- read_csv("./manuscript_files/data/fingerprint_coordinates.csv")

## Join
dets_all <- dets_all %>%
  left_join(cal_locs,
            by = "cp")


## Multilateration  ##########

## Make wide
cal_dets_w <- dets_all %>%
  
  ## Get mean of filtered values
  dplyr::group_by(tag,
                  cp) %>% 
  dplyr::mutate(nodes = n_distinct(node)) %>% 
  
  dplyr::ungroup() %>% 
  dplyr::select(-dets_node) %>% 
  dplyr::group_by(tag) %>% 
  tidyr::pivot_wider(names_from = "node",
                     values_from  = "mean_RSSI",
                     names_prefix =  "node_") %>% 
  ungroup()

## Apply multilateration to RSSI values ######

## Fingerprints and locations
fingerprints <-  cal_dets_w %>% 
  dplyr::filter(nodes >= 3) %>% 
  dplyr::select(contains("node_"))

positions <-  cal_dets_w %>% 
  dplyr::filter(nodes >= 3) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(x, y)

## Filter node pts
node_pts_f <- nodes_df %>% 
  dplyr::filter(paste0("node_",node) %in% names(cal_dets_w)) %>% 
  dplyr::arrange(match(paste0("node_",node), names(fingerprints))) %>% 
  dplyr::select(x,y)

## Define model
proxModel <- ipfProximity(bpos = node_pts_f, 
                          alpha = mod_alpha, 
                          rssirange = c(min_RSSI, 
                                        max_RSSI), 
                          norssi = -125, 
                          wapPow1 = mod_R0)


proxEstimation <- ipfEstimate(ipfmodel = proxModel, 
                              test_fgp = fingerprints, 
                              test_pos = positions)

## Summarise RSSI values
rssi_sum <- cal_dets_w %>% 
  dplyr::filter(nodes >= 3) %>% 
  select(cp, matches("[0-9*]")) %>% 
  mutate(mean_RSSI =  round(rowMeans(across(where(is.numeric)), na.rm = TRUE),0)) %>% 
  rowwise() %>% 
  mutate(max_RSSI = max(across(matches("[0-9*]")), na.rm = TRUE)) %>% 
  select(cp,
         mean_RSSI,
         max_RSSI)


## Get estimated values
ml_est <- data.frame(tag = cal_dets_w %>% 
                       dplyr::filter(nodes >= 3) %>% 
                       dplyr::pull(tag),
                     cp = cal_dets_w %>% 
                       dplyr::filter(nodes >= 3) %>% 
                       dplyr::pull(cp),
                     x = positions$x,
                     y = positions$y,
                     error = proxEstimation$errors[[1]],
                     x_est = proxEstimation$location$x,
                     y_est = proxEstimation$location$y)


## Save
readr::write_csv(ml_est,
                 "./manuscript_files/results/multilateration_estimates.csv")

