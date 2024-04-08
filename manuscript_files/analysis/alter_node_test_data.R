## Test removing varying amounts of nodes

## Housekeeping ######
library(tidyverse)

## Set percentages of nodes to remove
node_percentages = seq(0, 0.75, by = 0.05)

## Set replicates
reps = 50

## NA value
NA_value = -116

## Process detection data for KNN ######

## Read in detections
dets <- readr::read_csv("./manuscript_files/data/fingerprint_detections.csv") %>% 
  rename(node = receiver)

# Process for KNN - tags separately
dets_w_tag <- dets %>%
  dplyr::group_by(cp, node, tag) %>%
  dplyr::summarise(mean_RSSI = mean(rssi),
                   .groups = "keep") %>%
  tidyr::pivot_wider(names_from = "node",
                     values_from = "mean_RSSI",
                     names_prefix = "node_") 
## Reformat calibration point locations
cal_pt_location_df <- readr::read_csv("./manuscript_files/data/fingerprint_coordinates.csv")

## Nodes
nodes <- unique(dets$node)

## Replace NAs
dets_w_tag <- dets_w_tag %>% 
  dplyr::left_join(cal_pt_location_df,
                   by = "cp") %>% 
  dplyr::mutate_at(vars(contains("node_")), ~replace(., is.na(.), NA_value)) %>% 
  dplyr::ungroup()

## Remove varying amounts of test data from the nodes ######

## Empty data.frame for all outputs
output_df <- data.frame()

## Set progress bar
pb <- txtProgressBar(min = 0, max = length(node_percentages), style = 3)

## Drop varying number of grid points from 0 - 75%
for(gp_perc in node_percentages){
  
  # gp_perc <- 0.05
  
  ## Progress bar
  Sys.sleep(0.1)
  setTxtProgressBar(pb, which(node_percentages == gp_perc))
  
  ## Repeat 
  for(i in 1:reps){
    
    # i= 1
    
    ## Nodes to keep for testing
    if(gp_perc > 0){
      nodes_to_drop <- paste0("node_", sample(nodes, round(length(nodes)*(gp_perc))))  
    } else{
      nodes_to_drop <- "none"
    }
    
    knn_est <- data.frame()

    ## Test each point
    for (cp_f in unique(dets_w_tag$cp)){
      
      # ## Focal cp
      # cp_f <- unique(dets_w_tag$cp)[10]
      
      ## Calibration points to use for training
      x_tr <- dets_w_tag %>% 
        dplyr::filter(cp != cp_f) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(contains("node"))
      
      y_tr <- dets_w_tag %>% 
        dplyr::filter(cp != cp_f) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(x,
                      y)
      
      ## Tags at that calibration point
      tag_fs <- dets_w_tag[dets_w_tag$cp == cp_f,]$tag
      
      for(tag_f in tag_fs){
        
        ## Calibration point to use for testing
        x_ts <- dets_w_tag %>% 
          dplyr::filter(cp == cp_f & tag == tag_f) %>% 
          
          ## Change nodes to NA
          mutate(across(matches(nodes_to_drop), ~ifelse(is.numeric(.), NA_value, .))) %>%
          dplyr::select(names(x_tr))
        
        y_ts <- dets_w_tag %>% 
          dplyr::filter(cp == cp_f & tag == tag_f) %>% 
          dplyr::ungroup() %>% 
          dplyr::select(x,
                        y)
        
        ## Define model
        model <- ipft::ipfKnn(train_fgp = x_tr, 
                              train_pos = y_tr)
        
        ## Apply KNN
        knnEstimation <- ipft::ipfEstimate(ipfmodel = model, 
                                           test_fgp = x_ts, 
                                           test_pos = y_ts)
        
        ## Get error for calibraiton point
        est_cp_f_df <- data.frame(gp_perc = (1-gp_perc),
                                  rep = i,
                                  tag = tag_f,
                                  cp = cp_f,
                                  error = knnEstimation$errors,
                                  x_cal = y_ts$x,
                                  y_cal = y_ts$y,
                                  x_est = knnEstimation$location$x,
                                  y_est = knnEstimation$location$y)
        
        ## Bind rows
        knn_est <- bind_rows(knn_est,
                             est_cp_f_df)
        
      }
      
    }
    
    output_df <- bind_rows(output_df,
                           knn_est)
    
  }
}

## Check error
output_df_sum <- output_df %>% 
  group_by(gp_perc,rep) %>% 
  summarise(error = median(error),
            pts = n())

## Save outputs
readr::write_csv(output_df,
                 "./manuscript_files/results/altered_node_test_data.csv")
