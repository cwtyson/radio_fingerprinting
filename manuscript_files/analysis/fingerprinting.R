## KNN localization estimation

## Housekeeping ######
library(tidyverse)

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

## Summarise nodes
nodes_sum <- dets_w_tag %>%
  ungroup() %>% 
  # select(contains("node")) %>%
  rowwise() %>%
  mutate(nodes = sum(!is.na(c_across(contains("node")))))

dets_w_tag <- dets_w_tag %>% 
  dplyr::left_join(cal_pt_location_df,
                   by = "cp") %>% 
  dplyr::mutate_at(vars(contains("node_")), ~replace(., is.na(.), -116))

## Empty data frame for results
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
      dplyr::filter(cp == cp_f  & tag == tag_f) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(contains("node"))
    
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
    est_cp_f_df <- data.frame(tag = tag_f,
                              cp = cp_f,
                              error = knnEstimation$errors,
                              x_cal = y_ts$x,
                              y_cal = y_ts$y,
                              x_est = knnEstimation$location$x,
                              y_est = knnEstimation$location$y,
                              node_count = nodes_sum[nodes_sum$cp == cp_f & nodes_sum$tag == tag_f,]$nodes)
    
    ## Bind rows
    knn_est <- bind_rows(knn_est,
                         est_cp_f_df)
    
  }
}

median(knn_est$error)

## Save outputs
readr::write_csv(knn_est,
                 "./manuscript_files/results/fingerprinting_estimates.csv")



