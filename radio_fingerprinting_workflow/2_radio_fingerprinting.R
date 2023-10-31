## This script takes the detections processed in the script "1_process_raw_detections.csv" and creates the radio map from these detections
## 

## Housekeeping ######
library(tidyverse)

## Process detection data for radio fingerprinting ######

## Read in detections
dets <- readr::read_csv("./radio_fingerprinting_workflow/fingerprinting_files/fingerprint_detections.csv")

# Process for KNN 
dets <- dets %>%
  dplyr::group_by(cp, receiver, tag) %>%
  dplyr::summarise(mean_RSSI = mean(rssi),
                   .groups = "keep") %>%
  tidyr::pivot_wider(names_from = "receiver",
                     values_from = "mean_RSSI",
                     names_prefix = "node_") 


## Read in coordinates of fingerprints
pt_locations <- readr::read_csv("./radio_fingerprinting_workflow/fingerprinting_files/fingerprint_coordinates.csv")

## Add coordinates 
dets <- dets %>% 
  dplyr::left_join(pt_locations,
                   by = "cp") %>% 
  dplyr::mutate_at(vars(contains("node_")), ~replace(., is.na(.), -115))

## Empty data frame for results
knn_est <- data.frame()

## Predict locations
for (cp_f in unique(dets$cp)){
  
  ## Calibration points to use for training
  x_tr <- dets %>% 
    dplyr::filter(cp != cp_f) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(contains("node"))
  
  y_tr <- dets %>% 
    dplyr::filter(cp != cp_f) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(x,
                  y)
  
  ## Tags at that calibration point
  tag_fs <- dets[dets$cp == cp_f,]$tag
  
  for(tag_f in tag_fs){
    
    ## Point to use for testing
    x_ts <- dets %>% 
      dplyr::filter(cp == cp_f  & tag == tag_f) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(contains("node"))
    
    y_ts <- dets %>% 
      dplyr::filter(cp == cp_f & tag == tag_f) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(x,
                    y)
    
    ## Define model
    model <- ipft::ipfKnn(train_fgp = x_tr, 
                          train_pos = y_tr)
    
    ## Apply estimation
    knnEstimation <- ipft::ipfEstimate(ipfmodel = model, 
                                       test_fgp = x_ts, 
                                       test_pos = y_ts)
    
    ## Get error for calibration point
    est_cp_f_df <- data.frame(tag = tag_f,
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

## Save outputs
readr::write_csv(knn_est,
                 "/radio_fingerprinting_workflow/location_estimates.csv")
