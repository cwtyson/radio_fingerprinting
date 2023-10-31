## This script takes the detections processed in the script "1_process_raw_detections.csv" and creates the radio map from these detections.
## Here the locations of the fingerprints within the radio map are estimated using the remaining fingerprints
## The locations

## Housekeeping ######
library(tidyverse)
library(ipft)

## Process detection data for radio fingerprinting ######

## Read in detections
dets <- readr::read_csv("./radio_fingerprinting_workflow/fingerprinting_files/fingerprint_detections.csv")

# Process 
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

## A. Predict locations of test points ########

## Empty data frame for results
knn_est <- data.frame()

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


## B. Predict locations of new radio fingerprints ########

## Example of predicting new radio fingerprints 

## Read in detections
dets <- readr::read_csv("./radio_fingerprinting_workflow/fingerprinting_files/fingerprint_detections.csv")

# Process 
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

## New radio fingerprints
new_dets <- readr::read_csv("./radio_fingerprinting_workflow/fingerprinting_files/fingerprint_detections.csv")

# Process. Receiver/nodes columns must be in the same order as in the radio map dataset and all the same columns must be present
new_dets <- new_dets %>%
  dplyr::group_by(cp, receiver, tag) %>%
  dplyr::summarise(mean_RSSI = mean(rssi),
                   .groups = "keep") %>%
  tidyr::pivot_wider(names_from = "receiver",
                     values_from = "mean_RSSI",
                     names_prefix = "node_") %>% 
  dplyr::mutate_at(vars(contains("node_")), ~replace(., is.na(.), -115)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-cp,-tag) %>%
  
  ## Take random subset
  slice_sample(prop = 0.2)


## Create radio map of training data from full data: 

## Fingerprints
x_tr <- dets %>% 
  dplyr::ungroup() %>% 
  dplyr::select(contains("node"))

## Coordinates
y_tr <- dets %>% 
  dplyr::ungroup() %>% 
  dplyr::select(x,
                y)

## Define model based on radio map data
model <- ipft::ipfKnn(train_fgp = x_tr, 
                      train_pos = y_tr)


## Apply estimation
estimation <- ipft::ipfEstimate(ipfmodel = model, 
                                test_fgp = new_dets)

## Get estimated locations
estimated_locs <- data.frame(estimation$location)
  

