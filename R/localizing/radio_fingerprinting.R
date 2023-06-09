## KNN localization estimation

## Housekeeping ######
library(tidyverse)

## Process detection data for KNN ######

## Read in detections
dets <- readr::read_csv("./data/detections/field_cal_dets.csv")

# Process for KNN 
dets_w_tag <- dets %>%
  dplyr::rename(cp = calibration_point) %>%
  dplyr::group_by(cp, node, tag) %>%
  dplyr::summarise(mean_RSSI = mean(RSSI),
                   .groups = "keep") %>%
  tidyr::pivot_wider(names_from = "node",
                     values_from = "mean_RSSI",
                     names_prefix = "node_") 

## Test point locations
cal_pt_location <- sf::st_read("./data/field/WES_calibration_points.KML") %>%
  sf::st_set_crs(4326) %>% 
  sf::st_transform(28992) 

## Reformat calibration point locations
pt_locations <- cal_pt_location %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::mutate(cp = cal_pt_location$Name) %>%
  dplyr::select(cp,
                x = X,
                y = Y)

## Add coordinates of the test points
dets_w_tag <- dets_w_tag %>% 
  dplyr::left_join(pt_locations,
                   by = "cp") %>% 
  dplyr::mutate_at(vars(contains("node_")), ~replace(., is.na(.), -115))

## Empty data frame for results
knn_est <- data.frame()

## Predict locations
for (cp_f in unique(dets_w_tag$cp)){
  
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
    
    ## Point to use for testing
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
                 "./data/outputs/radio_fingerprint_localization_estimates.csv")
