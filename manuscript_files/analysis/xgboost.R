## Script to estimate calibration locations using random forests

## Housekeeping ######
library(tidyverse)
#library(ranger)  
library(xgboost)

## Plotting values
color_pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")


## Process detection data for random forests ######

## Read in detections
dets <- readr::read_csv("field_cal_dets.csv")

# Process for random forests
dets_w <- dets %>%
  dplyr::rename(cp = calibration_point) %>%
  dplyr::group_by(cp, tag, node) %>%
  dplyr::summarise(mean_RSSI = mean(RSSI)) %>%
  tidyr::pivot_wider(names_from = "node",
                     values_from = "mean_RSSI",
                     names_prefix = "node_") %>%
  dplyr::filter(!(grepl("WNC", cp))) %>% 
  dplyr::ungroup(cp)

## Read in calibration locations
cal_pts <- sf::read_sf("WES_calibration_points.KML") %>% 
  sf::st_transform(28992)

cal_pt_location_df <- cal_pts %>% 
  sf::st_coordinates() %>% 
  as.data.frame() %>% 
  dplyr::mutate(calibration_point = cal_pts$Name) %>% 
  dplyr::select(cp = calibration_point,
                x = X,
                y = Y)

## Join and add NA values
dets_w <- dets_w %>% 
  dplyr::left_join(cal_pt_location_df,
                   by = "cp") %>% 
  dplyr::mutate_at(vars(contains("node_")), ~replace(., is.na(.), -115)) %>% 
  ungroup()



## Empty data frame for results
output <- data.frame()

## Test 
for (cp_f in unique(dets_w$cp)){
  print(cp_f)
  # # ## Focal cp
  # cp_f <- "N_22"
  # cps_tr <- sample(dets_w$cp, length(dets_w$cp)*0.70)
  
  cp_f_dets <- dets_w[dets_w$cp == cp_f,]
  
  tags <- unique(cp_f_dets$tag)
  
  for(tag_f in tags){
    
    # tag_f = tags[1]
    
    ## Train and test
    dets_tr_x <- dets_w %>% 
      dplyr::filter(!(cp %in% cp_f))  %>% 
      dplyr::select(-cp,
                    -tag,
                    -y)
    
    dets_ts_x <- dets_w %>% 
      dplyr::filter((cp %in% cp_f) & tag == tag_f)  %>% 
      dplyr::select(-cp,
                    -tag,
                    -y)  
    
    #n_features <- length(names(dets_w))
    # rf_x <- ranger(
    #   x ~ ., 
    #   data = dets_tr_x,
    #   mtry = floor(n_features / 3),
    #   respect.unordered.factors = "order",
    #   seed = 123
    # )
    xgb_x <- xgboost::xgboost(
      as.matrix(dets_tr_x[,-ncol(dets_tr_x)]),
      dets_tr_x$x,
      nrounds = 100, max_depth = 6
    )
    ## Predict
    x_est <- predict(xgb_x,
                     as.matrix(dets_ts_x[,-ncol(dets_ts_x)]))
    
    ## Predict y 
    dets_tr_y <- dets_w %>% 
      dplyr::filter(!(cp %in% cp_f)) %>% 
      dplyr::select(-cp,
                    -tag,
                    -x)
    
    dets_ts_y <- dets_w %>% 
      dplyr::filter((cp %in% cp_f) & tag == tag_f) %>% 
      dplyr::select(-cp,
                    -tag,
                    -x)
    
    # n_features <- length(names(dets_w))
    # rf_y <- ranger(
    #   y ~ .,
    #   data = dets_tr_y,
    #   mtry = floor(n_features / 3),
    #   respect.unordered.factors = "order",
    #   seed = 123
    # )
    # 
    # 
    # y_est <- predict(rf_y,
    #                  dets_ts_y)$predictions
    xgb_y <- xgboost::xgboost(
      as.matrix(dets_tr_y[,-ncol(dets_tr_y)]),
      dets_tr_y$y,
      nrounds = 100, max_depth = 6
    )
    ## Predict
    y_est <- predict(xgb_y,
                     as.matrix(dets_ts_y[,-ncol(dets_ts_y)]))
    
    error <- raster::pointDistance(c(x_est, y_est), 
                                   c(dets_ts_x$x, dets_ts_y$y), 
                                   lonlat = F,
                                   allpairs = F)
    
    ## Get error for calibration point
    est_cp_f_df <- data.frame(cp = cp_f,
                              tag = tag_f,
                              x_cal = dets_ts_x$x,
                              y_cal = dets_ts_y$y,
                              x_est = x_est,
                              y_est = y_est,
                              error = error)
    
    ## Bind rows
    output <- bind_rows(output,
                        est_cp_f_df)
  }
  
  
  
}


median(output$error)
hist(output$error, breaks = 100)

## Save outputs
readr::write_csv(output,
                 "xgb_location_estimates.csv")


output <- readr::read_csv("./data/outputs/rf_location_estimates.csv")


## Read node locations
nodes_df <- readr::read_csv("wes-tracking-test-main/data/locations/WES_node_locations.csv")


## Map
output %>% 
  
  ggplot() +
  
  ## Estimated points
  geom_point(aes(x = x_est,
                 y = y_est)) +
  
  # ## Node points
  # geom_text(aes(x = x, 
  #               y = y,
  #               label = grid_point), 
  #           data = nodes_df,
  #           color = "red") +
  # 
  theme_minimal() +
  scale_colour_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous"))


ggplot(output) +
  
  geom_point(aes(x = x_est,
                 y = y_est),
             color = "red") +
  
  geom_point(aes(x = x_cal,
                 y = y_cal),
             color = "black") +
  
  geom_segment(aes(x = x_cal,
                   y = y_cal,
                   xend = x_est,
                   yend = y_est,
                   color = error)) +
  
  scale_colour_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) +
  theme_minimal()
