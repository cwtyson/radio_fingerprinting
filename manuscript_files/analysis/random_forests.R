## Script to estimate calibration locations using random forests

## Housekeeping ######
library(tidyverse)
library(ranger)  
library(caret)

## Plotting values
color_pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

## Process detection data for random forests ######

## Read in detections
dets <- readr::read_csv("./data/detections/field_cal_dets.csv")

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
cal_pts <- sf::read_sf("data/locations/WES_calibration_points.KML") %>% 
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

## Get sample
dets_w_sample <- dets_w %>% 
  slice_sample(prop =0.7)


## Hyperparameter tuning for x
grid <-  expand.grid(mtry = c(3:10),
                     splitrule = c("variance","extratrees","maxstat","beta"),
                     min.node.size = c(5:10))

fitControl <- trainControl(method = "CV",
                           number = 5,
                           verboseIter = TRUE)

## Train on subset
fit = caret::train(
  x = dets_w_sample %>% 
    dplyr::select(-cp,
                  -tag,
                  -x,
                  -y),
  y = dets_w_sample %>% 
    dplyr::select(x) %>% 
    data.frame() %>% 
    pull(),
  method = 'ranger',
  num.trees = 200,
  tuneGrid = grid,
  trControl = fitControl
)

## Best fit for hyperparameters
fit$bestTune

## Empty data frame for results
output <- data.frame()

## Test 
for (cp_f in unique(dets_w$cp)){
  
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
    
    n_features <- length(names(dets_w))
    rf_x <- ranger(
      x ~ ., 
      data = dets_tr_x,
      
      respect.unordered.factors = "order",
      
      ## Hyperparameters based on grid search
      mtry = 9,
      splitrule = "variance",
      min.node.size = 5,
      seed = 123
    )
    
    ## Predict
    x_est <- predict(rf_x,
                     dets_ts_x)$predictions
    
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
    
    n_features <- length(names(dets_w))
    rf_y <- ranger(
      y ~ ., 
      data = dets_tr_y,
      
      ## Hyperparameters based on grid search
      mtry = 9,
      splitrule = "variance",
      min.node.size = 5,
      seed = 123
    )
    
    
    y_est <- predict(rf_y,
                     dets_ts_y)$predictions
    
    
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

## Save outputs
readr::write_csv(output,
                 "./data/outputs/rf_location_estimates.csv")
