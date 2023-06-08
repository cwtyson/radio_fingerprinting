## Script to filter detections from calibration points

## Housekeeping ######
library(tidyverse)

## Process calibration detections ######
dets <- readr::read_csv("./data/detections/raw_dets.csv")

## Calibration tags
tags <- readr::read_csv("./data/field/test_tag_list.csv")

## Merge pole with dets
dets <- dets %>% 
  dplyr::left_join(tags, 
                   by = "tag")

## Field calibration log
field_log_intervals <- readr::read_csv("./data/field/WES_field_calibration_log.csv") %>% 
  janitor::clean_names() %>% 
  transmute(point,
            pole = tag_pole,
            ## Adjust start time by 1 minute
            start_time = lubridate::dmy_hms(paste(date, format(start_time, "%H:%M"))) + 60,
            end_time = lubridate::dmy_hms(paste(date, format(end_time, "%H:%M"))),
            interval = lubridate::interval(start_time, end_time)) %>% 
  na.omit() %>% 
  filter(lubridate::int_length(interval) >= 60)


## Calibration periods 
for(p in unique(field_log_intervals$pole)) {
  assign(paste0("calibration_periods_", p), as.list(field_log_intervals[field_log_intervals$pole == p,]$interval))
}

## For each pole number, filter detections based on calibration periods
field_cal_dets <- data.frame()
for(p in unique(field_log_intervals$pole)){
  period <- get(paste0("calibration_periods_", p))
  
  ## Filter
  dets_p <- dets  %>% 
    dplyr::filter((date_time %within% period) & (pole == p)) 
  
  ## Combine
  field_cal_dets <- bind_rows(field_cal_dets,
                      dets_p)
  
  rm(dets_p)
  
}

## Reformat
field_cal_dets <- field_cal_dets %>%
  dplyr::mutate(date_time_2 = date_time) %>% 
  dplyr::select(node, pole, tag, rssi, date_time, date_time_2) %>% 
  data.table::data.table() 

## For each pole, join calibration point based on date_time overlap
field_cal_dets_j <- data.frame()
for(p in unique(field_cal_dets$pole)){
  
  ## Subset pole  
  field_cal_dets_p <- field_cal_dets[field_cal_dets$pole == p,]
  
  ## Set up field log intervals as data.table to join
  field_log_intervals_p <- field_log_intervals[field_log_intervals$pole == p,] %>% 
    dplyr::select(calibration_point = point,
                  date_time = start_time,
                  date_time_2 = end_time) %>% 
    data.table::data.table()
  data.table::setkey(field_log_intervals_p, date_time, date_time_2)
  
  
  ## Join based on overlap
  field_cal_dets_p_j <- data.table::foverlaps(x = field_cal_dets_p, 
                                              y = field_log_intervals_p, 
                                              by.x = c("date_time", "date_time_2"), 
                                              by.y = c("date_time", "date_time_2"),
                                              type = "within", 
                                              mult = "all", 
                                              nomatch = 0L) %>% 
    
    ## Reformat
    dplyr::select(calibration_point,
                  pole,
                  node,
                  tag,
                  RSSI = rssi,
                  date_time = i.date_time) 
  
  ## Combine
  field_cal_dets_j <- dplyr::bind_rows(field_cal_dets_j, field_cal_dets_p_j)
  
}


## Reformat 
field_cal_dets <- field_cal_dets_j %>% 
  dplyr::mutate(cal_pt_f = factor(calibration_point, 
                            levels = stringr::str_sort(unique(calibration_point), numeric = T))) %>% 
  dplyr::arrange(cal_pt_f) %>% 
  dplyr::select(-cal_pt_f)


## Summary
field_cal_sum <- field_cal_dets %>%
  group_by(calibration_point) %>%
  summarise(dets = n(),
            tags = n_distinct(tag),
            nodes = n_distinct(node))

## Save 
readr::write_csv(field_cal_dets,
                 "./data/detections/field_cal_dets.csv")
