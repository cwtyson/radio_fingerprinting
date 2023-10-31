## Process raw detections associated with radio map fingerprint points

## This script provides an a template for associating detections associated with radio fingerprint points. 
## This script requires three files:
## 1) raw_dets.csv: The 'raw' detections (with minimally columns for the receiver/node, tag, rssi, and date time)
## 2) test_tag_list.csv: The tags used to create the radio map. If multiple tags are used at each test point, then a separate column that identifies groups of tags is also need, e.g. the 'pole' number.
## 3) field_log.csv: The log with the times when tags were located at the radio fingerprint (test) points. This file should have columns for the point name, the tag at that point or the tag group (pole), the start date/time for the radio fingerprint point, and the end date/time at that point. See the 'field_log.csv' for an example.

## Housekeeping ######
library(tidyverse)

## Process calibration detections ######

## Insert path to file 1 here
dets <- readr::read_csv() 

## Insert path to file 2 here
tags <- readxl::read_excel() 

## Merge tag group with dets
dets <- dets %>% 
  dplyr::left_join(tags, 
                   by = "tag")

## Insert path to file 3 here
field_log_intervals <- readr::read_csv() 

## Reformat
field_log_intervals <- field_log_intervals%>% 
  janitor::clean_names() %>% 
  transmute(point,
            pole,
            start_time = lubridate::dmy_hms(paste(date, format(start_time, "%H:%M:%S"))),
            end_time = lubridate::dmy_hms(paste(date, format(end_time, "%H:%M:%S"))),
            interval = lubridate::interval(start_time, end_time)) 


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


## Save 
readr::write_csv(field_cal_dets,
                 "./radio_fingerprinting_workflow/fingerprinting_files/fingerprint_detections.csv")


