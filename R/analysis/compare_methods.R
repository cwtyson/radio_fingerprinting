## Compare methods

# Housekeeping ------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(ggmap)

## Colors
color_pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")


# Summary -----------------------------------------------------------------
est_all <- readr::read_csv("./data/outputs/combined_location_estimates.csv")

## Total calibration points
n_distinct(est_all$cp)

## Calibration points per method
n_distinct(est_all[est_all$est == "ml",]$cp)
n_distinct(est_all[est_all$est == "knn",]$cp)

## Total localizations
nrow(est_all[est_all$est == "ml",])
nrow(est_all[est_all$est == "knn",])

# Overall error #######

est_all <- readr::read_csv("./data/outputs/combined_location_estimates.csv") %>% 
  mutate(est = case_when(est == "knn" ~ "Fingerprinting",
                         est == "ml" ~ "Multilateration"))

## Overall
median(est_all[est_all$est == "Multilateration",]$error)
median(est_all[est_all$est == "Fingerprinting",]$error)


# Error when not detected by nearest node ---------------------------------
median(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == FALSE,]$error)
nrow(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == FALSE,])
median(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == FALSE,]$error)
nrow(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == FALSE,])

## When detected by nearest node
median(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == TRUE,]$error)
median(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == TRUE,]$error)

## Error if detected by nearest node
ggplot() +
  
  geom_boxplot(aes(x = est,
                   y = error,
                   fill = nearest_gp_det),
               position = position_dodge2(width = 1),
               alpha = 0.7,
               notch = TRUE,
               data = est_all,
               outlier.alpha = 0.3) +
  geom_text(aes(x = est,
                y = 375,
                label = points),
            position = position_dodge2(width = 1),
            size = 7,
            data = est_all %>%
              group_by(est,nearest_gp_det) %>%
              summarise(points = n())) +
  labs(x = "Localization method",
       y = "Error (m)") +
  scale_fill_manual(name = "Detected by\nnearest receiver", values = color_pal[c(90,20)]) +
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = seq(0,400, by = 50), limits = c(0,400)) 


#  Error with varying amounts of node test data for RF method -------------

## Reduced nodes in RF method: file size exceeds GitHub maximum and must be downloaded separately from Dryad
rf_nodes <- readr::read_csv("./data/outputs/altered_node_test_percentages.csv")

## Detections
est_all <- readr::read_csv("./data/outputs/combined_location_estimates.csv") %>% 
  mutate(est = case_when(est == "knn" ~ "Fingerprinting",
                         est == "ml" ~ "Multilateration"))

est_all_sum <- est_all %>% 
  group_by(est, nearest_gp_det) %>% 
  summarise(median = median(error)) %>% 
  filter(est == "Multilateration")

(rf_nodes_plot <- rf_nodes %>%
    distinct(tag,gp_perc, rep, cp, .keep_all = T) %>% 
    group_by(est_method,
             gp_perc,
             rep) %>%
    summarise(median = median(error)) %>%
    mutate(est_method = case_when(est_method == "knn" ~ "Fingerprinting",
                                  est_method == "ml" ~ "Multilateration")) %>% 
    
    ggplot() +
    geom_hline(yintercept = est_all_sum[est_all_sum$nearest_gp_det == "TRUE",]$median, linetype = 5) +
    geom_hline(yintercept = est_all_sum[est_all_sum$nearest_gp_det == "FALSE",]$median, linetype = 2) +
    geom_jitter(aes(x = gp_perc,
                    y = median),
                height = 0,
                width = 0.01,
                alpha = 0.6) +
    stat_smooth(aes(x = gp_perc,
                    y = median),
                color = color_pal[90]) +
    labs(x = "Proportion of nodes in test data", y = "Median error (m)") +
    theme_classic(base_size = 16)) +
  scale_y_continuous(breaks = seq(0,450, by = 50), limits = c(0,450)) +
  scale_x_reverse() 


# Error when half of the receivers are used -------------------------------

## Altered nodes
alt_node_spacing <- readr::read_csv("./data/outputs/altered_node_spacing.csv") %>% 
  transmute(est_method,
            error,
            nodes = "Subset",
            node_count) %>% 
  mutate(node_count = ifelse(est_method == "ml", NA, node_count))

## All nodes
all_est <- readr::read_csv("./data/outputs/combined_location_estimates.csv") %>% 
  filter(!grepl("WNC", cp)) %>% 
  transmute(est_method = est,
            error,
            nodes = "All")  %>% 
  na.omit()

## Bind
node_spacing <- dplyr::bind_rows(alt_node_spacing,
                                 all_est) 

## Summarize
node_spacing_error <- node_spacing %>% 
  group_by(est_method, nodes) %>% 
  summarise(median = median(error))

## Plot
node_spacing %>%
  # mutate(nodes = ifelse(nodes %in% c("odd","even"), "Subset", "All")) %>% 
  mutate(est_method = case_when(est_method == "knn" ~ "Fingerprinting",
                                est_method == "ml" ~ "Multilateration")) %>% 
  group_by(est_method, nodes) %>%
  mutate(count = n(),
         median = median(error),
         sd = sd(error),
         nodes = stringr::str_to_title(nodes)) %>%
  
  ggplot() +
  
  geom_boxplot(aes(x = nodes,
                   y = error),
               notch = TRUE,
               alpha = 0.7,
               fill = color_pal[20]) +
  
  theme_classic(base_size = 16) +
  facet_grid(.~est_method) +
  xlab("Receivers") +
  ylab("Error (m)") +
  scale_y_continuous(breaks = seq(0,400, by = 50), limits = c(0,400)) 


# Distance from the nearest node center of the grid -----------------------

## Read in estimates
est_all <- readr::read_csv("./data/outputs/combined_location_estimates.csv") 

## Plot
est_all  %>% 
  
  ## Reformat for plotting
  mutate(est = case_when(est == "knn" ~ "Fingerprinting",
                         est == "ml" ~ "Multilateration")) %>% 
  
  ## Filter cases where nearest gp detected 
  filter(nearest_gp_det == TRUE) %>% 
  rename("Nearest receiver" = "distance_gp", 
         "Receiver grid center" = "dist_center") %>% 
  pivot_longer(cols = c("Nearest receiver", "Receiver grid center"),
               names_to = "metric",
               values_to = "distance") %>% 
  ggplot(aes(x = distance,
             y = error)) +
  
  geom_point(alpha = 0.2) +
  
  stat_smooth(method = "lm",
              linetype = 1,
              size = 1,
              color = color_pal[90]) +
  facet_grid(est~metric,
             scales = "free") +
  
  labs(x = "Distance (m)",
       y = "Error (m)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")


# Error relative to habitat -----------------------------------------------

## Read in combined estimates
est_all <- readr::read_csv("./data/outputs/combined_location_estimates.csv") %>% 
  filter(!grepl("WNC", cp)) 

## Field log with habitat information
field_log <- read.csv("./data/field/WES_field_calibration_log.csv") %>% 
  janitor::clean_names() %>% 
  transmute(cp = point,
            habitat) %>% 
  na.omit()

error_sum <- est_all %>% 
  left_join(field_log) %>% 
  mutate(est = case_when(est == "knn" ~ "Fingerprinting",
                         est == "ml" ~ "Multilateration")) %>% 
  
  filter(habitat %in% c(1,2,3,4,5)) %>% 
  
  mutate(habitat = case_when(habitat == 1 ~ 5,
                             habitat == 2 ~ 4,
                             habitat == 3 ~ 3,
                             habitat == 4 ~ 2,
                             habitat == 5 ~ 1))  

error_sum_table <- error_sum %>% 
  group_by(est,habitat) %>% 
  summarise(error_mean = median(error))


## Habitat and error
est_all %>% 
  left_join(field_log) %>% 
  mutate(est = case_when(est == "knn" ~ "Fingerprinting",
                         est == "ml" ~ "Multilateration")) %>% 
  
  group_by(est,habitat) %>% 
  mutate(median = median(error)) %>% 
  
  filter(habitat %in% c(1,2,3,4,5)) %>% 
  
  mutate(habitat = case_when(habitat == 1 ~ 5,
                             habitat == 2 ~ 4,
                             habitat == 3 ~ 3,
                             habitat == 4 ~ 2,
                             habitat == 5 ~ 1)) %>% 
  ggplot() +
  geom_boxplot(aes(x = as.factor(habitat), 
                   y = error),
               notch = TRUE,
               alpha = 0.7,
               fill = color_pal[20],
               outlier.alpha = 0.4) +
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = seq(0,400, by =50), limits = c(0,400)) +
  facet_grid(~est) +
  labs(y = "Error (m)",
       x = "Habitat score") +
  scale_fill_manual(values = color_pal[c(20,90)], name = "Method")


# Reduced time at each calibration point  ---------------------------------

## Check outputs
alt_cp_time <- readr::read_csv("./data/outputs/altered_calibration_point_times.csv") 

## Habitats
field_log <- read.csv("./data/field/WES_field_calibration_log.csv") %>% 
  janitor::clean_names() %>% 
  transmute(cp = point,
            habitat) %>% 
  na.omit() %>% 
  mutate(habitat = case_when(habitat == 1 ~ 5,
                             habitat == 2 ~ 4,
                             habitat == 3 ~ 3,
                             habitat == 4 ~ 2,
                             habitat == 5 ~ 1)) %>% 
  mutate(habitat = as.character(habitat))


## Join habitat types and plot detections by time
fig_rt_a <- alt_cp_time %>% 
  left_join(field_log) %>%
  filter(habitat %in% c(1,2,3,4,5)) %>% 
  
  ## Plot
  ggplot() +
  stat_smooth(aes(x = cp_time,
                  y = dets,
                  group = habitat,
                  color = habitat),
              fill = grey(0.7),
              formula = y ~ x,
              method = "lm") +
  scale_x_continuous(breaks = seq(5,60, by =5)) +
  labs(x = "Time (seconds) at test point ", 
       y = "Detections per tag") +
  scale_color_manual(name = "Habitat score", values = color_pal[seq(1,100, length.out = 5)]) +
  theme_minimal()  +
  theme_classic(base_size = 16) 

## Join habitat types and plot error by time
fig_rt_b <- alt_cp_time %>% 
  left_join(field_log) %>% 
  filter(habitat %in% c(1,2,3,4,5)) %>% 
  
  ## Plot
  ggplot() +
  stat_smooth(aes(x = cp_time,
                  y = error,
                  group = habitat,
                  color = habitat),
              fill = grey(0.8)) +
  scale_x_continuous(breaks = seq(5,60, by =5)) +
  labs(x = "Time (seconds) at test point ", y = "Median error (m)") +
  scale_color_manual(name = "Habitat score", values = color_pal[seq(1,100, length.out = 5)]) +
  theme_classic(base_size = 16) 

ggarrange(fig_rt_a,fig_rt_b, common.legend = TRUE, legend = "bottom",labels = c("A","B"))


# Compare rounds  ---------------------------------------------------------

## KNN estimates
br <- read_csv("./data/outputs/knn_rounds_combined.csv") %>% 
  filter(!grepl("WNC", cp))

median(br[br$round == "one",]$error)
median(br[br$round == "two",]$error)

# Compare ML methods  ----------------------------------------------------------

## Multilateration
est_all <- readr::read_csv("./data/outputs/combined_location_estimates.csv") %>% 
  mutate(method = case_when(est == "knn" ~ "Fingerprinting",
                            est == "ml" ~ "Multilateration")) %>% 
  filter(method == "Multilateration") %>% 
  select(method,
         error)

## Random forests
rf_df <- readr::read_csv("./data/outputs/rf_location_estimates.csv") %>% 
  mutate(method = "Random Forest")

## KNN
knn_df <- readr::read_csv("./data/outputs/knn_k_comparison_estimates.csv")  %>% 
  mutate(method = paste0("KNN (", k, ")")) %>% 
  filter(!grepl("WNC", cp))

## SVM
svm_df <- readr::read_csv("./data/outputs/svm_location_estimates.csv")  %>% 
  mutate(method = "SVM")

## xgboost
xgboost_df <- readr::read_csv("./data/outputs/xgb_location_estimates.csv")  %>% 
  mutate(method = "XGBoost")


## Combine
error_dfs <- bind_rows(rf_df,
                       knn_df,
                       svm_df,
                       xgboost_df,
                       est_all) %>% 
  select(-tag,
         -node_count)

## Summarise
error_sum <- error_dfs %>% 
  group_by(method) %>% 
  summarise(median = round(median(error))) %>% 
  arrange(median)

## Order
levels <- unique(error_sum$method)

error_dfs$method <- factor(error_dfs$method,
                           levels = levels)

## Plot
ggplot(error_dfs) +
  geom_boxplot(aes(x = method, 
                   y = error),
               notch = TRUE,
               alpha = 0.7,
               fill = color_pal[20],
               outlier.alpha = 0.4) +
  geom_text(aes(x = method,
                y = 300,
                label = median),
            position = position_dodge2(width = 1),
            size = 5,
            data = error_sum) +
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = seq(0,300, by =50), limits = c(0,300)) +
  labs(y = "Error (m)",
       x = "Localization method") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8))
