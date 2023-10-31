## Compile main results to report

## Housekeeping 
library(tidyverse)
library(ggpubr)
library(wesanderson)

# Summary #######

est_all <- readr::read_csv("./manuscript_files/data/combined_location_estimates.csv")

## Overall
median(est_all[est_all$est == "Multilateration",]$error)
length(est_all[est_all$est == "Multilateration",]$error)
median(est_all[est_all$est == "Fingerprinting",]$error)
length(est_all[est_all$est == "Fingerprinting",]$error)

# Figure 2. Error when not detected by nearest node #####

## Not detected
median(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == FALSE,]$error)
nrow(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == FALSE,])
median(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == FALSE,]$error)
nrow(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == FALSE,])

## Detected
median(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == TRUE,]$error)
length(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == TRUE,]$error)
median(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == TRUE,]$error)
length(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == TRUE,]$error)

## Fingerprinting 
wilcox.test(est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == TRUE,]$error,
            est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == FALSE,]$error)

## Multilateration 
wilcox.test(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == TRUE,]$error,
            est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == FALSE,]$error)

## Fingerprinting to Multilateration - detected
fw_true <- wilcox.test(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == TRUE,]$error,
                       est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == TRUE,]$error)

## Z value
qnorm(fw_true$p.value/2)

## P value
fw_true$p.value

## Fingerprinting to Multilateration - not detected
fw_false <- wilcox.test(est_all[est_all$est == "Multilateration" & est_all$nearest_gp_det == FALSE,]$error,
                        est_all[est_all$est == "Fingerprinting" & est_all$nearest_gp_det == FALSE,]$error)

## Z value
qnorm(fw_false$p.value/2)

## P value
fw_false$p.value

## Plot
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
  scale_fill_manual(name = "Detected by\nnearest receiver", values = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[c(90,20)]) +
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = seq(0,400, by = 50), limits = c(0,400)) 

## Save plot
ggsave("./figures/fig2.jpg",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)

rm(list=ls())

# Figure 3. Distance from the nearest node / center of the grid #####
est_all <- readr::read_csv("./manuscript_files/data/combined_location_estimates.csv") 

## Plot
est_all  %>% 
  
  ## Cases where nearest gp detected 
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
              color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[90]) +
  facet_grid(est~metric,
             scales = "free") +
  
  labs(x = "Distance (m)",
       y = "Error (m)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

## Save plot
ggsave("./figures/fig3.jpg",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)

rm(list=ls())

# Figure 4. Error relative to habitat ########

## Read in combined estimates
est_all <- readr::read_csv("./manuscript_files/data/combined_location_estimates.csv")

ggplot(est_all %>% 
         filter(habitat %in% c(1,2,3,4,5))) +
  
  geom_boxplot(aes(x = as.factor(habitat), 
                   y = error),
               notch = TRUE,
               alpha = 0.7,
               fill = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[20],
               outlier.alpha = 0.4) +
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = seq(0,400, by =50), limits = c(0,400)) +
  facet_grid(~est) +
  labs(y = "Error (m)",
       x = "Habitat score") +
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[c(20,90)], name = "Method")


## Save plot
ggsave("./figures/fig4.jpg",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)


rm(list=ls())

#	Figure 5. Error when half of the receivers are used  #######

## Altered nodes
alt_node_spacing <- readr::read_csv("./manuscript_files/data/altered_node_spacing.csv")

## All nodes
all_est <- readr::read_csv("./manuscript_files/data/combined_location_estimates.csv") %>% 
  transmute(est_method = est,
            error,
            nodes = "All") 

## Conbine
node_spacing <- dplyr::bind_rows(alt_node_spacing,
                                 all_est) 

## Plot
node_spacing %>%
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
               fill = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[20]) +
  
  theme_classic(base_size = 16) +
  facet_grid(.~est_method) +
  xlab("Receivers") +
  ylab("Error (m)") +
  scale_y_continuous(breaks = seq(0,400, by = 50), limits = c(0,400)) 

## Save plot
ggsave("./figures/fig5.jpg",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)

rm(list=ls())

# Figure 6. Error with varying amounts of node test data for RF method ##### 

## Reduced nodes in RF method
rf_nodes <- readr::read_csv( "./manuscript_files/data/altered_node_test_percentages.csv")

## Detections
est_all <- readr::read_csv("./manuscript_files/data/combined_location_estimates.csv") 

## Calculate median error for multilateration
est_all_sum <- est_all %>% 
  group_by(est, nearest_gp_det) %>% 
  summarise(median = median(error)) %>% 
  filter(est == "Multilateration")

## Plot 
rf_nodes_plot <- ggplot(rf_nodes) +
  geom_hline(yintercept = est_all_sum[est_all_sum$nearest_gp_det == "TRUE",]$median, linetype = 5) +
  geom_hline(yintercept = est_all_sum[est_all_sum$nearest_gp_det == "FALSE",]$median, linetype = 2) +
  geom_jitter(aes(x = gp_perc,
                  y = median),
              height = 0,
              width = 0.01,
              alpha = 0.6) +
  stat_smooth(aes(x = gp_perc,
                  y = median),
              color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[90]) +
  labs(x = "Proportion of nodes in test data", y = "Median error (m)") +
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = seq(0,450, by = 50), limits = c(0,450)) +
  scale_x_reverse() 

## Save plot
ggsave("./figures/fig6.jpg",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)

rm(list=ls())

# Figure 7. Reduced time at each calibration point #########

## Check results
alt_cp_time <- readr::read_csv("./manuscript_files/data/altered_calibration_point_times.csv") %>% 
  mutate(habitat = as.factor(habitat))

## Join habitat types and plot detections by time
(fig_rt_a <- ggplot(alt_cp_time) +
    
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
    scale_color_manual(name = "Habitat score", values = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1,100, length.out = 5)]) +
    theme_minimal()  +
    theme_classic(base_size = 16) )

## Plot error by time
(fig_rt_b <-  ggplot(alt_cp_time) +
    stat_smooth(aes(x = cp_time,
                    y = error,
                    group = habitat,
                    color = habitat),
                fill = grey(0.8)) +
    scale_x_continuous(breaks = seq(5,60, by =5)) +
    labs(x = "Time (seconds) at test point ", y = "Median error (m)") +
    scale_color_manual(name = "Habitat score", values = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1,100, length.out = 5)]) +
    theme_classic(base_size = 16))

ggarrange(fig_rt_a,fig_rt_b, common.legend = TRUE, legend = "bottom",labels = c("A","B"))
ggsave("./figures/fig7.jpg",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)

rm(list=ls())

# Supplemental: ML comparisons #######

## Read in errors from the different methods
ml_methods <- readr::read_csv("./manuscript_files/data/ml_comparison.csv")

## Summarise
error_sum <- ml_methods %>% 
  group_by(method) %>% 
  summarise(median = round(median(error))) %>% 
  arrange(median)

## Order
levels <- unique(error_sum$method)

ml_methods$method <- factor(ml_methods$method,
                            levels = levels)

## Plot
ggplot(ml_methods) +
  geom_boxplot(aes(x = method, 
                   y = error),
               notch = TRUE,
               alpha = 0.7,
               fill = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[20],
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

ggsave("./figures/figS1",
       device = "jpg",
       units = "px",
       width = 2304,
       height = 1529)

rm(list=ls())

