##################################################
## Project: nltksoil
## Script purpose:  plot figures from model analysis
## Date: Tue May 26 12:56:37 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
library(tidyverse)
library(ggplot2)
rm(list = ls())
outputFolder <- "output/"
source("code/nltksoil_functions.R")
## Section B: Read Data 
##################################################

custom_stats_ <- list.files(path =  outputFolder, pattern = "customstats.*")
stats <- list.files(path =  outputFolder, pattern = "_stats.*")
roc <- list.files(path =  outputFolder, pattern = "roc.*")
pr <- list.files(path =  outputFolder, pattern = "pr.*")
term <- "customstats"
patternname <- paste(term, ".*", sep = "")
files <- list.files(path = outputFolder, pattern = patternname)
files_list <- list()
for(i in 1:length(files)){
  files_list[[i]] <- read_csv(paste(outputFolder, files[i], sep = ""))
}
names(files_list) <- files



#custom stats
custom_stats_data <- bind_rows(files_list, .id = ".property") %>%
  pivot_longer(cols = c("TP", "FP", "TN", "FN", "PPV", "NPV", "sensitivity", "specificity",
                             "detection_rate", "detection_prevalence", "balance_accuracy"), 
                             names_to = ".statistic", values_to = ".value") %>%
  mutate(.property = gsub("rf_|_envo_biome_2_customstats.csv|[0-9]", "", .property)) %>%
  filter(.statistic != "FP" & .statistic != "TP" & .statistic != "FN" & .statistic != "TN") %>%
  mutate(.variable = factor(.variable, levels = c("average", "anthropogenic_terrestrial", "desert", "forest", "grassland", "shrubland", "tundra"))) 
p1 <- ggplot(filter(custom_stats_data, .property != "abundance"), aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = CI_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free"); p1
ggsave(plot = p1, "output/custom_stats_CI_comparison.svg", height = 8, width = 8)

p2 <- ggplot(filter(custom_stats_data), aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = CI_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free"); p2
ggsave(plot = p2, "output/custom_stats_CI_all.svg", height = 8, width = 8)

p3 <- ggplot(filter(custom_stats_data, .property != "abundance"), aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = CI2_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free"); p3
ggsave(plot = p3, "output/custom_stats_CI2_comparison.svg", height = 8, width = 8)

p4 <- ggplot(filter(custom_stats_data), aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = CI2_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free"); p4
ggsave(plot = p5, "output/custom_stats_CI2_all.svg", height = 8, width = 8)

p5 <- ggplot(filter(custom_stats_data, .property != "abundance"), aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = quantiles_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free"); p5
ggsave(plot = p5, "output/custom_stats_quantile_comparison.svg", height = 8, width = 8)

p6 <- ggplot(filter(custom_stats_data), aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = quantiles_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free"); p6
ggsave(plot = p6, "output/custom_stats_quantile_all.svg", height = 8, width = 8)

#custom stats
custom_stats <- bind_rows(files_list, .id = ".property") %>%
  pivot_longer(cols = c("TP", "FP", "TN", "FN", "PPV", "NPV", "sensitivity", "specificity",
                        "detection_rate", "detection_prevalence", "balance_accuracy"), 
               names_to = ".statistic", values_to = ".value") %>%
  mutate(.property = gsub("rf_|_envo_biome_2_customstats.csv|[0-9]", "", .property)) %>%
  filter(.statistic != "FP" & .statistic != "TP" & .statistic != "FN" & .statistic != "TN") %>%
  #filter(.property != "abundance") %>%
  mutate(.variable = factor(.variable, levels = c("average", "anthropogenic_terrestrial", "desert", "forest", "grassland", "shrubland", "tundra"))) %>%
  ggplot(aes(x= .property, y= .value, fill = .property)) +
  stat_summary(fun.data = CI_95, geom="boxplot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.statistic~.variable, scales = "free")
custom_stats
ggsave(plot = custom_stats, "output/custom_stats_CI_all.svg", height = 8, width = 8)

## Section C: Plot Figures  
##################################################
term <- "roc"
patternname <- paste(term, ".*", sep = "")
files <- list.files(path = outputFolder, pattern = patternname)
files_list <- list()
for(i in 1:length(files)){
  files_list[[i]] <- read_csv(paste(outputFolder, files[i], sep = ""))
}
names(files_list) <- files
roc_data <- bind_rows(files_list, .id = ".property") 
r1 <-roc_data %>%
  filter(.set == "test") %>%
  mutate(.fold = as.factor(.fold)) %>%
  mutate(.property = gsub("rf_|_envo_biome_2_roc.csv|[0-9]", "", .property)) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = .fold, type = .sample)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.property~.level)
ggsave(plot = r1, "output/roc.svg", width = 6, height = 3)

term <- "pr"
patternname <- paste(term, ".*", sep = "")
files <- list.files(path = outputFolder, pattern = patternname)
files_list <- list()
for(i in 1:length(files)){
  files_list[[i]] <- read_csv(paste(outputFolder, files[i], sep = ""))
}
names(files_list) <- files
pr_data <- bind_rows(files_list, .id = ".property") 
r2 <- pr_data %>%
  filter(.set == "test") %>%
  mutate(.fold = as.factor(.fold)) %>%
  mutate(.property = gsub("rf_|_envo_biome_2_roc.csv|[0-9]", "", .property)) %>%
  ggplot(aes(x = recall, y = precision, color = .fold, type = .sample)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.property~.level); r2
ggsave(plot = r2, "output/pr.svg", width = 6, height = 3)
