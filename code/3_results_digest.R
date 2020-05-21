##################################################
## Project: nlktsoil
## Script purpose: Reshape output from model analysis
## Date: Wed May 20 08:54:48 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
library(tidyverse)
library(reshape2)
rm(list = ls())
outputFolder <- "output/"
project <- "biome"
title = "word_pca"

fileFolder <- paste(outputFolder, project, sep = "")
csvs_files <- list.files(fileFolder, pattern = "_stats.*.csv")
custom_csvs_files <- list.files(fileFolder, pattern = "customstats.*.csv")

csvs <- list()
custom_csvs <- list()
for(i in 1:length(csvs_files)){
  csvs[[i]] <- read_csv(paste(fileFolder,"/",csvs_files[i], sep = ""))
  custom_csvs[[i]] <- read_csv(paste(fileFolder,"/",custom_csvs_files[i], sep = ""))
}
names(csvs) <- csvs_files
names(custom_csvs) <- custom_csvs_files
res <- bind_rows(lapply(csvs, as.data.frame), .id = "id") %>%
  mutate(id = gsub(".csv", "", id)) %>%
  mutate(.metric_estimator = paste(.metric, .estimator, sep = "_")) %>%
  select(c(id, .metric_estimator, .estimate)) %>%
  pivot_wider(names_from = .metric_estimator,values_from =  .estimate) 

write_csv(res[c(2,4,1,3),], paste(outputFolder, title, "_stats.csv", sep = ""))
res %>% filter(grepl("test",id)) %>%
  pivot_longer(cols = contains("_"), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = metric, y = value, fill = factor(id))) + geom_bar(stat = "identity", position = "dodge")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(outputFolder, title, "_barplot_comparison.svg", sep = ""),  width = 9, height = 6)  

names(custom_csvs) <- custom_csvs_files
custom_res <- bind_rows(lapply(custom_csvs, as.data.frame), .id = "id") %>%
  rename(biome = X1) %>%
  mutate(id = gsub("_.*_", "", id )) %>%
  mutate(id = gsub(".csv", "", id)) %>%
  mutate(biome = as.factor(biome)) %>%
  mutate(biome = gsub("average", "z_average", biome)) %>%
  filter(grepl("test", id)) %>%
  pivot_longer(cols = colnames(custom_csvs[[1]])[2:9], names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = biome, y = value, fill = as.factor(biome))) + 
  geom_bar(stat = "identity") + 
  facet_grid(id~metric,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
custom_res
ggsave(paste(outputFolder, title, "_custom_barplot_comparison.svg", sep = ""),  width = 9, height = 6)  

