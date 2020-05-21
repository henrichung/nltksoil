##################################################
## Project: NLKTsoil
## Script purpose:  tidymodels version
## Date: Mon May 18 15:35:07 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
library(tidyverse)
library(tidymodels)
library(reshape2)
#clear environment
rm(list = ls())
#custom functions
source("code/nltksoil_functions.R")

#parameters
dataFolder <- "data/"
outputFolder <- "output/"
soil_data <- readRDS(paste(dataFolder, "2_soil_data.RDS", sep = ""))
rf_word0 <- readRDS(paste(dataFolder, "2_rf_word.RDS", sep = ""))
rf_pca0 <- readRDS(paste(dataFolder, "2_rf_pca.RDS", sep = ""))
rf_abundance0 <- (readRDS(paste(dataFolder, "2_rf_abundance.RDS", sep = ""))) %>%
  as.data.frame() %>%
  rownames_to_column("X.SampleID")

soil_data <- filter(soil_data, X.SampleID %in% rf_word0$X.SampleID)
{
title <- "word"
#title <- "pca"
#title <- "abundance"
  #join with metadata
rf <- rf_word0 %>%
  left_join(soil_data, by = "X.SampleID") %>%
  mutate(response = envo_biome_2) %>%
  select( matches("V[0-9]{2}|response")) %>%
  mutate(response = gsub(" ", "_", response)) %>%
  mutate(response = gsub("_biome", "", response)) %>%
  mutate(response = as.factor(response))
  
#split into setaside, training, testing
split <-initial_split(rf, prop = 0.9) 
setaside <- testing(split) #10% set aside
current <- initial_split(training(split), prop = 8/9)
training <- training(current) #80% training
testing <- testing(current) #10% testing

#fit random forest model
model <-  rand_forest(mode = "classification") %>%
  set_engine("randomForest") %>%
  fit(response ~ ., data = training) 

#predict on training set
m_train <- predict(model, training) %>% bind_cols(training)
#predict on testing set
m_test <- predict(model,testing) %>% bind_cols(testing)

#probability predictions
mp_test <- predict(model, testing, type = "prob") %>%
  bind_cols(m_test) %>%
  select(-matches("V.{2}")) %>%
  rename(estimate = .pred_class) 

mp_train <- predict(model, training, type = "prob") %>%
  bind_cols(m_train) %>%
  select(-matches("V.{2}")) %>%
  rename(estimate = .pred_class)
##
multi_metric <- metric_set(ppv, npv, sens, spec, precision, recall, accuracy,kap, pr_auc, roc_auc)
stats_train <- multi_metric(mp_train, truth = response, estimate = estimate,... = matches(".pred.*"));stats_train
write_csv(stats_train, paste(outputFolder, title, "_stats_train.csv", sep = ""))
stats_test <- multi_metric(mp_test, truth = response, estimate = estimate, ... = matches(".pred.*")); stats_test
write_csv(stats_test, paste(outputFolder, title, "_stats_test.csv", sep = ""))

custom_stats_train <- custom_evaluate(conf_mat(mp_train, response, estimate)$table); custom_stats_train
write.csv(custom_stats_train, paste(outputFolder, title, "_customstats_train.csv", sep = ""))
custom_stats_test <- custom_evaluate(conf_mat(mp_test, response, estimate)$table); custom_stats_test
write.csv(custom_stats_test, paste(outputFolder, title, "_customstats_test.csv", sep = ""))

##
roc_curve(mp_train, ... = matches(".pred.*"), truth = response) %>% autoplot()
ggsave(paste(outputFolder, title, "_train_roc.svg", sep = ""), width = 9, height = 6)
pr_curve(mp_train, ... = matches(".pred.*"), truth = response) %>% autoplot()
ggsave(paste(outputFolder, title, "_train_pr.svg", sep = ""),  width = 9, height = 6)
roc_curve(mp_test, ... = matches(".pred.*"), truth = response) %>% autoplot()
ggsave(paste(outputFolder, title, "_test_roc.svg", sep = ""),  width = 9, height = 6)
pr_curve(mp_test, ... = matches(".pred.*"), truth = response) %>% autoplot()
ggsave(paste(outputFolder, title, "_test_pr.svg", sep = ""),  width = 9, height = 6)
##
mp_test %>% conf_mat(response, estimate) %>% autoplot(type = "heatmap") 
ggsave(paste(outputFolder, title, "_train_heatmap.svg", sep = ""),  width = 9, height = 9)
mp_train %>% conf_mat(response, estimate) %>% autoplot(type = "heatmap") 
ggsave(paste(outputFolder, title, "_test_heatmap.svg", sep = ""), width = 9, height = 9)
##
}


