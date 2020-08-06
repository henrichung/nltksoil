##################################################
## Project: nltksoil
## Script purpose:  Test out random forest predictions
## Date: Tue Jun 23 08:43:24 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("tidyverse","tidymodels",  "reshape2", "furrr", "future", "doParallel", 
                      "text2vec", "stringi")#
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))
rm(list = ls())
source("code/nltksoil_functions.R") #Load functions

outputFolder <- "../../../../work/idoerg/hchung/outputs/soil"
##random forest and 50 dimensions
results_files <- list.files(outputFolder)
results_filenames <- paste(outputFolder, results_files, sep = "/")
results_list_tune <- list(); 
results_list_test <- list(); 
results_list_models <- list()
results_list_best <- list()
results_filenames
for(i in 1:length(results_files)){
  temp <- readRDS(results_filenames[i])
  results_list_tune[[i]] <- temp$tune
  results_list_test[[i]] <- temp$test
  results_list_models[[i]] <- temp$model
  results_list_best[[i]] <- temp$best
  if(i %% 5 == 0){message(i)}
}
names(results_list_tune) <- results_files; names(results_list_test) <- results_files
names(results_list_best) <- results_files; names(results_list_models) <- results_files
results_tune_raw <- bind_rows(results_list_tune, .id = "filename")
results_test_raw <- bind_rows(results_list_test, .id = "filename")
results_best_raw <- bind_rows(results_list_best, .id = "filename")

saveRDS(results_list_models, "../../../../work/idoerg/hchung/outputs/soil/results_models.RDS")

error_multi_metric <- function (x) {
  multi_metric <- metric_set(accuracy, bal_accuracy, f_meas, kap, npv, ppv, precision, recall, sens, pr_auc, roc_auc)
  levels = paste(as.character(levels(x$response)), collapse = "|")
  res <- x %>%
  group_by(mtry, min_n) %>% 
  multi_metric(truth = "response", estimate = ".pred_class", ... = matches(levels), na_rm = TRUE, estimator = "macro_weighted")
  #res <- multi_metric(x, truth = "response", estimate = ".pred_class", ... = matches(levels), na_rm = TRUE, estimator = "macro_weighted")
  return(res)
}

error_multi_metric2 <- function (x) {
  multi_metric <- metric_set(accuracy, bal_accuracy, f_meas, kap, npv, ppv, precision, recall, sens)
  levels = paste(as.character(levels(x$response)), collapse = "|")
  res <- multi_metric(x, truth = "response", estimate = "estimate", ... = matches(levels), na_rm = TRUE, estimator = "macro_weighted")
  return(res)
}


plan(multiprocess)

results_tune <- results_tune_raw %>%
  mutate(filename = gsub("envo_biome_2", "envobiome2", filename)) %>%
  separate(filename, c(".taxa", ".response", ".model", ".no_components", ".method")) %>%
  mutate(.method = factor(.method, levels = c("base", "abundance", "transformpca", "transformword", "pca", "word"))) %>%
  mutate(.no_components = factor(.no_components, levels = c(50, 100, 150, 200, 250, 300))) %>%
  mutate(.taxa = factor(.taxa, levels = c("species", "genus","family", "order", "class", "phylum"))) %>%
  mutate(.metrics = future_map(.predictions, error_multi_metric)) %>%
  mutate(.custom = future_map(.predictions, custom_evaluate)) %>%
  mutate(.pr_curve = future_map(.predictions, custom_pr_curve)) %>%
  mutate(.roc_curve = future_map(.predictions, custom_roc_curve)) %>%
  mutate(.null = future_map(.null, ~mutate(., estimate = as.factor(estimate)))) %>%
  mutate(.null = future_map(.null, ~refactor(., "response", "estimate"))) %>%
  mutate(.nullmetrics = future_map(.null, error_multi_metric2))

saveRDS(results_tune, "results_tune.RDS")
 

results_best <- results_best_raw %>%
  mutate(filename = gsub("envo_biome_2", "envobiome2", filename)) %>%
  separate(filename, c(".taxa", ".response", ".model", ".no_components", ".method")) %>%
  mutate(.method = factor(.method, levels = c("base", "abundance", "transformpca", "transformword", "pca", "word"))) %>%
  mutate(.no_components = factor(.no_components, levels = c(50, 100, 150, 200, 250, 300))) %>%
  mutate(.taxa = factor(.taxa, levels = c("species", "genus","family", "order", "class", "phylum"))) 
saveRDS(results_best, "results_best.RDS")

results_best2 <- results_best_raw %>% group_by(filename) %>% slice(which.max(mean))

results_list_train <- list()
for(i in 1:length(results_list_tune)){
  temp <- results_list_tune[[i]]
  bmtry = results_best2[i,] %>% pull(mtry)
  bmin_n = results_best2[i,] %>% pull(min_n)
  results_list_train[[i]] <-  temp %>% mutate(.predictions = future_map(.predictions, ~filter(., mtry == bmtry, min_n == bmin_n)))
}
names(results_list_train) <- results_files
results_train_raw <- bind_rows(results_list_train, .id = "filename")

results_train <- results_train_raw %>%
  mutate(filename = gsub("envo_biome_2", "envobiome2", filename)) %>%
  separate(filename, c(".taxa", ".response", ".model", ".no_components", ".method")) %>%
  mutate(.method = factor(.method, levels = c("base", "abundance", "transformpca", "transformword", "pca", "word"))) %>%
  mutate(.no_components = factor(.no_components, levels = c(50, 100, 150, 200, 250, 300))) %>%
  mutate(.taxa = factor(.taxa, levels = c("species", "genus","family", "order", "class", "phylum"))) %>%
  mutate(.metrics = future_map(.predictions, error_multi_metric)) %>%
  mutate(.custom = future_map(.predictions, custom_evaluate)) %>%
  mutate(.pr_curve = future_map(.predictions, custom_pr_curve)) %>%
  mutate(.roc_curve = future_map(.predictions, custom_roc_curve)) %>%
  mutate(.null = future_map(.null, ~mutate(., estimate = as.factor(estimate)))) %>%
  mutate(.null = future_map(.null, ~refactor(., "response", "estimate"))) %>%
  mutate(.nullmetrics = future_map(.null, error_multi_metric2))
saveRDS(results_train, "results_train.RDS")


results_test <- results_test_raw %>%
  mutate(filename = gsub("envo_biome_2", "envobiome2", filename)) %>%
  separate(filename, c(".taxa", ".response", ".model", ".no_components", ".method")) %>%
  mutate(.method = factor(.method, levels = c("base", "abundance", "transformpca", "transformword", "pca", "word"))) %>%
  mutate(.no_components = factor(.no_components, levels = c(50, 100, 150, 200, 250, 300))) %>%
  mutate(.taxa = factor(.taxa, levels = c("species", "genus","family", "order", "class", "phylum"))) %>%
  mutate(.predictions = future_map(.predictions, ~mutate(., response = droplevels(response), .pred_class = droplevels(.pred_class)))) %>%
  mutate(.predictions = future_map(.predictions, ~refactor(., "response", ".pred_class"))) %>%
  mutate(.metrics = future_map(.predictions, error_multi_metric)) %>%
  mutate(.custom =future_map(.predictions, custom_evaluate)) %>%
  mutate(.pr_curve = future_map(.predictions, custom_pr_curve)) %>%
  mutate(.roc_curve = future_map(.predictions, custom_roc_curve)) %>%
  mutate(.null = future_map(.null, ~mutate(., estimate = as.factor(estimate)))) %>%
  mutate(.null = future_map(.null, ~refactor(., "response", "estimate"))) %>%
  mutate(.nullmetrics = future_map(.null, error_multi_metric2))
saveRDS(results_test, "results_test.RDS")
