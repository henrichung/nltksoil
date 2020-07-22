##################################################
## Project: nltksoil
## Script purpose:  Test out random forest predictions
## Date: Tue Jun 23 08:43:24 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("tidyverse", "parallel", "parsnip", "rsample", "yardstick", "reshape2", 
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel", 
                      "text2vec", "Rcpp", "profvis","Matrix", "themis", "stringi", "data.table", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))
source("code/nltksoil_functions.R") #Load functions

outputFolder <- "../../../../work/idoerg/hchung/outputs/soil"
##random forest and 50 dimensions
results_files <- list.files(outputFolder)
results_filenames <- paste(outputFolder, results_files, sep = "/")
results_list_validation <- list(); results_list_tune <- list()
for(i in 1:length(results_files)){
	temp <- readRDS(results_filenames[i])
	results_list_tune[[i]] <- temp$tuning %>% select(-c("splits"))
	results_list_validation[[i]] <- temp$final %>% select(-c("splits"))
	if(i %% 5 == 0){message(i)}
}
names(results_list_tune) <- results_files; names(results_list_validation) <- results_files
results_tune_raw <- bind_rows(results_list_tune, .id = "filename")
results_validation_raw <- bind_rows(results_list_validation, .id = "filename")

multi_metric <- metric_set(accuracy, bal_accuracy, detection_prevalence, f_meas, j_index, kap, mcc, npv, ppv, precision, recall, sens, spec)
plan(multiprocess)
results_validation <- results_validation_raw %>%
	mutate(filename = gsub("envo_biome_2", "envobiome2", filename)) %>%
	separate(filename, c(".taxa", ".response", ".model", ".no_components", ".method")) %>%
	unnest(.metrics) %>%
	mutate(.method = factor(.method, levels = c("base", "abundance", "transformpca", "transformword", "pca", "word"))) %>%
  	mutate(.no_components = factor(.no_components, levels = c(50, 100, 150, 200, 250, 300))) %>%
  	mutate(.taxa = factor(.taxa, levels = c("species", "genus","family", "order", "class", "phylum"))) %>%
  	unnest(.extracts) %>%
  	rename(.oob_error = .extracts) %>%
  	#unnest(.oob_error) %>%
  	mutate(.custom = future_map(.predictions, multi_metric, estimate = ".pred_class", truth = "response", estimator = "macro_weighted")) %>%
  	select(-c(".metric", ".estimator", ".estimate")) %>%
  	unnest(.custom) %>%
	mutate(.custom = future_map(.predictions, custom_evaluate))
saveRDS(results_validation, "results_validation.RDS")

results_tune <- results_tune_raw %>%
	mutate(filename = gsub("envo_biome_2", "envobiome2", filename)) %>%
	separate(filename, c(".taxa", ".response", ".model", ".no_components", ".method")) %>%
	unnest(.metrics) 
saveRDS(results_tune, "results_tune.RDS")



