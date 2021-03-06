##################################################
## Project: NLKTsoil
## Script purpose:  tidymodels version
## Date: Mon May 18 15:35:07 2020
## Author: Henri C Chung
## Revisions: 
##################################################
#also look at microbial diversity 
#TRY OUT SUPPORT VECTOR MACHINES AND GRADIENT BOOSTING JUST FOR FUN 

## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("tidyverse", "parallel", "parsnip", "rsample", "yardstick", "reshape2", 
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

rm(list = ls()) #clear environment
source("code/nltksoil_functions.R") #Load functions
args = commandArgs(trailingOnly=TRUE)
#set up folder structure
dataFolder <- "data/"; dir.create(dataFolder, showWarnings = FALSE)
curr <- format(Sys.time(), "%Y%m%d/") #Set today's date as output
outputFolder <- paste("output/", curr, sep = ""); dir.create(outputFolder, showWarnings = FALSE)

#Model Parameters
nfolds = 10
nrepeats = 10
PAR = TRUE
embedding <- args[1]
filename <- args[2]
file <- readRDS(paste("data/rds/", filename, ".RDS", sep = ""))    
taxonomy <- readRDS(paste(dataFolder, "rds/taxonomy.RDS", sep = ""))

training = training(file)
testing = testing(file)
message("loaded files")
## Section 
##################################################
if(is.numeric(training$response)){NUM = TRUE}else{NUM = FALSE}
message(paste("NUM is ", NUM))
if(NUM == TRUE){
  #regression model
  tune_spec <- rand_forest(
    mtry = tune(), 
    min_n = tune()
  ) %>%
    set_mode("regression") %>%
    set_engine("randomForest")
  
  #define metric evaluations for regression
  multi_metric <- metric_set(mae, mase, rmse, rsq,ccc)
  
  #custom function for macro_weighted
  error_multi_metric <- function (x) {
    return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), 
                                 na_rm = TRUE, estimator = "macro_weighted"), error=function(e) NULL))
  }
  #
  cv <- vfold_cv(training, v = nfolds, repeats = nrepeats)
  switch(embedding, 
         "word" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_word(all_predictors()) 
         },
         "pca" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_cpca(all_predictors()) 
         },
         "abundance" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_abundance(all_predictors())
         },
         "wordsym" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_wordsym(all_predictors()) 
         })
}else{
  #classification model
  tune_spec <- rand_forest(
    mtry = tune(), 
    #trees = tune(,)tidy
    min_n = tune()
  ) %>%
    set_mode("classification") %>%
    set_engine("randomForest") 
  #define metric evaluations for classification
  multi_metric <- metric_set(accuracy, bal_accuracy, f_meas, kap, npv,
                             ppv, sensitivity, specificity, pr_auc, roc_auc)
  
  #custom function for macro_weighted
  error_multi_metric <- function (x) {
    return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), 
                                 na_rm = TRUE, estimator = "macro_weighted"), error=function(e) NULL))
  }
  
  cv <- vfold_cv(training, v = nfolds, repeats = nrepeats, strata = response)
  
  switch(embedding, 
         "word" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_word(all_predictors()) %>%
             step_mutate(response = droplevels(response))
         },
         "pca" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>%  
             step_cpca(all_predictors()) %>%
             step_mutate(response = droplevels(response))
         },
         "abundance" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_abundance(all_predictors()) %>%
             step_mutate(response = droplevels(response))
         },
         "wordsym" = {
           rec <- recipe(head(training)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")  %>% 
             step_wordsym(all_predictors()) %>%
             step_mutate(response = droplevels(response))
         })
}

message("defined model")


#just for reference
prep <- prep(rec)
juiced <- juice(prep)
#define tuning model workflow
tune_rf <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(tune_spec)
#define tuning grid parameters
rf_grid <- grid_regular(
  mtry(range = c(10, 50)),
  min_n(range = c(2, 10)),
  levels = 5)
message("set tuning parameters")
if(PAR == TRUE){
ncores <- detectCores() - 1
cl <- makeCluster(12)
registerDoParallel()
message(paste("registered", ncores, "cores in parallel"))
}
message("tunnning...")
tune_results <- tune_grid(tune_rf, resamples = cv, grid = rf_grid, metrics = multi_metric, 
                control = control_grid(save_pred = TRUE, 
                                       allow_par = TRUE, 
                                       verbose = TRUE,
                                       extract =  function (x) extract_model(x))) 
message("tuning finished")
if(PAR == TRUE){
registerDoSEQ()
stopCluster(cl)
message("stopped clusters")
}

###finalize model on testing
message("fitting final model")
if(NUM == TRUE){
  param_final <- tune_results %>%
    select_best(metric = "rmse")
  rf_workflow <- tune_rf %>%
    finalize_workflow(param_final)
  
  fit <- rf_workflow %>%
    last_fit(file)
}else{
  param_final <- tune_results %>%
    select_best(metric = "f_meas")
  rf_workflow <- tune_rf %>%
    finalize_workflow(param_final)
  
  fit <- rf_workflow %>%
    last_fit(file)
}
if(!NUM){
tune_results <- tune_results %>%
  mutate(.metrics2 = map(.predictions, custom_evaluate)) %>%
  mutate(.fold = gsub("Repeat|Fold", "", paste(id2, id, sep = "."))) %>%
  mutate(pr_curves = map(.predictions, custom_pr_curve)) %>%
  mutate(roc_curves = map(.predictions, custom_roc_curve)) %>%
  select(-c(splits))
}
message("saving results")
saveRDS(tune_results, paste("tuning_", filename,embedding, ".RDS", sep = ""))
saveRDS(fit, paste(filename,embedding, ".RDS", sep = ""))