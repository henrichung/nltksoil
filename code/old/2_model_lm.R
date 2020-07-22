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


filename <- 'currentph'
file <- readRDS(paste("data/rds/", filename, ".RDS", sep = ""))    
taxonomy <- readRDS(paste(dataFolder, "rds/taxonomy.RDS", sep = ""))

training = training(file)
testing = testing(file)
rm(file)
message("loaded files")


## Section 
##################################################
if(is.numeric(training$response)){NUM = TRUE}else{NUM = FALSE}
  if(NUM == TRUE){
    #regression model
    model <- linear_reg() %>%
      set_engine("lm") %>%
      set_mode("regression")
    
    #define metric evaluations for regression
    multi_metric <- metric_set(mae, mase, rmse, rsq,ccc)
    
    #custom function for macro_weighted
    error_multi_metric <- function (x) {
      return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), 
                                   na_rm = TRUE, estimator = "macro_weighted"), error=function(e) NULL))
    }
    #
    cv <- vfold_cv(training, v = nfolds, repeats = nrepeats)
    rec <-  recipe(value ~ 1, data = training)
    
}

  #just for reference
prep <- prep(rec)
juiced <- juice(prep)
#define tuning model workflow
rf_workflow <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(model)
  
  
if(PAR == TRUE){
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel()
}
  
results <- fit_resamples(rf_workflow, cv, metrics = multi_metric, 
                           control = control_grid(save_pred = TRUE, 
                                                  allow_par = TRUE, 
                                                  verbose = TRUE,
                                                  extract =  function (x) extract_model(x))) 
  
if(PAR == TRUE){
  registerDoSEQ()
  stopCluster(cl)
}
  

message("saving results")
saveRDS(results, paste(filename,"nullph.RDS", sep = ""))

