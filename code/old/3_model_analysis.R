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
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel", 
                      "text2vec", "Rcpp", "profvis","Matrix", "themis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))

rm(list = ls()) #clear environment
source("code/nltksoil_functions.R") #Load functions
#invisible(sourceCpp("code/matrix.cpp")) #

args = commandArgs(trailingOnly=TRUE)

#set up folder structure
dataFolder <- "data/"; dir.create(dataFolder, showWarnings = FALSE)
curr <- format(Sys.time(), "%Y%m%d/furrrtest/") #Set today's date as output

#Model Parameters
nfolds = 2
nrepeats = 2
PAR = FALSE
ncores <- 2
filename <- "soilenvo_biome_2"#as.character(args[1])
SMOTE = FALSE# as.logical(args[2])

if(SMOTE){outputFolder <- paste("outputs/", curr, "smote/", sep = ""); dir.create(outputFolder, showWarnings = FALSE)}else{
  outputFolder <- paste("outputs/", curr, sep = ""); dir.create(outputFolder, showWarnings = FALSE)
}



model_space <- c("randomforest", "xgboost","svm")
embedding_space <- c("word", "wordsym", "pca", "abundance")

message(Sys.time()," ", toupper(filename))
file <- readRDS(paste("data/rds/", filename, ".RDS", sep = ""))    
training = training(file)
testing = testing(file)
rm(file)

#is response numerical?
if(is.numeric(training$response)){NUM = TRUE}else{NUM = FALSE}

#taxonomy <- readRDS(paste(dataFolder, "rds/taxonomy.RDS", sep = ""))
message(Sys.time(), " LOADED FILES")
message(Sys.time(), " SMOTE IS ", SMOTE)
#Parallelize workflow
if(PAR){
  cl <- makeForkCluster(ncores, outfile = paste(outputFolder,filename,".txt", sep = ""))
  registerDoParallel(cl)
  message(Sys.time()," REGISTERED ", ncores, " CORES")
}

#for(j in 1:length(model_space)){
  j = 1
  modeltype <- model_space[j]

## Section 
##################################################
#Define Model
switch(modeltype,
  "randomforest" = {model <- rand_forest() %>% set_engine("randomForest")},
  "xgboost" = {model <- boost_tree() %>% set_engine("xgboost")},
  "linear" = {model <- linear_reg() %>% set_engine("lm")},
  "svm" = {model <- svm_poly() %>% set_engine("kernlab")}
  )

#Response specific settings
#define mode
#set multimetric function
if(NUM){
  #regression model
    model <- set_mode(model, "regression")
    #define metric evaluations for regression
    multi_metric <- metric_set(mae, mase, rmse, rsq,ccc)
    #custom function for macro_weighted
    cv <- vfold_cv(training, v = nfolds, repeats = nrepeats)
}else{
    #classification model
    model <- set_mode(model, "classification")
    #define metric evaluations for classification
    multi_metric <- metric_set(accuracy, bal_accuracy, f_meas, kap, npv, ppv, sensitivity, specificity, pr_auc, roc_auc)  
    #custom function for macro_weighted
    cv <- vfold_cv(training, v = nfolds, repeats = nrepeats, strata = response)
}
#
#Define error multi_metric to avoid ppv error
error_multi_metric <- function (x){
  return(tryCatch(
    multi_metric(x, 
      truth = response, 
      estimate = estimate, 
      ... = matches(".pred.*"), 
      na_rm = TRUE, 
      estimator = "macro_weighted"), 
  error=function(e) NULL))}

#loop over embedding spaces
#for(i in 1:length(embedding_space)){
  i = 3
  #choose embedding
  embedding <- embedding_space[[i]]
  #current embedding space
  message(Sys.time()," ", modeltype, ": CURRENTLY ON ", toupper(embedding))
  #define recipe
  rec <- recipe(head(training, n = 20)) %>%
             update_role(everything()) %>%
             update_role(response, new_role = "outcome")

  #define embedding space
  switch(embedding, 
         "word" = {rec <- step_word(rec, all_predictors())},
         "pca" = {rec <- step_cpca(rec, all_predictors())},
         "abundance" = {rec <- step_abundance(rec, all_predictors())},
         "wordsym" = {rec <- step_wordsym(rec, all_predictors())}
         )
  #drop levels if categorical
  if(!NUM){
    rec <- step_mutate(rec, response = droplevels(response))
  }
  if(SMOTE){
    rec <- step_smote(rec, response)
  }

  #special case for null linear model
  if(model == "linear" && NUM){
     rec <-  recipe(value ~ 1, data = head(training))
  }
  #define tuning model workflow
  rf_workflow <- 
    workflow() %>%
    add_recipe(rec) %>%
    add_model(model)

  message(Sys.time()," ", modeltype, ": FITTING RESAMPLES...")
  #fit resamples


  future::plan(multicore)
results <- cv %>%
  mutate(.ana = future_map(splits, analysis),
       .ass = future_map(splits, assessment),
       .rec = future_map(.ana, ~prep(rec, training = .x)),
       .ana = future_map(.rec, juice),
       .ass = future_map2(.rec, .ass, ~bake(.x, new_data = .y)), #check!
       .fit = future_map(.ana, ~fit(model, response ~. , data = .x)),
       .fit = future_map(.fit, ~stripRF(.fit)),
       .pred = future_map2(.fit, .ass, ~predict(.x$finalModel, new_data = .y)),
       .prob = future_map2(.fit, .ass, ~predict(.x$finalModel, new_data = .y, type = "prob")))
  select(-c(fit))
  rename(.pred = response)
  #results <- fit_resamples(rf_workflow, cv, metrics = multi_metric, 
  #  control = control_resamples(
  #    save_pred = TRUE, 
  #    allow_par = TRUE,
  #    verbose = TRUE))#,
      #extract =  function (x) extract_model(x))) 
  
  #additional stats for categorical responses
  if(NUM){ results <- results %>%
      message(Sys.time()," ", modeltype, ": CALCULATING BONUS STATS...")
    results <- results %>%
      mutate(.metrics2 = future_map(.predictions, custom_evaluate)) %>%
      mutate(.fold = gsub("Repeat|Fold", "", paste(id2, id, sep = "."))) %>%
      mutate(pr_curves = future_map(.predictions, custom_pr_curve)) %>%
      mutate(roc_curves = future_map(.predictions, custom_roc_curve))
  }
  #save results
  message(Sys.time()," ", modeltype, ": SAVING RESULTS")
  saveRDS(results, file =  paste(outputFolder, "/", filename,"_", modeltype, "_",embedding, ".RDS", sep = ""))
} #end of embeddings loop
} #end of model space loop
#Dregister Parallel
