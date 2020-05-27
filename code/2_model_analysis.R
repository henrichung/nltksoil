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
library(tidyverse)
library(parsnip)
library(rsample)
library(yardstick)
library(reshape2)
#clear environment
rm(list = ls())
#custom functions
source("code/nltksoil_functions.R")

#args = commandArgs(trailingOnly=TRUE)
#rdsfile <- args[1]
#resp = args[2]

rdsfile <- "rf_pca"
resp = "envo_biome_2"
title = paste(rdsfile, resp, sep = "_")
dataFolder <- "data/"
outputFolder <- "output/"
nfolds = 10
nsamples = 10
## Section B: Reshape data
##################################################
file <- readRDS(paste(dataFolder, rdsfile, ".RDS", sep = ""))
soil_data <- readRDS(paste(dataFolder, "soil_data.RDS", sep = ""))
soil_data <- filter(soil_data, X.SampleID %in% file$X.SampleID)

rf <- file %>%
  left_join(soil_data, by = "X.SampleID") %>%
  mutate_(response = resp) %>%
  select( matches("[PC|V][0-9]{1,2}|response")) %>%
  mutate(response = gsub(" ", "_", response)) %>%
  mutate(response = gsub("_biome", "", response)) %>%
  mutate(response = as.factor(response)) %>%
  drop_na()
    
#split into setaside, training, testing
split <-initial_split(rf, prop = 9/10) 
setaside <- testing(split); dim(setaside) #10% set aside
current <- initial_split(training(split), prop = 8/9); dim(current)
training <- training(current); dim(training) #80% training
testing <- testing(current); dim(testing) #10% testing

#split into folds and repeat samples     
cv_train <- mutate(vfold_cv(training, v = nfolds, repeats = nsamples),
                   df_ana = map (splits,  analysis),
                   df_ass = map (splits,  assessment))
    
my_rf <- function(training, testing){
  #remove empty classes
  training <- mutate(training, response = droplevels(response))
  testing <- mutate(testing, response = droplevels(response))
      
  #set model
  model <-  rand_forest(mode = "classification") %>%
    set_engine("randomForest") %>%
    fit(response ~ ., data = training) 
      
  #predict on training set
  m_train <- predict(model, training) %>% bind_cols(training) %>%
    mutate(estimate = .pred_class) %>%
    select(c(estimate, response)) #classification
  mp_train <- predict(model, training, type = "prob") %>%
    bind_cols(m_train) #probabilities
      
  #predict on testing set
  m_test <- predict(model,testing) %>% bind_cols(testing)%>%
    mutate(estimate = .pred_class) %>%
    select(c(estimate, response)) #classification
  #probability predictions
  mp_test <- predict(model, testing, type = "prob") %>%
    bind_cols(m_test) #probabilities
  
  res <- list()
  res[[1]] <- mp_train
  res[[2]] <- mp_test
  names(res) <- c("train", "test")
  return(res)
}
a <- as.list(cv_train$df_ana)
b <- as.list(cv_train$df_ass)
temp_names<- expand.grid(unique(names(a)), 1:(length(names(a))/length(unique(names(a))))) %>% 
  unite(name, c(Var1, Var2), sep = ".") %>%
  dplyr::pull(name)
names(a) <- temp_names;names(b) <- temp_names
cv_rf0 <- mapply(training = a, testing = b, my_rf, SIMPLIFY = F)
cv_rf <- unlist(cv_rf0, recursive = F, use.names = T)
###
multi_metric <- metric_set(ppv, npv, sens, spec,accuracy,kap, pr_auc, roc_auc)

error_multi_metric <- function (x) {
  return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), na_rm = TRUE), error=function(e) NULL))
}

stats <- lapply(cv_rf, error_multi_metric) %>%
  bind_rows(.id = ".fold") %>%
  separate(.fold, c(".sample",".fold",".set"), sep = "\\.")
    
custom_stats <- lapply(cv_rf, custom_evaluate) %>%
  bind_rows(.id = ".fold") %>%
  separate(.fold, c(".sample",".fold",".set"), sep = "\\.")
write_csv(stats, paste(outputFolder, title, "_stats.csv", sep = ""))
write_csv(custom_stats, paste(outputFolder, title, "_customstats.csv", sep = ""))
###

roc_curves <- lapply(cv_rf, custom_roc_curve) %>%
  bind_rows(.id = ".fold") %>%
  separate(.fold, c(".sample",".fold",".set"), sep = "\\.")
write_csv(roc_curves,paste(outputFolder, title, "_roc.csv", sep = ""))

pr_curves <- lapply(cv_rf, custom_pr_curve) %>%
  bind_rows(.id = ".fold") %>%
  separate(.fold, c(".sample",".fold",".set"), sep = "\\.")
write_csv(pr_curves,paste(outputFolder, title, "_pr.csv", sep = ""))
    


