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

rf <- rf_word %>%
  left_join(soil_data, by = "X.SampleID") %>%
  mutate(response = envo_biome_2) %>%
  select( matches("V[0-9]{2}|response")) %>%
  mutate(response = gsub(" ", "_", response)) %>%
  mutate(response = gsub("_biome", "", response)) %>%
  mutate(response = as.factor(response))

#split into setaside, training, testing
split <-initial_split(rf, prop = 9/10) 
setaside <- testing(split); dim(setaside) #10% set aside
current <- initial_split(training(split), prop = 8/9); dim(current)
training <- training(current); dim(training) #80% training
testing <- testing(current); dim(testing) #10% testing


cv_train <- mutate(vfold_cv(training, v = 3, repeats = 10),
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

roc_curves %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = .fold, type = .sample)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.set~.level)

pr_curves <- lapply(cv_rf, custom_pr_curve) %>%
  bind_rows(.id = ".fold") %>%
  separate(.fold, c(".sample",".fold",".set"), sep = "\\.")
pr_curves %>%
  ggplot(aes(x = recall, y = precision, color = .fold, type = .sample)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(.set~.level)

