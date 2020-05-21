##################################################
## Project: NLKT soil
## Script purpose:  Use embedding space to train random Forest models
## Date: Mon May 18 12:26:55 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
#load packages
list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel", "pROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
source("code/1_microbiome_data_functions.R")
rf_word <- readRDS("data/2_rf_word.RDS")
rf_pca <- readRDS("data/2_rf_pca.RDS")
rf_abundance <- readRDS("data/2_rf_abundance.RDS")
predict <- readRDS("data/2_rf_abundance_predict.RDS")


## Section B: Train Random Forest Models
##################################################
#Define training samples
train_holdout <- 0.80
n = nrow(rf_word) * train_holdout #n = the number of samples * our training holdout
train <- base::sample(nrow(rf_pca),n) #randomly sample n bacteria-samples
#Subset data
rf_word_train <- rf_word[train,]
rf_word_train <- mutate(rf_word_train, feature = factor(feature)) #after subsetting, we have to re-factorized the feature column to account for now missing factor levels. 
rf_word_test <- rf_word[-c(train),]
#
rf_pca_train <- rf_pca[train,]
rf_pca_train <- mutate(rf_pca_train, feature = factor(feature)) 
rf_pca_test <- rf_pca[-c(train),]
#
rf_abundance_train <- rf_abundance[train,]
rf_abundance_test <- rf_abundance[-c(train),]

#train models
word_model <- randomForest::randomForest(feature ~ ., data = rf_word_train)
pca_model <- randomForest::randomForest(feature~., data = rf_pca_train)
#(takes several minutes)
abundance_model <- randomForest::randomForest(y = as.factor(predict$feature[train]), x = rf_abundance_train) 
#check model performance on training data
pred_word_train <- predict(word_model, rf_word_train, type = "class")
pred_pca_train <- predict(pca_model, rf_pca_train, type = "class")
pred_abundance_train <- predict(abundance_model, rf_abundance_train, type = "class")
#Accuracy, or TP+FP/N
mean(pred_word_train == rf_word_train$feature)
mean(pred_pca_train == rf_pca_train$feature)
mean(pred_abundance_train == (predict$feature[train]))
