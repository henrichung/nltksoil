##################################################
## Project: NLTK Soil
## Script purpose:  take in microbiome
## Date: Thu May 07 16:34:10 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################

#biocmanager packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomformat");library(biomformat)

#normal packages
list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
#custom functions
source("code/1_microbiome_data_functions.R")

#parameters
dataFolder <- "data/"
OTU_filename <- "emp_cr_gg_13_8.subset_2k.biom"
sample_filename <- "emp_qiime_mapping_subset_2k.tsv"
columns_filename <- "out.csv"
train_holdout <- 0.66
###BIG DATA ADJUST
#Alternative filenames for complete dataset
#OTU_filename <- "table.from_biom.txt"
#sample_filename <- "emp_qiime_mapping_qc_filtered.tsv"

#load the data
OTU_file <- paste(dataFolder, OTU_filename, sep = "")
sample_file <- paste(dataFolder, sample_filename, sep = "")
column_file <- paste(dataFolder, columns_filename, sep = "")

#convert database OTU to co-occurence matrix
biom <- biomformat::read_biom(OTU_file)
OTU_table <- biomformat::biom_data(biom)
sample_data <- read_tsv(sample_file)

###BIG DATA ADJUST
#column_names <- readr::read_csv(column_file, col_names = c("temp")) %>%
#  pull(temp)
#reshape data for complete data
#OTU_table <- base::read.table(OTU_file, sep = "\t")
#colnames(OTU_table)[1] <- "temp" #taxa row names is read as first column
#rownames(OTU_table) <- OTU_table$temp #set rownames to that column
#OTU_table <- (OTU_table[,-1]) #remove temp column
#colnames(OTU_table) <- column_names #add column names

int_data <- dplyr::filter(sample_data, grepl('soil', Description)) %>% #filter data into description
  mutate(feature = envo_biome_2)
int_filter <- dplyr::pull(int_data, "#SampleID") #pull sample names with soil in description
###BIG DATA ADJUST
#int_filter <- int_filter[int_filter %in% column_names] #subset names with soil with existing column names in data

OTU_int <- OTU_table[,int_filter]
OTU_sub <- OTU_int[!Matrix::rowSums(OTU_int) == 0 ,]
###BIG DATA ADJUST
#OTU_int <- Reduce(cbind2, parallel:mclapply(OTU_int, Matrix, sparse = TRUE))

###
#train word embedding
cooccur_table <- Dtm2Tcm(t(OTU_sub))
glove = GlobalVectors$new(rank = 50, x_max = 10)
wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.01, n_threads = 8)
wv_context = glove$components
dim(wv_context)
word_vectors = wv_main + t(wv_context)

#PCR embeddings
pcr_embeddings <- prcomp(t(OTU_sub),rank = 50)$rotation

#subset data by taxa data of interest into train and test
predict <- select(int_data, c("#SampleID", "feature"))

rf_word <- add_predict(x = word_vectors, y = OTU_sub, predict = predict)
rf_pca <- add_predict(x = pcr_embeddings, y = OTU_sub, predict = predict)
#relative abundance
rf_abundance <- apply(OTU_sub, MARGIN = 2, FUN = function(x){x/sum(x)}) %>%
  as.matrix() %>%
  t()
#
n = nrow(rf_pca) * train_holdout
train <- base::sample(nrow(rf_pca),n)
###
rf_word_train <- rf_word[train,]
rf_word_train <- mutate(rf_word_train, feature = factor(feature)) 
rf_word_test <- rf_word[-c(train),]
###
rf_pca_train <- rf_pca[train,]
rf_pca_train <- mutate(rf_pca_train, feature = factor(feature)) 
rf_pca_test <- rf_pca[-c(train),]
###
rf_abundance_train <- rf_abundance[train,]
rf_abundance_test <- rf_abundance[-c(train),]

word_model <- randomForest::randomForest(feature ~ ., data = rf_word_train)
pca_model <- randomForest::randomForest(feature~., data = rf_pca_train)
abundance_model <- randomForest::randomForest(y = as.factor(predict$feature[train]), x = rf_abundance_train)

pred_word_train <- predict(word_model, rf_word_train, type = "class")
pred_pca_train <- predict(pca_model, rf_pca_train, type = "class")
pred_abundance_train <- predict(abundance_model, rf_abundance_train, type = "class")

mean(pred_word_train == rf_word_train$feature)
mean(pred_pca_train == rf_pca_train$feature)
mean(pred_abundance_train == (predict$feature[train]))
####
pred_word_test <- predict(word_model, rf_word_test, type = "class")
pred_pca_test <- predict(pca_model, rf_pca_test, type = "class")
pred_abundance_test <- predict(abundance_model, rf_abundance_test, type = "class")

mean(pred_word_test == rf_word_test$feature)
mean(pred_pca_test == rf_pca_test$feature)
mean(pred_abundance_test == predict$feature[-c(train)])

###
resTrain <- table(predTrain, rf_train$feature) 
mean(predTrain == rf_train$feature)
predTest <- predict(model1, rf_test, type = "class")
resTest <- table(predTest, rf_test$feature) 
mean(as.character(predTest) == rf_test$feature)
