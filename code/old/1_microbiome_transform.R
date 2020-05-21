##################################################
## Project: NLKT soil
## Script purpose:  Read in microbiome data and transform into embedding space
## Date: Mon May 18 09:37:10 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries
##################################################
#load packages
list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel", "pROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
source("code/1_microbiome_data_functions.R")

#parameters
dataFolder <- "data/"
OTU_filename <- "emp_cr_gg_13_8.subset_2k.biom"
sample_filename <- "emp_qiime_mapping_subset_2k.tsv"

#OTU_filename <- "table.from_biom.txt"
#sample_filename <- "emp_qiime_mapping_qc_filtered.tsv"
#columns_filename <- "out.csv"
big = TRUE
OTU_file <- paste(dataFolder, OTU_filename, sep = "")
sample_file <- paste(dataFolder, sample_filename, sep = "")
##soil sample extract

if(big == FALSE){
  OTU_sub <- custom_import_small(x = OTU_file, y = sample_file, query = "soil", column = "Description")
}else{
  OTU_sub <- custom_import_big(x = OTU_file, y = sample_file, z = columns_filname, query = "soil", column = "Description")
}
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
saveRDS(rf_word, "data/2_rf_word.RDS")
saveRDS(rf_pca, "data/2_rf_pca.RDS")
saveRDS(rf_abundance, "data/2_rf_abundance.RDS")
saveRDS( as.factor(predict$feature), "data/2_rf_abundance_predict.RDS")