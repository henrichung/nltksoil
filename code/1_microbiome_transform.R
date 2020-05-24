##################################################
## Project: NLKT soil
## Script purpose: Read in microbiome RDS from 0_ file and create property embeddings
## Date: Mon May 18 09:37:10 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries
##################################################
#load packages here
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel", "pROC", "biomformat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
#custom functions
source("code/nltksoil_functions.R")

#parameters
dataFolder <- "data/"
OTU_sub <- readRDS(paste(dataFolder, "OTU_sub.RDS", sep = ""))

#sample file names for EMP database (after python processing)
sample_filename <- "emp_qiime_mapping_qc_filtered.tsv"
sample_file <- paste(dataFolder, sample_filename, sep = "")

#filter to soil samples
sample_data <- read.csv(sample_file, sep = "\t")
int_data <- dplyr::filter(sample_data, grepl("soil", pull(sample_data["Description"]))) #filter data into description
int_filter <- dplyr::pull(int_data, "X.SampleID")

## Section B: Train Word Embeddings
##################################################
cooccur_table <- Dtm2Tcm(t(OTU_sub)) #create cooccurence table
glove = GlobalVectors$new(rank = 50, x_max = 10) #hyperparameters
wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.01, n_threads = 8) #run glove algorithm
wv_context = glove$components
dim(wv_context)
word_vectors = wv_main + t(wv_context)

rf_word <- crossprod(x = word_vectors, y = OTU_sub) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "X.SampleID")

## Section C: PCR Word Embeddings
##################################################
pcr_embeddings <- prcomp(t(OTU_sub),rank = 50)$rotation
rf_pca <- crossprod(x = pcr_embeddings, y = OTU_sub) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "X.SampleID")

## Section D: Relative Abundance 
##################################################
rf_abundance <- apply(OTU_sub, MARGIN = 2, FUN = function(x){x/sum(x)}) %>%
  as.matrix() %>%
  t()

saveRDS(rf_word, "data/rf_word.RDS")
saveRDS(rf_pca, "data/rf_pca.RDS")
saveRDS(rf_abundance, "data/rf_abundance.RDS")
