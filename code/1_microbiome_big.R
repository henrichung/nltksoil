##################################################
## Project: NLKT soil
## Script purpose:  Read in microbiome data and transform into embedding space
## Date: Mon May 18 09:37:10 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries
##################################################
#load packages here
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomformat")

list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel", "pROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
source("code/1_microbiome_data_functions.R")

#parameters
dataFolder <- "data/"

OTU_filename <- "table.from_biom.txt"
sample_filename <- "emp_qiime_mapping_qc_filtered.tsv"
columns_filename <- "out.csv"
big = TRUE
OTU_file <- paste(dataFolder, OTU_filename, sep = "")
sample_file <- paste(dataFolder, sample_filename, sep = "")
##soil sample extract
#read with fread for speed
OTU_table <- data.table::fread(OTU_file)
#reshape taxa rows
colnames(OTU_table)[1] <- "temp" #taxa row names is read as first column
rownames(OTU_table) <- OTU_table$temp #set rownames to that column
OTU_table <- (OTU_table[,-1]) #remove temp column
#filter to soil samples
sample_data <- read.csv(sample_file, sep = "\t")
int_data <- dplyr::filter(sample_data, grepl("soil", pull(sample_data["Description"]))) %>% #filter data into description
int_filter <- dplyr::pull(int_data, "X.SampleID") #pull sample names with soil in description
#########
column_file <- paste(dataFolder, columns_filename, sep = "")
column_names <- readr::read_csv(column_file, col_names = c("temp")) %>%
  pull(temp)
colnames(OTU_table) <- column_names #add column names
int_filter <- int_filter[int_filter %in% column_names] #subset names with soil with existing column names in data

###
OTU_table <- to_sparse(OTU_table)
OTU_int <- OTU_table[,int_filter]
OTU_sub <- OTU_int[!Matrix::rowSums(OTU_int) == 0 ,]
#train word embedding
cooccur_table <- Dtm2Tcm(t(OTU_sub))
glove = GlobalVectors$new(rank = 50, x_max = 10)
wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.01, n_threads = 8)
wv_context = glove$components
dim(wv_context)
word_vectors = wv_main + t(wv_context)
#PCR embeddings
pcr_embeddings <- prcomp(t(OTU_sub),rank = 50)$rotation


soil_data <- dplyr::filter(sample_data, grepl("soil", pull(sample_data["Description"]))) 
#word embeddings
rf_word <- crossprod(x = word_vectors, y = OTU_sub) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "X.SampleID")
#pca
rf_pca <- crossprod(x = pcr_embeddings, y = OTU_sub) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "X.SampleID")
#relative abundance
rf_abundance <- apply(OTU_sub, MARGIN = 2, FUN = function(x){x/sum(x)}) %>%
  as.matrix() %>%
  t()

saveRDS(rf_word, "data/2_rf_word.RDS")
saveRDS(rf_pca, "data/2_rf_pca.RDS")
saveRDS(rf_abundance, "data/2_rf_abundance.RDS")
saveRDS(soil_data, "data/2_soil_data.RDS")