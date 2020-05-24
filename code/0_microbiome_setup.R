##################################################
## Project: nltksoil
## Script purpose:  Setup up working samples of soil microbiome data
## Date: Sun May 24 09:44:29 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel", "pROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
source("code/nltksoil_functions.R")

#parameters
dataFolder <- "data/"

##Use the following lines with the biom python package to proprocess the .biome file into a tsv for R
#table = load_table("emp_cr_gg_13_8.qc_filtered.biom") 
#a = table.ids()  
#b = a.astype("U") 
#b.tofile("out.csv", sep="\n", format="%s")   

#file names for EMP database (after python processing)
OTU_filename <- "table.from_biom.txt"
sample_filename <- "emp_qiime_mapping_qc_filtered.tsv"
columns_filename <- "out.csv"

OTU_file <- paste(dataFolder, OTU_filename, sep = "")
sample_file <- paste(dataFolder, sample_filename, sep = "")

#read in with data.table::fread for speed
OTU_table <- data.table::fread(OTU_file)

#reshape taxa rows
colnames(OTU_table)[1] <- "temp" #taxa row names is read as first column
rownames(OTU_table) <- OTU_table$temp #set rownames to that column
OTU_table <- (OTU_table[,-1]) #remove temp column

#filter to soil samples
sample_data <- read.csv(sample_file, sep = "\t")
int_data <- dplyr::filter(sample_data, grepl("soil", pull(sample_data["Description"]))) #filter data into description
int_filter <- dplyr::pull(int_data, "X.SampleID") #pull sample names with soil in description
saveRDS(int_data, "data/soil_data.RDS")
#add back column (sample) names
column_file <- paste(dataFolder, columns_filename, sep = "")
column_names <- readr::read_csv(column_file, col_names = c("temp")) %>%
  pull(temp)
colnames(OTU_table) <- column_names #add column names
int_filter <- int_filter[int_filter %in% column_names] #subset names with soil with existing column names in data

#Reshape
OTU_table <- to_sparse(OTU_table) #convert to sparse matrix object
OTU_int <- OTU_table[,int_filter] #filter to soil samples
OTU_sub <- OTU_int[!Matrix::rowSums(OTU_int) == 0 ,] #remove 0 rowSums
saveRDS(OTU_sub, "data/OTU_sub.RDS") #save

#create working subsample for debugging on local machine
set.seed(12345) #set random seed
n <- sample(1:ncol(OTU_sub), 100) #sample 100 items
OTU_working <- OTU_sub[,n] #subset to 100 samples
OTU_working <- OTU_working[!Matrix::rowSums(OTU_working) == 0,] #remove empty rows
saveRDS(OTU_working, "data/OTU_working.RDS")
