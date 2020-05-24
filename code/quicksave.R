list.of.packages <- c("textmineR", "tidyverse", "randomForest", "text2vec", "parallel", "pROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

#clear environment
rm(list = ls())
source("code/nltksoil_functions.R")

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
int_data <- dplyr::filter(sample_data, grepl("soil", pull(sample_data["Description"]))) #filter data into description
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
saveRDS(OTU_sub, "data/OTU_sub.RDS")
