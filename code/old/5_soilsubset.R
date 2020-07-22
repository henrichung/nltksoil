##################################################
## Project: nltksoil
## Script purpose:  Test out random forest predictions
## Date: Tue Jun 23 08:43:24 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("tidyverse", "parallel", "parsnip", "rsample", "yardstick", "reshape2", 
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel", 
                      "text2vec", "Rcpp", "profvis","Matrix", "themis", "stringi", "data.table", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))

rm(list = ls()) #clear environment
source("code/nltksoil_functions.R") #Load functions

#filenames
dataFolder <- "../../../../work/idoerg/hchung/"
outputFolder <- "../../../../work/idoerg/hchung/outputs/soil"
otuFilename_test <- "emp_deblur_90bp.subset_2k_w_taxonomy2k.RDS"
otuFilename_train <- "emp_deblur_90bp.qc_filtered.RDS"
sampleFilename <- "emp_qiime_mapping_qc_filtered.tsv"

#read sample data
sampledata <- read.csv(paste(dataFolder, sampleFilename, sep = ""), sep = "\t") %>% 
  filter(grepl("soil", Description))

#read in RDS data
otu_raw <- readRDS(paste(dataFolder, otuFilename_train, sep = ""))%>%
  select(-c("#OTU ID"))
#set aside taxonomy
taxonomy_list <- otu_raw$taxonomy
soil_ind <- (colnames(otu_raw) %in% sampledata$X.SampleID)
soil_raw <- otu_raw[,..soil_ind]
set.seed(123)

soil_raw_tidy <- t(soil_raw)
n <- c(sample(1:(ncol(soil_raw)-1), ncol(soil_raw)/10), ncol(soil_raw))
soil_setaside <- soil_raw[,..n] %>%
  mutate(taxonomy = taxonomy_list)
soil_working <- soil_raw %>%
  select(!contains(colnames(soil_setaside[,1:(ncol(soil_setaside)-1)]))) %>%
  t() #
samplenames <- rownames(soil_working)
soil_split <- soil_working  %>%
  as.data.frame() %>%
  mutate(X.SampleID = samplenames) %>%
  initial_split(prop = 0.9)
soil_train <- training(soil_split)
rownames(soil_train) <- NULL
soil_train <- soil_train %>%
  column_to_rownames("X.SampleID") %>%
  t() %>%
  as.data.frame() %>%
  mutate(taxonomy = taxonomy_list)

soil_test <- testing(soil_split)
rownames(soil_test) <- NULL
soil_test <- soil_test %>%
  column_to_rownames("X.SampleID") %>%
  t() %>% 
  as.data.frame()  %>%
  mutate(taxonomy = taxonomy_list)

rm(otu_raw); gc()

saveRDS(soil_test, paste(dataFolder, "soil_test.RDS", sep = "/"))
saveRDS(soil_train, paste(dataFolder, "soil_train.RDS", sep = "/"))
saveRDS(soil_setaside, paste(dataFolder, "soil_setaside.RDS", sep = "/"))