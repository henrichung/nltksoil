##################################################
## Project: nltksoil
## Script purpose:  Test out random forest predictions
## Date: Tue Jun 23 08:43:24 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("tidyverse","tidymodels",  "reshape2", "furrr", "future", "doParallel", 
                      "text2vec", "stringi", "data.table", "janitor")#, "Rcpp", "profvis","Matrix", "themis", "stringi", "data.table", "tictoc")
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


#read in RDS data
otu_raw <- readRDS(paste(dataFolder, otuFilename_train, sep = ""))

#Save OTU and taxonomy data
otu_id <-  otu_raw %>%
  select(c("#OTU ID", "taxonomy")) %>%
  rename(otuid= `#OTU ID`)
saveRDS(otu_id, paste(dataFolder, "emp_deblur_90bp.qc_filtered_otu_taxonomy_list.RDS", sep = "/")) 

#read sample data, subset to soil
sampledata <- read.csv(paste(dataFolder, sampleFilename, sep = ""), sep = "\t") %>%
	janitor::clean_names()
colnames(sampledata) <- paste("metadata__", colnames(sampledata), sep = "")
soil_sampledata <- sampledata  %>% 
  filter(grepl("soil", sampledata[["metadata__description"]]))

#get row indices of soil data
soil_ind <- colnames(otu_raw) %in% soil_sampledata[["metadata__x_sample_id"]]
#subset raw otu data by soil_ind and tranpose to sample = row, otu = columns
soil_raw <- otu_raw[,..soil_ind] %>% t()
#rename columns with OTU ID
colnames(soil_raw) <- otu_raw[["#OTU ID"]] #

#############################
#convert matrix soil_raw to dataframe to add columns
soil_raw <- as.data.frame(soil_raw)
#convert rownames to column
soil_raw[["metadata__x_sample_id"]] <- rownames(soil_raw)


custom_prev_filter <- function(input,prevalence){
  prev = prevalence/100
  input <- as.matrix(input)
  bin <- input > 0
  index <- Matrix::colSums(bin)/(nrow(bin))
  index <- index > prev
  res <- input[,index]
  return(as_tibble(res))
}

soil_raw <- custom_prev_filter(soil_raw, 1)

#convert to data table for faster merge
soil_raw <- as.data.table(soil_raw)
soil_sampledata <- as.data.table(soil_sampledata)
#set SQL style key
setkey(soil_raw, "metadata__x_sample_id")
setkey(soil_sampledata, "metadata__x_sample_id")
#merge
soil_data <- merge(soil_raw, soil_sampledata, all.x = TRUE)

#split into training and testing sets.
soil_split <- initial_split(soil_data, prop = 0.9)
soil_setaside <- testing(soil_split) #set aside data
soil_split <- initial_split(training(soil_split), prop = 0.9)
saveRDS(soil_split, paste(dataFolder, "soil_split.RDS", sep = "/")) 
saveRDS(soil_setaside, paste(dataFolder, "soil_setaside.RDS", sep = "/"))
