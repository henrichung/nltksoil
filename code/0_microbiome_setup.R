##################################################
## Project: nltksoil
## Script purpose:  Setup up working samples of soil microbiome data
## Date: Sun May 24 09:44:29 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
list.of.packages <- c("tidyverse", "tidytext", "quanteda", "rsample", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

rm(list = ls()) #clear environment
source("code/nltksoil_functions.R") #Load functions

#set up folder structure
dataFolder <- "data/"; dir.create(dataFolder, showWarnings = FALSE)
curr <- format(Sys.time(), "%Y%m%d/") #Set today's date as output
outputFolder <- paste("output/", curr, sep = ""); dir.create(outputFolder, showWarnings = FALSE)

#Parameters
cdebug = TRUE#If using the debugging dataset instead of full (needs a lot of RAM)
resp = "envo_biome_2" #response variable to use 
#resp = "ph"
query = "soil" #sample information to filter from
#file locations
sample_filename <- paste(dataFolder, "emp_qiime_mapping_release1.tsv" , sep = "")
OTU_filename <- paste(dataFolder, "table.from_biom.txt", sep = "")

## Section B: Read in metadata 
##################################################

#filter to soil samples
sample_data <- read.csv( sample_filename, sep = "\t") %>% 
  dplyr::filter(grepl(query, Description)) %>% #filter based on Description
  mutate(response = !!as.name(resp)) %>% 
  select(c("response", "X.SampleID")) %>%#filter data into description
  mutate(X.SampleID = as.character(X.SampleID))
int_filter <- dplyr::pull(sample_data, "X.SampleID") #pull sample names with soil in description
saveRDS(sample_data, paste(dataFolder, "rds/sample_data.RDS", sep = ""))

if(is.numeric(sample_data$response)){NUM = TRUE}else{NUM = FALSE}

#Read Taxonomy Data
taxon <- read_table2(paste(dataFolder, "97_otu_taxonomy.txt", sep = ""), col_names = F) %>%
  mutate_all(funs(gsub("[a-z]__|;", "", .))) %>%
  mutate_all(list(~na_if(.,""))) 
colnames(taxon) <- c("sample", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
saveRDS(taxon, paste(dataFolder, "rds/taxonomy.RDS", sep = ""))

#file names for EMP database (after python processing)
if(cdebug == FALSE){
  #read in with data.table::fread for speed
  OTU_table <- data.table::fread(OTU_filename) %>%
    as_tibble()
  #reshape taxa rows
  OTU_table <- column_to_rownames(OTU_table,var =  "#OTU ID")
  taxa <- rownames(OTU_table)
  
  }else{
  OTU_table <- readRDS(paste(dataFolder, "rds/OTU_working.RDS", sep = ""))
}

## Section C: Subset OTU dataset by samples of interest 
##################################################
int_filter <- int_filter[int_filter %in% colnames(OTU_table)]

OTU_sub <- select(OTU_table, all_of(int_filter)) #filter to soil samples
OTU_sub <- OTU_sub[!Matrix::rowSums(OTU_sub) == 0 ,] #remove 0 rowSum
OTU_sub <- rownames_to_column(as.data.frame(t(OTU_sub)), "X.SampleID") %>%
  left_join(sample_data, by = "X.SampleID") %>%
  select(-c("X.SampleID")) %>%
  drop_na(response)

#SPLIT 
if(NUM){
  split <-initial_split(OTU_sub, prop = 9/10) 
  setaside <- testing(split); dim(setaside) #10% set aside
  current <- initial_split(training(split), prop = 8/9); dim(current)
}else{
split <-initial_split(OTU_sub, prop = 9/10, strata = "response") 
setaside <- testing(split); dim(setaside) #10% set aside
current <- initial_split(training(split), prop = 8/9, strata = "response"); dim(current)
}

saveRDS(setaside, paste(dataFolder, "rds/setaside", resp, ".RDS", sep = ""))
saveRDS(current, paste(dataFolder, "rds/current", resp, ".RDS", sep = ""))




