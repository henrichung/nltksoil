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
cdebug = FALSE#If using the debugging dataset instead of full (needs a lot of RAM)
resp_space = c("ph", "envo_biome_2", "principal_investigator") #response variable to use 
query = "soil" #sample information to filter from
prevalence = 0.05

#taxonomy
taxon <- read_table2(paste(dataFolder, "97_otu_taxonomy.txt", sep = ""), col_names = F) %>%
  mutate_all(funs(gsub("[a-z]__|;", "", .))) %>%
  mutate_all(list(~na_if(.,""))) 
colnames(taxon) <- c("sample", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
saveRDS(taxon, paste(dataFolder, "rds/taxonomy.RDS", sep = ""))

#file locations
sample_filename <- paste(dataFolder, "emp_qiime_mapping_release1.tsv" , sep = "")
sample_data <- read.csv( sample_filename, sep = "\t") 
OTU_filename <- "../../../../work/idoerg/hchung/table.from_biom.txt"

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
#create filter using sample names
int_filter <- dplyr::filter(sample_data, grepl(query, Description)) %>% 
  dplyr::pull("X.SampleID") #pull sample names with soil in description
int_filter <- int_filter[int_filter %in% colnames(OTU_table)]



OTU_int <- select(OTU_table, all_of(int_filter))  %>% 
	t() %>% as.data.frame()
OTU_notint <- select(OTU_table, !all_of(int_filter)) %>% 
	t() %>% as.data.frame()

#split soil data
OTU_int_split <- initial_split(OTU_int, prop = 9/10)
OTU_int_train <- training(OTU_int_split)
OTU_int_test <- testing(OTU_int_split)

#split not soil data
OTU_notint_split <- initial_split(OTU_notint, prop = 9/10)
OTU_notint_train <- training(OTU_notint_split)
OTU_notint_test <- testing(OTU_notint_split)

training <- list(OTU_int_train, OTU_notint_train)
saveRDS(training, "data/rds/training.RDS")


testing <- list(OTU_int_test, OTU_notint_test)
saveRDS(testing, "data/rds/testing.RDS")