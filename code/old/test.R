##################################################
## Project: nltksoil
## Script purpose:  Test out random forest predictions
## Date: Tue Jun 23 08:43:24 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################

#load required packages
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
outputFolder <- "../../../../work/idoerg/hchung/outputs"
otuFilename_train <- "emp_deblur_90bp.subset_2k_w_taxonomy2k.RDS"
otuFilename_test <- "emp_deblur_90bp.qc_filtered.RDS"
sampleFilename <- "emp_qiime_mapping_qc_filtered.tsv"

#read in RDS data
otu_train <- readRDS(paste(dataFolder, otuFilename_train, sep = "")) %>%
  select(-c("#OTU ID"))
otu_test_raw <- readRDS(paste(dataFolder, otuFilename_test, sep = ""))%>%
  select(-c("#OTU ID"))
set.seed(123)
n <- sample(1:(ncol(otu_test_raw)-1), 2000)
otu_setaside <- otu_test_raw[,..n]
otu_test <- otu_test_raw %>%
  select(!contains(colnames(otu_train[,1:(ncol(otu_train)-1)]))) %>%
  select(!contains(colnames(otu_setaside[,1:(ncol(otu_setaside)-1)])))

#read sample data
sampledata <- read.csv(paste(dataFolder, sampleFilename, sep = ""), sep = "\t")

#model Parameters
nrepeats <- 10
prevalence <- 0.05
ncores <- 32
PAR = TRUE

#iterations
responses <- c("principal_investigator", "empo_3", "empo_2", "env_feature", "env_material", "envo_biome_1", "envo_biome_2", "envo_biome_3")
taxa_levels <- c("genus" , "family" , "order", "class", "phylum", "species")
model_space <- c("randomforest", "xgboost")
embedding_space <- c("word", "wordsym", "pca", "abundance", "base", "word50", "word100", "word150", "word200")

#set up parallel processing
if(PAR){
  cl <- makeForkCluster(ncores, outfile ="out.txt")
  registerDoParallel(cl)
  message(Sys.time()," REGISTERED ", ncores, " CORES")
}

#filter data by prevalence
filtered_test <- custom_prev_filter(otu_test, prevalence)
filtered_train <- otu_train[otu_train$taxonomy %in% filtered_test$taxonomy]
message(Sys.time()," FILTERED DATA")
## Section B: Reshape data
##################################################
source("code/nltksoil_functions.R")
for(i in 1:length(taxa_levels)){
taxa_level <- taxa_levels[i]
custom_summarize_taxa <- function(x, y){
  remove <- paste(substr(y,1,1), "__", sep = "")
  taxa <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  taxa <- stri_remove_empty(str_remove(taxa, y))
  res <- x %>%
    separate(taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "; ") %>%
    select(-c(all_of(taxa))) %>%
    drop_na(!!y)
  res <- setDT(res)[, lapply(.SD, sum, na.rm = TRUE), by = y]
  #res <- aggregate(as.formula(paste(".~", y, sep = "")), res, sum)
  res <- res %>%
    mutate_(taxonomy = y) %>%
    mutate(taxonomy = gsub(remove, "", taxonomy)) %>%
    filter(taxonomy != "")
}
summarized_train <- custom_summarize_taxa(filtered_train, taxa_level)
summarized_test <- custom_summarize_taxa(filtered_test, taxa_level)
message(Sys.time()," GROUP TAXA BY ", toupper(taxa_level))


#reshape data function; Remove taxonomy, make rownames as SampleID
custom_reshape1 <- function(x){
  res <- x %>%
    select(-c("taxonomy")) %>% 
    column_to_rownames(taxa_level) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("X.SampleID")
  return(res)
}

##section here for change relative abundance
message(Sys.time()," RESHAPING DATA")
reshaped_train <- custom_reshape1(summarized_train)
reshaped_test <- custom_reshape1(summarized_test)

## Section C: Response selection and model tuning 
##################################################

for(j in 1:length(responses)){

resp = responses[j]
message(Sys.time()," CURRENTLY PREDICTING ", toupper(resp))
#reshape sample data
sampledata2 <- sampledata %>%
  dplyr::rename(response = !!(resp)) %>%
  select(c("X.SampleID", "response")) 

custom_add_sampledata <- function(x, min){
  res <- x %>%
    left_join(sampledata2, by = "X.SampleID") %>% #add sampledata by ID
    dplyr::select(-c("X.SampleID")) %>% #remove ID column
    group_by(response) %>% #group by response factors
    filter(n() > min) %>% #filter minimum response 
    ungroup() %>% 
    filter(response != "") %>% #remove blanks
    drop_na(response) %>% #remove NA
    mutate(response = droplevels(response)) #Drop missing factors
  return(res)
}
#add sample data
model_train <- custom_add_sampledata(reshaped_train, min = 0)
model_test <- custom_add_sampledata(reshaped_test, min = 0)

#model
for(k in 1:length(model_space)){
modeltype <- model_space[k]

switch(modeltype, 
       "randomforest" = {
         #model
         model <- rand_forest(mtry = tune()) %>% set_engine("ranger") %>% set_mode("classification"); 
         #search tune parameters
         search_grid <- grid_regular(
          mtry(range = c(10, 100)),
          levels = 5)
         },
       "xgboost" = {
         #model
         model <- boost_tree(mtry = tune()) %>% set_engine("xgboost") %>% set_mode("classification"); 
         #search tune parameters
         search_grid <- grid_regular(
          mtry(range = c(10, 100)),
          levels = 5)
         }
      )

#recipe
remove <- paste(substr(taxa_level,1,1), "__", sep = "")
base_rec <- recipe(head(model_train)) %>%
  update_role(starts_with(remove), new_role = "predictor") %>%
  update_role(response, new_role = "outcome") %>%
  step_naomit(response)
#rec <- rec %>% step_smote(response)
for(l in 1:length(embedding_space)){
source("code/nltksoil_functions.R") #Load functions

embedding <- embedding_space[[l]]

#current embedding space
message(Sys.time()," ", toupper(modeltype), ": CURRENTLY ON ", toupper(embedding))

#define embedding space
switch(embedding, 
        "base"  = {},
        "abundance" = {rec <- step_abundance(base_rec, all_predictors)},
        "pca" = {rec <- step_cpca(base_rec, all_predictors())},
        "word50" = {rec <- step_word(base_rec, all_predictors(), options = list("rank" = 50, names = TRUE))},
        "word100" = {rec <- step_word(base_rec, all_predictors(), options = list("rank" = 100, names = TRUE))},
        "word150" = {rec <- step_word(base_rec, all_predictors(), options = list("rank" = 150, names = TRUE))},
        "word200" = {rec <- step_word(base_rec, all_predictors(), options = list("rank" = 200, names = TRUE))},        
#        "word" = {rec <- step_word(base_rec, all_predictors(), rank = tune())},
        "wordsym" = {rec <- step_wordsym(base_rec, all_predictors())}
        )

#tuning work flow
tune_workflow <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(model)

#cross validation for training
folds_train <- vfold_cv(model_train, v = nrepeats, strata = response)

#check tune results
message(Sys.time()," TUNING MODEL")
tune_results <- tune_grid(tune_workflow, 
                          resamples = folds_train, 
                          grid = search_grid, 
                          metrics = metric_set(accuracy),
                          control = control_grid(verbose = TRUE))

#DONT FORGET RANKINGS ARE LIMITED TO 50 = WILL RESULT IN MTRY ERROR
message(Sys.time()," SELECTING BEST MODEL")
#what are the best parameters
message(show_best(tune_results, metric = "accuracy",n = 5))
#pick out the best parameters based on a performance metric
best_parameters <- select_best(tune_results, metric = "accuracy")

## Section D: final model 
##################################################
#add the best parameters to the work flow we described before
message(Sys.time()," FITTING FINAL MODEL")
final_workflow <- finalize_workflow(tune_workflow, best_parameters)
folds_test <- vfold_cv(model_test, v = nrepeats, strata = response)
temp_results <- fit_resamples(final_workflow, folds_test, metrics = metric_set(accuracy))

message(Sys.time()," SAVING RESULTS")
final_results <- list(tune_results, temp_results)
names(final_results) <- c("tuning", "final")
saveRDS(final_results, paste(outputFolder, taxa_level,"_", resp,"_", modeltype, "_", embedding, ".RDS", sep = ""))
   }
  }
 }
}
