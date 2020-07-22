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
total_samples <- ncol(soil_raw)
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
###
#model Parameters
nrepeats <- 10
prevalence <- 0.05
ncores <- 10
PAR = TRUE

#iterations
#responses <- c("principal_investigator", "empo_3", "empo_2", "env_feature", "env_material", "envo_biome_1", "envo_biome_2", "envo_biome_3")
responses = "envo_biome_2"
taxa_space <- c("species", "genus", "family" , "order", "class", "phylum")
model_space <- c("randomforest") #xgboost
rank_levels <- c(50, 100, 150, 200, 250, 300)
embedding_space <- c( "transformword", "transformpca", "abundance", "base", "word", "pca")
#set up parallel processing
if(PAR){
  cl <- makeForkCluster(ncores, outfile ="out.txt")
  registerDoParallel(cl)
  message(Sys.time()," REGISTERED ", ncores, " CORES")
}

#filter data by prevalence
filtered_train <- custom_prev_filter(soil_train, prevalence)
filtered_test <- soil_train[as.numeric(rownames(filtered_train)),]
message(Sys.time()," FILTERED DATA")
## Section B: Reshape data
## Section B: Reshape data
##################################################

for(i in 1:length(taxa_space)){
taxa_level <- taxa_space[i]
custom_summarize_taxa <- function(x, y){
      remove <- paste(substr(y,1,1), "__", sep = "")
        taxa <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
        taxa <- stri_remove_empty(str_remove(taxa, y))

  if(y == "species"){
      res <- x %>%
        separate(taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "; ") %>%
        select(-c(all_of(taxa)))

      res <- res %>%
        mutate_(taxonomy = y) %>%
        mutate(taxonomy = gsub(remove, "", taxonomy)) %>%
        mutate(species = paste("s__", 1:nrow(res), sep = ""))

    }else{

      res <- x %>%
        separate(taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "; ") %>%
        select(-c(all_of(taxa))) %>%
        drop_na(!!y)

      res <- setDT(res)[, lapply(.SD, sum, na.rm = TRUE), by = y]

      res <- res %>%
        mutate_(taxonomy = y) %>%
        mutate(taxonomy = gsub(remove, "", taxonomy)) %>%
        filter(taxonomy != "")
      }
    return(res)
}
summarized_train <- custom_summarize_taxa(filtered_train, taxa_level)
summarized_test <- custom_summarize_taxa(filtered_test, taxa_level)
message(Sys.time()," GROUP TAXA BY ", toupper(taxa_level))
message(Sys.time()," NUMBER OF UNIQUE ", toupper(taxa_level), " IN TRAINING SET : ", length(unique(summarized_train$taxonomy)))
message(Sys.time()," NUMBER OF UNIQUE ", toupper(taxa_level), " IN TESTING SET : ", length(unique(summarized_test$taxonomy)))
#reshape data function; Remove taxonomy, make rownames as SampleID
custom_reshape1 <- function(x){
  rownames(x) = NULL
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
for(n in 1:length(rank_levels)){
rank <- rank_levels[n]
word_matrix <- try(create_word_embedding(reshaped_train[,-1], rank = rank, learning_rate = 0.1), silent = TRUE)
pca_matrix <- create_pca_embedding(reshaped_train[,-1], rank = rank)
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
          levels = 10)
         },
       "xgboost" = {
         #model
         model <- boost_tree(mtry = tune()) %>% set_engine("xgboost") %>% set_mode("classification"); 
         #search tune parameters
         search_grid <- grid_regular(
          mtry(range = c(10, 100)),
          levels = 10)
         }
      )

#recipe
remove <- paste(substr(taxa_level,1,1), "__", sep = "")
base_rec <- recipe(head(model_train)) %>%
  update_role(starts_with(remove), new_role = "predictor") %>%
  update_role(response, new_role = "outcome") %>%
  step_naomit(response)

for(l in 1:length(embedding_space)){
embedding <- embedding_space[l]
#current embedding space
message(Sys.time()," ", toupper(modeltype), ": CURRENTLY ON ", toupper(rank), " ",toupper(embedding))

#define embedding space
switch(embedding, 
        "base"  = {},
        "abundance" = {rec <- step_abundance(base_rec, all_predictors())},
        "pca" = {rec <- step_cpca(base_rec, all_predictors(), options = list("rank" = rank, names = TRUE))},
        "word" = {rec <- step_word(base_rec, all_predictors(), options = list("rank" = rank, names = TRUE))}, 
        "transformword" = {rec <- step_transform(base_rec, all_predictors(), options = list("transform_embedding" = word_matrix, names = TRUE))},
        "transformpca" = {rec <- step_transform(base_rec, all_predictors(), options = list("transform_embedding" = pca_matrix, names = TRUE))}      
        )
if(class(word_matrix) == "try-error" && embedding == "transformword" | ncol(pca_matrix) < rank && embedding == "transformpca"){message("NO EMBEDDING: SKIP"); next}
#drop unused factor levels
rec <- step_mutate(rec, response = droplevels(response))

#tuning work flow
tune_workflow <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(model)

#cross validation for training
folds_train <- vfold_cv(model_train, v = nrepeats, strata = response)

#check tune results
message(Sys.time()," TUNING MODEL")
options(tidymodels.dark = TRUE)
multi_metric <- metric_set(accuracy, bal_accuracy, detection_prevalence, f_meas, j_index, kap, mcc, npv, ppv, precision, recall, sens, spec)
tune_results <- tune_grid(tune_workflow, 
                          resamples = folds_train, 
                          grid = search_grid, 
                          metrics = multi_metric,
                          control = control_grid(
                            verbose = TRUE,
                            allow_par = TRUE,
                            extract = function (x) {xtract_model(x)$prediction.error},
                            save_pred = TRUE))

#DONT FORGET RANKINGS ARE LIMITED TO 50 = WILL RESULT IN MTRY ERROR
message(Sys.time()," SELECTING BEST MODEL")
#what are the best parameters
show_best(tune_results, metric = "accuracy",n = 5)
#pick out the best parameters based on a performance metric
best_parameters <- select_best(tune_results, metric = "accuracy")

## Section D: final model 
##################################################
#add the best parameters to the work flow we described before
message(Sys.time()," FITTING FINAL MODEL")
final_workflow <- finalize_workflow(tune_workflow, best_parameters)
folds_test <- vfold_cv(model_test, v = nrepeats, strata = response)
temp_results <- fit_resamples(final_workflow, resamples = folds_test, 
                                              metrics = multi_metric,
                                              control = control_resamples(
                                                 verbose = TRUE,
                                                 allow_par = TRUE,
                                                 extract = function (x){extract_model(x)$prediction.error},
                                                 save_pred = TRUE))

message(Sys.time()," SAVING RESULTS")
final_results <- list(tune_results, temp_results)
names(final_results) <- c("tuning", "final")
saveRDS(final_results, paste(outputFolder,"/", taxa_level,"_", resp,"_", modeltype, "_", rank, "_", embedding, ".RDS", sep = ""))
    }
   }
  }
 }
}
