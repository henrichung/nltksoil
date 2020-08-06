#R --max-ppsize=500000
## Section A: Set up libraries and load data
##################################################
list.of.packages <- c("tidyverse","tidymodels",  "reshape2", "furrr", "future", "doParallel", 
                      "text2vec", "stringi")#, "Rcpp", "profvis","Matrix", "themis", "stringi", "data.table", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))

rm(list = ls()) #clear environment
args <- commandArgs(trailingOnly = TRUE)
source("code/nltksoil_functions.R") #Load functions
#filenames
dataFolder <- "data"
outputFolder <- "../../../../work/idoerg/hchung/outputs/soil"
sampleFilename <- "emp_qiime_mapping_qc_filtered.tsv"

#read sample data
sampledata <- read.csv(paste(dataFolder, sampleFilename, sep = "/"), sep = "\t") %>% 
  filter(grepl("soil", Description))

#read in RDS data
dataFolder <- "../../../../work/idoerg/hchung/"
otu_id <- readRDS(paste(dataFolder, "emp_deblur_90bp.qc_filtered_otu_taxonomy_list.RDS", sep = "/"))
soil_split <- readRDS(paste(dataFolder, "soil_split.RDS", sep = "/"))
soil_train <- training(soil_split)
soil_test <- testing(soil_split)


#model Parameters
nrepeats <- 10 #define the number of folds to use during cross-validation
PAR = TRUE #We speed up our model runs on each fold by parallizing across multiple cores.
ncores <- 36 #Set the number of cores to use here

#tuning parameters
#set number of parameters to test for each model
rf_search_levels = 20
xgb_search_levels = 30

#iterations
#model variables to iterate over
resp <- "envo_biome_2"
taxa_level <- as.character(args[1])
modeltype <-as.character(args[2])
rank<- as.numeric(args[3]) #c(50, 100, 150, 200, 250, 300)
embedding = as.character(args[4])


#taxa_level = "phylum"
#modeltype = "randomforest"
#rank = 50
#embedding = "base"
message(toupper(paste(resp, taxa_level, modeltype, rank, embedding)))

#set up parallel processing
if(PAR){
  cl <- makeForkCluster(ncores, outfile ="out.txt")
  registerDoParallel(cl)
  message(Sys.time()," REGISTERED ", ncores, " CORES")
}


otu_train <- select(soil_train, -contains("metadata"))
if(embedding == "pca" | embedding == "transformpca" ){message("CREATING PCA MATRIX") ; pca_matrix <- create_pca_embedding(otu_train, rank = rank)}
if(embedding == "word" | embedding == "transformword" ){message("CREATING WORD MATRIX") ; word_matrix <- create_word_embedding(otu_train, rank = rank)}
#detach("package:text2vec", unload=TRUE)

#n <- ncol(soil_train)
#m <- n - 1000
#test <- soil_train[1:300, m:n]

message("FITTING RECIPE")
response = paste("metadata", resp, sep = "__")
base_rec <- recipe(head(soil_train)) %>%
  update_role(contains("metadata"), new_role = "id") %>% #select column that start with the taxa prefix ($__), let those be predictors
  update_role(-contains("metadata"), new_role = "predictor") %>%
  update_role(matches(response), new_role = "outcome") %>% #response variable is our outcome
  step_naomit(all_outcomes(), skip = TRUE) %>%
  step_string2factor(all_outcomes(), skip = TRUE) %>%
  step_factor2string(all_outcomes(), skip = TRUE) %>%
  step_aggregate(all_predictors(), options = list("taxa_level" = taxa_level, "otu_key" = otu_id, names = TRUE))

switch(modeltype, 
       "randomforest" = {
         #model
         model <- rand_forest(
          trees = 1000,
          mtry = tune(),
          min_n = tune()
          ) %>% set_engine("ranger") %>% set_mode("classification")
         
         #search tune parameters
         search_grid <- grid_regular(
          mtry(range = c(10, 100)),
          min_n(range = c(1, 20)),
          levels = rf_search_levels)
         },
       "xgboost" = {
         #model 
         model <- boost_tree(
          trees = 1000, 
          tree_depth = tune(), 
          min_n = tune(), 
          loss_reduction = tune(),                     
          sample_size = tune(), 
          mtry = tune(),                               
          learn_rate = tune()                          
          ) %>% set_engine("xgboost") %>% set_mode("classification"); 
         #search tune parameters
          search_grid <- grid_latin_hypercube(
            tree_depth(),
            min_n(),
            loss_reduction(),
            sample_size = sample_prop(),
            finalize(mtry(), soil_train),
            learn_rate(),
            size = xgb_search_levels)
        }
      )


#define embedding space
switch(embedding, 
        "base"  = {rec <- base_rec},
        "abundance" = {rec <- step_abundance(base_rec, all_predictors())},
        "pca" = {rec <- step_cpca(base_rec, all_predictors(), options = list("rank" = rank, names = TRUE))},
        "word" = {rec <- step_word(base_rec, all_predictors(), options = list("rank" = rank, names = TRUE))}, 
        "transformword" = {rec <- step_transform(base_rec, all_predictors(), options = list("transform_embedding" = word_matrix, names = TRUE))},
        "transformpca" = {rec <- step_transform(base_rec, all_predictors(), options = list("transform_embedding" = pca_matrix, names = TRUE))}      
       )

#Define the tuning work flow
tune_workflow <- 
  workflow() %>% #workflow objects
  add_recipe(rec) %>% #add recipe
  add_model(model) #add model

#cross validation for training
#we stratify samples across the sample response to make sure that each fold has an even distribution of classes.
train_folds <- vfold_cv(soil_train, v = nrepeats, strata = response)

#check tune results
options(tidymodels.dark = TRUE)
#Define a vector of metric parameters to evaluate model accuracy
multi_metric <- metric_set(accuracy, bal_accuracy, f_meas, kap, npv, ppv, precision, recall, sens, pr_auc, roc_auc)
message("TUNING MODEL")
tune_results <- tune_grid(tune_workflow,  #tune a parameter grid using the tune_workflow
                          resamples = train_folds,  #apply the workflow to folds created by vfold_cv
                          grid = search_grid,  #search of provided list of parameters
                          #metrics = multi_metric, #calculate performance using metrics in metric set
                          control = control_grid( #control parameters to tuning process
                            verbose = TRUE, #print progress of tuning to console
                            allow_par = TRUE, #allow parallel processing if parallel backend registers
                            #extract = function (x) {extract_model(x)}, #extract model specificiations for each tune iterations
                            save_pred = TRUE)) #save model predictions 

#what are the best parameters
best <- show_best(tune_results, metric = "accuracy", n = 5) #pick top n results from tuning based on metric
best_parameters <- select_best(tune_results, metric = "accuracy")

final_workflow <- tune_workflow %>%
  finalize_workflow(best_parameters) 


#train_results <- fit_resamples(final_workflow, train_folds, control = control_grid( #control parameters to tuning process
#                            verbose = TRUE, #print progress of tuning to console
#                            allow_par = FALSE, #allow parallel processing if parallel backend registers
#                            extract = function (x) {extract_model(x)}, #extract model specificiations for each tune iterations
#                            save_pred = TRUE))

#final_fit <- final_workflow %>%
#  last_fit(soil_split)
message("FITTING FINAL MODEL")
final_model <- parsnip::fit(final_workflow, training(soil_split))


#create custom model prediction function based off proportion of classes within data (equivalent to guessing)
custom_nullsample <- function(x, response = "blank"){ #take in a splits object which defines the data splits created by vfold_cv
  train <- analysis(x)[[response]] #take out class categories from training data
  test <- assessment(x)[[response]] #take out class categories from testing data
  preds = sample(size = length(test), x = names(table(train)), prob = as.numeric(table(train)), replace = T) #make n draws, where n is the length of responses in the test set, from the classes in the training set with equal probability to their proportions in the testing set with replacement
  res <- tibble(estimate = preds, response = test) #merge test responses with test predictions
  return(res) #return results
}

tune_results <- tune_results %>% 
  mutate(.null = future_map(splits, custom_nullsample, response = response)) %>% 
  mutate(.predictions = future_map(.predictions, ~rename_(., response = response))) %>%#create null guess
  select(-c("splits")) #remove splits object which takes a majority of the data within the tune_results objects

message("FITTING TEST")
soil_bootstraps <- bootstraps(soil_test,times = nrepeats, strata = response)

test_results <- soil_bootstraps %>%
  mutate(.probabilities = future_map(splits, ~predict(final_model, new_data =  analysis(.), type = "prob"))) %>%
  mutate(.predictions = future_map(splits, ~predict(final_model, analysis(.)))) %>%
  mutate(.predictions = future_map2(splits, .predictions, ~cbind(analysis(.x)[[response]], .y))) %>%
  mutate(.predictions = future_map(.predictions, ~rename(., response = "analysis(.x)[[response]]"))) %>%
  mutate(.predictions = future_map2(.predictions, .probabilities, cbind)) %>%
  mutate(.predictions = future_map(.predictions, ~mutate(., response = droplevels(response)))) %>%
  select(-c(.probabilities))

train = soil_train[[response]]

custom_test_null <- function(x, y, response){
  test <- analysis(x)[[response]]
  preds <- sample(size = length(test), x = names(table(y)), prob = as.numeric(table(y)), replace = T)
  res <- tibble(estimate = preds, response = test)
}


test_results <- test_results %>% 
  mutate(.null = future_map(splits, custom_test_null, y = train, response = response)) %>% #create null guess
  select(-c("splits")) 

message("SAVING RESULTS")
final_results <- list(tune_results, best, test_results, final_model)
names(final_results) <- c("tune", "best", "test", "model")
#save final results as object along with model parameters as name.
saveRDS(final_results, paste(outputFolder,"/", taxa_level,"_", resp,"_", modeltype, "_", rank, "_", embedding, ".RDS", sep = ""))
