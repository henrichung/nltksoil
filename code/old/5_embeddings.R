## Section A: Set up libraries and load data
list.of.packages <- c("tidyverse", "parallel", "parsnip", "rsample", "yardstick", "reshape2", 
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel", 
                      "text2vec", "Rcpp", "profvis","Matrix", "themis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))


rm(list = ls()) #clear environment
source("code/nltksoil_functions.R") #Load functions
PAR = TRUE
ncores = 36
#set up folder structure
dataFolder <- "data/"; dir.create(dataFolder, showWarnings = FALSE)
curr <- format(Sys.time(), "%Y%m%d/") #Set today's date as output
outputFolder <- paste("output/", curr, sep = ""); dir.create(outputFolder, showWarnings = FALSE)
##read samplename
sample_filename <- paste(dataFolder, "emp_qiime_mapping_release1.tsv" , sep = "")
sample_data <- read.csv( sample_filename, sep = "\t") 
#read training/testing
current <- readRDS("data/rds/current.RDS")

if(PAR){
  cl <- makeForkCluster(ncores, outfile = "out.txt")
  registerDoParallel(cl)
  message(Sys.time()," REGISTERED ", ncores, " CORES")
}


training <- current$notint
cooccur_table <- cofxn(training, FALSE)
message("COOCCURENCE CALCULATED")
glove = GlobalVectors$new(rank = 50, x_max = 10, learning_rate = 0.1) #hyperparameters
message("ALGORITHM START")
wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.1, n_threads = 36) #run glove algorithm
message("EXTRACT COMPONENTS")
wv_context = glove$components
message("GET CONTEXT COMPONENTS")
word_vectors = wv_main + t(wv_context)
message("SAVING WORD VECTORS")
saveRDS(word_vectors, "word_vectors.RDS")

if(PAR){
  registerDoSEQ()
  stopCluster(cl)
  message(Sys.time()," ", modeltype, ": DEREGISTERED CORES")
  message(Sys.time()," ", modeltype, ": ", pryr::mem_used())
}
