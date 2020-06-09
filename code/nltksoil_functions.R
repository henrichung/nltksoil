
dtm2tcm <- function (dtm) 
{
  dtm_binary <- dtm > 0
  result <- Matrix::t(dtm_binary) %*% dtm
  result
}

custom_evaluate <- function(x, weighted = FALSE){
  res <- list()
  if(sum(colnames(x) %in% "estimate") == 0){x <- rename(x, estimate =  .pred_class)}
  temp <- x %>%
    select(c(estimate, response)) %>%
    table() %>%
    as.data.frame() %>%
    mutate(Freq = as.numeric(Freq))
  colnames(temp) <- c("predict", "observed", "Freq")
  temp <- temp %>%
    mutate(predict = as.character(predict), observed = as.character(observed))
  categories <- unique(temp$predict)
  for(i in 1:length(categories)){
    subject = categories[i]
    TP <- filter(temp, observed == subject & predict == subject) %>% 
      summarize(n = sum(Freq))
    FP <- filter(temp, observed != subject & predict == subject) %>% 
      summarize(n = sum(Freq))
    TN <- filter(temp, observed != subject & predict != subject) %>% 
      summarize(n = sum(Freq))
    FN <- filter(temp, observed == subject & predict != subject) %>% 
      summarize(n = sum(Freq))
    PPV <- TP/(TP+FP) #same as precision
    NPV <- TN/(FN+TN)
    sensitivity <- TP/(TP+FN) #recall
    specificity <- TN/(FP+TN)
    accuracy <- (TP+TN)/(TP+FN+FP+TN)
    balanced_accuracy <- (sensitivity+specificity)/2
    f_meas <- (2*PPV*sensitivity)/(PPV+sensitivity)
    if(weighted == TRUE){fraction <- (filter(temp, observed == subject) %>% nrow()) / (nrow(temp))}else{fraction = 1}
    temp_res_1 <- as.numeric(c(TP, FP, TN, FN)) 
    temp_res_2 <- as.numeric(c(PPV, NPV, sensitivity, specificity, accuracy, balanced_accuracy, f_meas)) * fraction
    temp_res <- c(temp_res_1, temp_res_2)
    names(temp_res) <- c("TP", "FP", "TN", "FN", "PPV", "NPV", "sensitivity", "specificity",  "accuracy", "balance_accuracy", "f_meas")
    res[[i]] <- temp_res
  }
  res <- data.frame(do.call(rbind, res))
  if(weighted == TRUE){res <- rbind(res, (apply(res, 2, sum, na.rm = T)))}else{res <- rbind(res, (apply(res, 2, mean, na.rm = T)))}
  #res <- rbind(res, (apply(res, 2, sum, na.rm = T)))
  rownames(res) <- c(categories, "average")
  res <- rownames_to_column(res, ".variable")
  if(weighted == TRUE){
    wa <- res %>%
      filter(.variable == "average") %>% 
      mutate(.variable = gsub("average", "average_weighted", .variable))
    res <- bind_rows(res, wa)
    }
  return(res)
}

##EMBEDDING RECIPES 
##################################################

word_embedding <- function(x, sym = FALSE, rank = 50, x_max = 10){
  library(text2vec)
  if(sym){cooccur_table <- dtm2tcm(x)}else{cooccur_table <- t(as.matrix(x)) %*% as.matrix(x)}
  glove = GlobalVectors$new(rank = rank, x_max = x_max) #hyperparameters
  wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.01, n_threads = 8) #run glove algorithm
  wv_context = glove$components
  word_vectors = wv_main + t(wv_context)
  res <- crossprod(x = word_vectors, y = t(x)) %>%
    t() %>%
    as.data.frame()
  return(res)
}

step_word_new <- function(terms, role, trained, skip, columns, id) {
  step(subclass = "word",  terms = terms, role = role, 
       trained = trained, skip = skip, columns = columns, id = id)
}

step_word<-function(recipe, ..., role = "predictor", trained = FALSE, skip = FALSE,  
                    columns = NULL, id = rand_id("word")) {
  terms = ellipse_check(...)
  add_step(recipe, 
           step_word_new(terms = terms, role = role, trained = trained,  
                         skip = skip, columns = columns, id = id))
}

prep.step_word <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  step_word_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    id = x$id
  )
}

bake.step_word <- function(object, new_data, ...) {
  predictors <- word_embedding(dplyr::select(new_data, object$columns))
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}
##############
#PCR Word Embeddings
pca_embedding <- function(x_){
  pca <- prcomp(x_,rank = 50)$rotation
  res <- crossprod(x = pca, y = t(x_)) %>%
    t() %>% 
    as.data.frame()# %>%
  #tibble::rownames_to_column(var = "X.SampleID")
  return(res)
}

step_cpca_new <- function(terms, role, trained, skip, columns, id) {
  step(subclass = "cpca",  terms = terms, role = role, 
       trained = trained, skip = skip, columns = columns, id = id)
}

step_cpca<-function(recipe, ..., role = "predictor", trained = FALSE, skip = FALSE,  
                    columns = NULL, id = rand_id("cpca")) {
  terms = ellipse_check(...)
  add_step(recipe, 
           step_cpca_new(terms = terms, role = role, trained = trained,  
                         skip = skip, columns = columns, id = id))
}

prep.step_cpca <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  step_cpca_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    id = x$id
  )
}

bake.step_cpca <- function(object, new_data, ...) {
  predictors <- pca_embedding(dplyr::select(new_data, object$columns))
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}
####################

##Relative Abundance 
abundance_embedding <- function(x_, taxonomy = readRDS(paste(dataFolder, "rds/taxonomy.RDS", sep = "")), taxa = "Phylum"){
  res <- apply(x_, MARGIN = 1, FUN = function(x){x/sum(x)}) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("X.SampleID") %>%
    reshape2::melt(id.vars = "X.SampleID") %>%
    rename(sample = variable) %>%
    left_join(taxonomy, by = "sample") %>%
    mutate(taxa  = gsub("\\[|\\]", "", !!as.name(taxa))) %>%
    group_by(X.SampleID, !!as.name(taxa)) %>%
    summarise(n = sum(value, na.rm = T)) %>%
    ungroup() %>%
    pivot_wider(names_from = !!as.name(taxa), values_from = n) %>%
    drop_na() %>%
    column_to_rownames("X.SampleID")
  return(res)
}

step_abundance_new <- function(terms, role, trained, skip, columns, id) {
  step(subclass = "abundance",  terms = terms, role = role, 
       trained = trained, skip = skip, columns = columns, id = id)
}

step_abundance<-function(recipe, ..., role = "predictor", trained = FALSE, skip = FALSE,  
                   columns = NULL, id = rand_id("abundance")) {
  terms = ellipse_check(...)
  add_step(recipe, 
           step_abundance_new(terms = terms, role = role, trained = trained,  
                         skip = skip, columns = columns, id = id))
}

prep.step_abundance <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  step_abundance_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    id = x$id
  )
}

bake.step_abundance <- function(object, new_data, ...) {
  predictors <- abundance_embedding(dplyr::select(new_data, object$columns))
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}



#######################

rf_predict_cat <- function(x){
  estimate <- predict(fit_workflow, analysis(x)) %>% pull()
  response <- analysis(x)[["response"]]
  probs <- predict(fit_workflow, analysis(x), type = "prob")
  res <- tibble(estimate = estimate, response = response)
  res <- cbind(res, probs)
  return(res)
}

rf_predict_num <- function(x){
  estimate <- predict(fit_workflow, analysis(x)) %>% pull()
  response <- analysis(x)[["response"]]
  res <- tibble(estimate = estimate, response = response)
  return(res)
}

error_multi_metric <- function (x) {
  return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), na_rm = TRUE), error=function(e) NULL))
}
custom_roc_curve <- function(x){
  x <- mutate(x, response = droplevels(response))
  levels <- paste(as.character(levels(x$response)), collapse = "|")
  res <- return(tryCatch(roc_curve(x, ... = matches(levels), truth = response), error=function(e) NULL))
  return(res)
}

custom_pr_curve <- function(x){
  x <- mutate(x, response = droplevels(response))
  levels <- paste(as.character(levels(x$response)), collapse = "|")
  res <- return(tryCatch(pr_curve(x, ... = matches(levels), truth = response), error=function(e) NULL))
  return(res)
}
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


refactor <- function(a,x,y){
  both <- union(levels(a[[x]]), levels(a[[y]]))
  a[[x]] <- factor(a[[x]], levels=both)
  a[[y]] <- factor(a[[y]], levels=both)
  return(a)
}


#Python stuff
#table = load_table("emp_cr_gg_13_8.qc_filtered.biom") 
#a = table.ids()  
#b = a.astype("U") 
#b.tofile("out.csv", sep="\n", format="%s")   


