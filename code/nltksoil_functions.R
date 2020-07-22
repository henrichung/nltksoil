custom_nullsample <- function(x){
  train <- analysis(x)$response
  test <- assessment(x)$response
  preds = sample(size = length(test), x = names(table(train)), prob = as.numeric(table(train)), replace = T)
  res <- tibble(estimate = preds, response = test)
  return(res)
}


#prevalence filter function
custom_prev_filter <- function(x,y){
  bin <- x > 0
  index <- Matrix::rowSums(bin)/(ncol(bin)-1)
  index <- index > y 
  res <- x[index,]
  return(res)
}

create_word_embedding <- function(input, sym = FALSE, rank = 50){
  cooccur_table <- cofxn((input), sym)
  glove = GlobalVectors$new(rank = rank, x_max = 1, learning_rate = 0.1) #hyperparameters
  wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.1, n_threads = 1)
  wv_context = glove$components
  word_vectors = wv_main + t(wv_context)
  colnames(word_vectors) <- paste(remove, 1:ncol(word_vectors), sep = "")
  return(word_vectors)
}


create_pca_embedding <- function(x_, rank = 50){
  #temp_name <- paste(substr(colnames(x_)[1],1,1), "__", sep = "")
  x_ <- (as.matrix(x_))
  pca <- prcomp(x_,rank = rank, center=TRUE)$rotation
  colnames(pca) <- paste(remove, 1:ncol(pca), sep = "")
  return(pca)
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
  rownames(res) <- c(categories, "_average")
  res <- rownames_to_column(res, ".variable")
  if(weighted == TRUE){
    wa <- res %>%
      filter(.variable == "_average") %>% 
      mutate(.variable = gsub("_average", "_average_weighted", .variable))
    res <- bind_rows(res, wa)
  }
  res <- reshape2::melt(res, id.vars = ".variable")
  colnames(res) = c(".class", ".metric2", ".value")
  return(as_tibble(res))
}

cofxn <- function(x, sym){
  x <- as.matrix(x)
  x <- as(x, "dgTMatrix")
  y <- x > 0
    if(sym){
      z <- Matrix::crossprod(y)
    }else{
      z <- Matrix::crossprod(x,y)
      #diag(z) <- 0
    }
  return(z)
}

custom_split <- function(x, prop = 1/5){
  n <- nrow(x) * prop
  ind <- sample(1:nrow(x),n)
  res <- list(x[ind,], x[-ind,])
  names(res) <- c("training", "testing")
  return(res)
}
##EMBEDDING RECIPES 
##################################################
#WORD EMBEDDINGS
word_embedding <- function(input, sym = FALSE, rank = 50){
  temp_name <- paste(substr(colnames(input)[1],1,1), "__", sep = "")
  cooccur_table <- cofxn((input), sym)
  glove = GlobalVectors$new(rank = rank, x_max = 1, learning_rate = 0.1) #hyperparameters
  wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.1, n_threads = 1)
  wv_context = glove$components
  word_vectors = wv_main + t(wv_context)
  res <- crossprod(x = word_vectors, y = t(as.matrix(input))) %>%
    t() %>%
    as.data.frame()
  names(res) <- paste(temp_name, 1:ncol(res), sep = "")
  return(res)
}
step_word_new <- function(terms, role, trained, skip, columns, options, id) {
  step(subclass = "word",  
       terms = terms, 
       role = role, 
       trained = trained, 
       skip = skip, 
       columns = columns, 
       options = options, 
       id = id)
}

step_word<-function(
  recipe, 
  ..., 
  role = "predictor", 
  trained = FALSE, 
  skip = FALSE,  
  columns = NULL, 
  options = list(rank = 50, names = TRUE), 
  id = rand_id("word")) {
  
  terms = ellipse_check(...)
  
  add_step(recipe, 
           step_word_new(
             terms = terms,
             role = role, 
             trained = trained,  
             skip = skip, 
             columns = columns,
             options = options,
             id = id
             )
    )
  }

prep.step_word <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  
  if (!any(names(x$options) == "rank")) {
    x$options$rank <- 50
  } 
  
  step_word_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    options = x$options,
    id = x$id
  )
}

bake.step_word <- function(object, new_data, ...) {
  predictors <- word_embedding(input = dplyr::select(new_data, object$columns), rank = object$options$rank)#
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}

#SYMMETRICAL WORDSYM
wordsym_embedding <- function(input, sym = TRUE, rank = 50){
  temp_name <- paste(substr(colnames(input)[1],1,1), "__", sep = "")
  cooccur_table <- cofxn((input), sym)
  glove = GlobalVectors$new(rank = rank, x_max = 1, learning_rate = 0.1) #hyperparameters
  wv_main = glove$fit_transform(cooccur_table, n_iter = 10, convergence_tol = 0.1, n_threads = 1)
  wv_context = glove$components
  word_vectors = wv_main + t(wv_context)
  res <- crossprod(x = word_vectors, y = t(as.matrix(input))) %>%
    t() %>%
    as.data.frame()
  colnames(res) <- paste(temp_name, 1:ncol(res), sep = "")
  return(res)

}

step_wordsym_new <- function(terms, role, trained, skip, columns, options, id) {
  step(subclass = "word",  
       terms = terms, 
       role = role, 
       trained = trained, 
       skip = skip, 
       columns = columns, 
       options = options, 
       id = id)
}

step_wordsym<-function(
  recipe, 
  ..., 
  role = "predictor", 
  trained = FALSE, 
  skip = FALSE,  
  columns = NULL, 
  options = list(rank = 50, names = TRUE), 
  id = rand_id("wordsym")) {
  
  terms = ellipse_check(...)
  
  add_step(recipe, 
           step_wordsym_new(
             terms = terms,
             role = role, 
             trained = trained,  
             skip = skip, 
             columns = columns,
             options = options,
             id = id
             )
    )
  }

prep.step_wordsym <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  
  if (!any(names(x$options) == "rank")) {
    x$options$rank <- 50
  } 
  
  step_wordsym_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    options = x$options,
    id = x$id
  )
}

bake.step_wordsym <- function(object, new_data, ...) {
  predictors <- wordsym_embedding(input = dplyr::select(new_data, object$columns), rank = object$options$rank)#
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}
####################################
##############
#PCR Word Embeddings
pca_embedding <- function(x_, rank = 50){
  temp_name <- paste(substr(colnames(x_)[1],1,1), "__", sep = "")
  x_ <- as.matrix(x_)
  pca <- prcomp(x_,rank = rank, center=TRUE)$rotation
  res <- crossprod(x = pca, y = t(x_)) %>%
    t() %>% 
    as.data.frame()
  colnames(res) <- paste(temp_name, 1:ncol(res), sep = "")
  return(res)

}
step_cpca_new <- function(terms, role, trained, skip, columns,options,  id) {
  step(subclass = "cpca",  
       terms = terms, 
       role = role, 
       trained = trained, 
       skip = skip, 
       columns = columns, 
       options = options, 
       id = id)
}

step_cpca<-function(
  recipe, 
  ..., 
  role = "predictor", 
  trained = FALSE, 
  skip = FALSE, 
  columns = NULL,
  options = list(rank = 50, names = TRUE),
  id = rand_id("cpca")) {
  
  terms = ellipse_check(...)
  
  add_step(recipe, 
           step_cpca_new(
             terms = terms, 
             role = role, 
             trained = trained,  
             skip = skip, 
             columns = columns, 
             options = options, 
             id = id))
}

prep.step_cpca <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  if (!any(names(x$options) == "rank")) {
    x$options$rank <- 50
  } 
  step_cpca_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    options = x$options,
    id = x$id
  )
}

bake.step_cpca <- function(object, new_data, ...) {
  predictors <- pca_embedding(dplyr::select(new_data, object$columns), rank = object$options$rank)
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}
####################
#PCR Word Embeddings
is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
abundance_embedding <- function(x_, rank = 50){
  temp_name <- paste(substr(colnames(x_)[1],1,1), "__", sep = "")
  res <- apply(x_, MARGIN = 1, FUN = function(x){x/sum(x, na.rm = T)}) 
  res[is.nan(res)] <- 0
  res <- as.data.frame(t(res))
  names(res) <- paste(temp_name, 1:ncol(res), sep = "")
  return(res)
}

step_abundance_new <- function(terms, role, trained, skip, columns,options,  id) {
  step(subclass = "abundance",  
       terms = terms, 
       role = role, 
       trained = trained, 
       skip = skip, 
       columns = columns, 
       options = options, 
       id = id)
}

step_abundance<-function(
  recipe, 
  ..., 
  role = "predictor", 
  trained = FALSE, 
  skip = FALSE, 
  columns = NULL,
  options = list(rank = 50, names = TRUE),
  id = rand_id("abundance")) {
  
  terms = ellipse_check(...)
  
  add_step(recipe, 
           step_abundance_new(
             terms = terms, 
             role = role, 
             trained = trained,  
             skip = skip, 
             columns = columns, 
             options = options, 
             id = id))
}

prep.step_abundance <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  if (!any(names(x$options) == "rank")) {
    x$options$rank <- 50
  } 
  step_abundance_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    options = x$options,
    id = x$id
  )
}

bake.step_abundance <- function(object, new_data, ...) {
  predictors <- abundance_embedding(dplyr::select(new_data, object$columns), rank = object$options$rank)
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}
####

custom_transform <- function(input, transform_embedding){
    res <- crossprod(x = transform_embedding, y = t(as.matrix(input))) %>%
    t() %>% 
    as.data.frame()
  return(res)
}

step_transform_new <- function(terms, role, trained, skip, columns,options,  id) {
  step(subclass = "transform",  
       terms = terms, 
       role = role, 
       trained = trained, 
       skip = skip, 
       columns = columns, 
       options = options, 
       id = id)
}
##
step_transform<-function(
  recipe, 
  ..., 
  role = "predictor", 
  trained = FALSE, 
  skip = FALSE, 
  columns = NULL,
  options = list(transform_embedding = 0, names = TRUE),
  id = rand_id("transform")) {
  
  terms = ellipse_check(...)
  
  add_step(recipe, 
           step_transform_new(
             terms = terms, 
             role = role, 
             trained = trained,  
             skip = skip, 
             columns = columns, 
             options = options, 
             id = id))
}
#
prep.step_transform <- function(x, training, info = NULL, ...) {
  col_names <- terms_select(terms = x$terms, info = info)
  if (!any(names(x$options) == "transform_embedding")) {
    message("Error: No transformation matrix provided")
  } 
  step_transform_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    skip = x$skip,
    columns = col_names,
    options = x$options,
    id = x$id
  )
}

bake.step_transform <- function(object, new_data, ...) {
  predictors <- custom_transform(dplyr::select(new_data, object$columns), transform_embedding = object$options$transform_embedding)
  new_data[, object$columns] <- NULL
  bind_cols(new_data, predictors)
}

#####

error_multi_metric <- function (x) {
  return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), na_rm = TRUE), error=function(e) NULL))
}
custom_roc_curve <- function(x){
  #x <- mutate(x, response = droplevels(response))
  levels <- paste(as.character(levels(x$response)), collapse = "|")
  res <- return(tryCatch(roc_curve(x, ... = matches(levels), truth = response), error=function(e) NULL))
  return(res)
}

custom_pr_curve <- function(x){
  #x <- mutate(x, response = droplevels(response))
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

#############