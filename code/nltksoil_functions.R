
to_sparse <- function(d_table){
  
  i_list <- lapply(d_table, function(x) which(x != 0))
  counts <- unlist(lapply(i_list, length), use.names = F)
  
  sparseMatrix(
    i = unlist(i_list, use.names = F),
    j = rep(1:ncol(d_table), counts),
    x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F),
    dims = dim(d_table),
    dimnames = list(NULL, names(d_table)))
}

custom_evaluate <- function(x){
  res <- list()
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
    PPV <- TP/(TP+FP)
    NPV <- TN/(FN+TN)
    sensitivity <- TP/(TP+FN)
    specificity <- TN/(FP+TN)
    detection_rate <- TP/(TP+FN+FP+TN)
    detection_prevalence <- (TP+TN)/(TP+FN+FP+TN)
    balanced_accuracy <- (sensitivity+specificity)/2
    temp_res <- as.numeric(c(TP, FP, TN, FN, PPV, NPV, sensitivity, specificity, detection_rate, detection_prevalence, balanced_accuracy))
    names(temp_res) <- c("TP", "FP", "TN", "FN", "PPV", "NPV", "sensitivity", "specificity", "detection_rate", "detection_prevalence", "balance_accuracy")
    res[[i]] <- temp_res
  }
  res <- data.frame(do.call(rbind, res))
  res <- rbind(res, (apply(res, 2, mean, na.rm = T)))
  rownames(res) <- c(categories, "average")
  res <- rownames_to_column(res, ".variable")
  return(res)
}

error_multi_metric <- function (x) {
  return(tryCatch(multi_metric(x, truth = response, estimate = estimate, ... = matches(".pred.*"), na_rm = TRUE), error=function(e) NULL))
}
custom_roc_curve <- function(x){
  levels <- paste(as.character(unique(x$response)), collapse = "|")
  res <- return(tryCatch(roc_curve(x, ... = matches(levels), truth = response), error=function(e) NULL))
  return(res)
}

custom_pr_curve <- function(x){
  levels <- paste(as.character(unique(x$response)), collapse = "|")
  res <- return(tryCatch(pr_curve(x, ... = matches(levels), truth = response), error=function(e) NULL))
  return(res)
}
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

CI_95 <- function(x) {
  ymin = min(x)
  lower = lower_ci(mean = mean(x), se = (sd(x)/ sqrt(length(x))),n =  length(x))
  middle = mean(x)
  upper = upper_ci(mean = mean(x), se = (sd(x)/ sqrt(length(x))),n =  length(x))
  ymax = max(x)
  r <- c(ymin, lower, middle, upper, ymax)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}  

CI2_95 <- function(x) {
  ymin = mean(x)
  lower = lower_ci(mean = mean(x), se = (sd(x)/ sqrt(length(x))),n =  length(x))
  middle = mean(x)
  upper = upper_ci(mean = mean(x), se = (sd(x)/ sqrt(length(x))),n =  length(x))
  ymax = mean(x)
  r <- c(ymin, lower, middle, upper, ymax)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}  


lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

#Python stuff
#table = load_table("emp_cr_gg_13_8.qc_filtered.biom") 
#a = table.ids()  
#b = a.astype("U") 
#b.tofile("out.csv", sep="\n", format="%s")   