
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
    as.data.frame() %>%
    mutate(Freq = as.numeric(Freq))
  sum <- sum(temp$Freq)
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
    temp_res <- as.numeric(c(TP/sum, FP/sum, TN/sum, FN/sum, PPV, NPV, sensitivity, specificity))
    names(temp_res) <- c("TP", "FP", "TN", "FN", "PPV", "NPV", "sensitivity", "specificity")
    res[[i]] <- temp_res
  }
  res <- data.frame(do.call(rbind, res))
  res <- rbind(res, (apply(res, 2, mean, na.rm = T)))
  rownames(res) <- c(categories, "average")
  return(res)
}
#Python stuff
#table = load_table("emp_cr_gg_13_8.qc_filtered.biom") 
#a = table.ids()  
#b = a.astype("U") 
#b.tofile("out.csv", sep="\n", format="%s")   