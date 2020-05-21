#######
#TODO
#REPLACE RANDOM FOREST WITH CARET
#SEPARATE EMBEDDING AND MODEL STEPS

custom_import_small <- function(x,y, query, column, feature = "envo_biome_2"){
  temp <- biomformat::read_biom(x)
  OTU_table <- biomformat::biom_data(temp)
  
  sample_data <- readr::read_tsv(y)
  int_data <- dplyr::filter(sample_data, grepl(query, sample_data[column])) %>% #filter data into description
    mutate(feature = pull(sample_data[feature])) #define feature of interest as envo_biome_2
  int_filter <- dplyr::pull(int_data, "#SampleID") #pull sample names with soil in description
  #convert to spare
  sub <- OTU_table[,int_filter] #subset to only soil samples
  sub <- sub[!Matrix::rowSums(sub) == 0 ,] #remove rows with no taxa
  return(sub)
}
custom_import_big <- function(x,y,z, query, column, feature = "envo_biome_2"){
  #read with fread for speed
  OTU_table <- data.table::fread(x)
  OTU_table <- to_sparse(OTU_table)
  #reshape taxa rows
  colnames(OTU_table)[1] <- "temp" #taxa row names is read as first column
  rownames(OTU_table) <- OTU_table$temp #set rownames to that column
  OTU_table <- (OTU_table[,-1]) #remove temp column
  #filter to soil samples
  sample_data <- read_table2(y)
  int_data <- dplyr::filter(sample_data, grepl(query, sample_data[column])) %>% #filter data into description
    mutate(feature = pull(sample_data[feature])) #define feature of interest as envo_biome_2
  int_filter <- dplyr::pull(int_data, "#SampleID") #pull sample names with soil in description
  #########
  column_file <- paste(dataFolder, z, sep = "")
  column_names <- readr::read_csv(column_file, col_names = c("temp")) %>%
    pull(temp)
  colnames(OTU_table) <- column_names #add column names
  int_filter <- int_filter[int_filter %in% column_names] #subset names with soil with existing column names in data
  ###
  OTU_table <- to_sparse(OTU_table)
  OTU_int <- OTU_table[,int_filter]
  OTU_sub <- OTU_int[!Matrix::rowSums(OTU_int) == 0 ,]
}


#Add prediction variable for RF training to embeddings
add_predict <- function(x, y, predict){
  res <- crossprod(x = x, y = y) %>%
    t() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "#SampleID") %>%
    left_join(predict, by = "#SampleID") %>%
    select(-c("#SampleID")) %>%
    mutate(feature = as.factor(feature))
  return(res)
}

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
#Python stuff
#table = load_table("emp_cr_gg_13_8.qc_filtered.biom") 
#a = table.ids()  
#b = a.astype("U") 
#b.tofile("out.csv", sep="\n", format="%s")   