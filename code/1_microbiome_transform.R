##################################################
## Project: NLKT soil
## Script purpose: Read in microbiome RDS from 0_ file and create property embeddings
## Date: Mon May 18 09:37:10 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries
##################################################
#load packages here
library("tidyverse")
rm(list = ls())
dataFolder <- "data/"

OTU_sub <- readRDS(paste(dataFolder, "OTU_sub.RDS", sep = ""))
taxon <- read_table2(paste(dataFolder, "97_otu_taxonomy.txt", sep = ""), col_names = F) %>%
  mutate_all(funs(gsub("[a-z]__|;", "", .))) %>%
  mutate_all(list(~na_if(.,""))) 

taxas <- readRDS(paste(dataFolder, "taxa.RDS", sep = ""))
colnames(taxon) <- c("taxa", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

rdsfile <- "rf_abundance_post"
file <- readRDS(paste(dataFolder, rdsfile, ".RDS", sep = ""))
file2 <- as.matrix(file[,-1])
file2[file2 < 0.20] <- 0
file2 <- apply(file2, MARGIN = 1, FUN = function(x)return(x/sum(x))) %>% t() %>% as.data.frame()
file2 <- file2[,colSums(file2, na.rm = T)!=0]
file2$response <- file$response
test <- file2 %>%
  drop_na() %>%
  rownames_to_column("sample") %>%
  reshape2::melt(id.vars = c("response", "sample")) %>%
  ggplot(aes(x = sample, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~response, scales = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(title = "Bacteria Abundance by Envo_biome x Phylum")
ggsave("output/comparison.svg", height = 4, width = 8)
