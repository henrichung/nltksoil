list.of.packages <- c("tidyverse", "parallel", "parsnip", "rsample", "yardstick", "reshape2", 
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel", 
                      "text2vec", "Rcpp", "profvis","Matrix", "themis", "stringi", "data.table", 
                      "tictoc", "gridExtra", "ggpubr", "rstatix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))
rm(list = ls())
source("code/nltksoil_functions.R")
results_tune <- readRDS("data/results_tune2.RDS")
#results_validation <- readRDS("data/results_validation2.RDS")

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
#read taxonomy data
taxonomy <- read.table("data/97_otu_taxonomy.txt", fill = TRUE) %>%
  as_tibble()
#read sample data
sampledata <- read.csv("data/emp_qiime_mapping_qc_filtered.tsv", sep = "\t") %>% 
  filter(grepl("soil", Description))


#Response class distributions
p0 <- sampledata %>%
  pull(envo_biome_2) %>%
  table() %>% 
  as.data.frame() %>%
  filter(Freq > 0) %>% 
  ggplot(aes(x = ., y = Freq, fill = .)) + geom_bar(stat = "identity") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(title = "Response class distributions (Envo_biome_2)",
       subtitle = "soil samples only, n = 1728") +
  ylab("OTU Count") 
p0
ggsave(p0, filename = "outputs/p0.svg", height = 6, width = 8)


#METHOD COMPARISON 50 COMPONENTS BALANCED ACCURACY
p1 <- results_tune %>% 
  filter(.taxa == "species" & .model == "randomforest") %>%
  unnest(.metrics) %>%
  filter(.metric == "f_meas") %>%
  filter(.no_components == "50") %>%
  ggplot(aes(x = as.factor(min_n), y = .estimate, fill = .method )) + 
  geom_boxplot() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
p1

p1b <- results_tune %>% 
  filter(.taxa == "species" & .model == "randomforest") %>%
  unnest(.extracts) %>% unnest(.extracts) %>%
  filter(.no_components == "50") %>%
  ggplot(aes(x = as.factor(min_n), y = .extracts, fill = .method )) + 
  geom_boxplot() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
p1b

p2 <- results_tune %>% 
  filter(.taxa == "species" & .model == "randomforest") %>%
  unnest(.metrics) %>%
  filter(.metric == "f_meas") %>%
  filter(.no_components == "50") %>%
  ggplot(aes(x = as.factor(try), y = .estimate, fill = .method )) + 
  geom_boxplot() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
p2

p2b <- results_tune %>% 
  filter(.taxa == "species" & .model == "randomforest") %>%
  unnest(.extracts) %>% unnest(.extracts) %>%
  filter(.no_components == "50") %>%
  ggplot(aes(x = as.factor(mtry), y = .extracts, fill = .method )) + 
  geom_boxplot() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
p2b

p3 <- results_tune %>% 
  filter(.taxa == "species" & .model == "randomforest") %>%
  unnest(.metrics) %>%
  filter(.metric == "f_meas") %>%
  filter(.no_components == "50") %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method )) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p3

p4 <- results_tune %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.model == "randomforest") %>%
  filter(.no_components == "50") %>%
  unnest(.metrics) %>%
  filter(.metric == "f_meas") %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
p4

p5 <- results_tune %>%
  filter(.method == "transformword" | .method == "transformpca") %>%
  filter(.model == "randomforest") %>%
  unnest(.metrics) %>%
  filter(.metric == "f_meas") %>%
  ggplot(aes(x = .no_components, y = .estimate, fill = .method)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
p5

p5b <- results_tune %>%
  filter(.method == "transformword" | .method == "transformpca") %>%
  filter(.model == "randomforest") %>%
  unnest(.extracts) %>% unnest(.extracts) %>%
  ggplot(aes(x = .no_components, y = .extracts, fill = .method)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
p5b

x1 <- results_tune %>% 
  filter(.taxa == "species" & .model == "xgboost") %>%
  unnest(.metrics) %>%
  filter(.metric == "f_meas") %>%
  filter(.no_components == "50") %>%
  ggplot(aes(x = as.factor(min_n), y = .estimate, fill = .method )) + 
  geom_boxplot() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.method)
x1

test <- head(results_tune) %>% unnest(.predictions)
