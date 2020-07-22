list.of.packages <- c("tidyverse", "parallel", "parsnip", "rsample", "yardstick", "reshape2", 
                      "dials", "tune", "furrr", "future", "workflows", "recipes", "doParallel", 
                      "text2vec", "Rcpp", "profvis","Matrix", "themis", "stringi", "data.table", 
                      "tictoc", "gridExtra", "ggpubr", "rstatix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))
rm(list = ls())
source("code/nltksoil_functions.R")
#results_tune <- readRDS("data/results_tune.RDS")
results_validation <- readRDS("data/results_validation.RDS")

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
p1data <- results_validation %>% filter(.model == "randomforest" & .metric == "f_meas" & .taxa == "species" & .no_components == 50) %>% distinct()
p1 <- p1data %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison", 
       subtitle = "50 components, KAP",
       caption = "95% CI")
p1
ggsave(p1, filename =  "outputs/p1kap.svg", height = 6, width = 10)

p1b <- results_validation %>%
  filter(.model == "randomforest" & .metric == "kap" & .taxa == "species" & .no_components == 50) %>%
  filter(.method != "pca" & .method != "word") %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison", 
       subtitle = "50 components, KAP",
       caption = "95% CI")
p1b
ggsave(p1b, filename =  "outputs/p1bkap.svg", height = 6, width = 10)

p2 <- results_validation %>%
  filter(.model == "randomforest" & .metric == "kap" & .taxa == "species") %>%
  filter(.method != "pca" & .method != "word") %>%
  distinct() %>%
  ggplot(aes(x = .no_components, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  facet_wrap(~.method) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison by # of components", 
       subtitle = "kap",
       caption = "95% CI")
p2
ggsave(p2, filename = "outputs/p2kap.svg", height = 6, width = 8)

p2b <- results_validation %>%
  filter(.model == "randomforest" & .metric == "kap" & .taxa == "species") %>%
  filter(.method == "pca" | .method == "word") %>%
  distinct() %>%
  ggplot(aes(x = .no_components, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  facet_wrap(~.method) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison by # of components", 
       subtitle = "kap",
       caption = "95% CI")
p2b
ggsave(p2b, filename = "outputs/p2bkap.svg", height = 6, width = 8)

p3 <- results_validation %>%
  filter(.model == "randomforest" & .metric == "bal_accuracy" & .taxa == "species") %>%
  filter(.method != "pca" & .method != "word") %>%
  distinct() %>%
  ggplot(aes(x = .no_components, y = .oob_error, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  facet_wrap(~.method) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison", 
       subtitle = "OOB Error",
       caption = "95% CI")
p3
ggsave(p3, filename = "outputs/p3.svg", height = 6, width = 8)

p4 <- results_validation %>%
  filter(.model == "randomforest" & .metric == "kap" & .taxa == "species" & .no_components == 200) %>%
  filter(.method != "word") %>%
  unnest(.oob_error) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison", 
       subtitle = "200 components, kap",
       caption = "95% CI")
p4
ggsave(p4, filename = "outputs/p4.svg", height = 6, width = 8)

p6 <- results_validation %>%
  filter(.model == "randomforest" & .taxa == "species" & .no_components == 200) %>%
  filter(.method != "word") %>%
  filter(.metric %in% c("kap", "accuracy", "bal_accuracy", "f_meas", "npv", "ppv", "sens", "spec")) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~.metric) +
  labs(title = "Methods Comparison", 
       subtitle = "200 components, all metrics",
       caption = "95% CI")
p6
ggsave(p6, filename = "outputs/p6.svg", height = 6, width = 8)

p6b <- results_validation %>%
  filter(.model == "randomforest" & .taxa == "species" & .no_components == 200) %>%
  filter(.metric %in% c("kap", "accuracy", "bal_accuracy", "f_meas", "npv", "ppv", "sens", "spec")) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~.metric) +
  labs(title = "Methods Comparison", 
       subtitle = "200 components, all metrics",
       caption = "95% CI")
p6b
ggsave(p6b, filename = "outputs/p6b.svg", height = 6, width = 8)

p7 <- results_validation %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.metric == "kap") %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison",
       subtitle = "kap")
p7
ggsave(p7, filename = "outputs/p7kap.svg", height = 6, width = 8)

p8 <- results_validation %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.metric == "kap") %>%
  unnest(.oob_error) %>%
  distinct() %>%
  ggplot(aes(x = .taxa, y = .oob_error, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison",
       subtitle = "OOB_error")
p8
ggsave(p8, filename = "outputs/p8.svg", height = 6, width = 8)
p9 <- results_validation %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.metric %in% c("accuracy", "bal_accuracy", "f_meas", "npv", "ppv", "sens", "spec", "kap")) %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison") +
  facet_wrap(~.metric) 
p9
ggsave(p9, filename = "outputs/p9kap.svg", height = 6, width = 8)

p10 <- results_validation %>%
  filter(.metric == "kap") %>%
  filter(.taxa != "species") %>%
  distinct() %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  facet_wrap(.~.no_components) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison with all methods")
p10
ggsave(p10, filename = "outputs/p10kap.svg", height = 6, width = 8)
###########


p9 <- results_validation %>%
  filter(.model == "randomforest" & .taxa == "species" & .no_components == 200) %>%
  filter(.method != "pca" & .method != "word") %>%
  filter(.metric == "accuracy") %>%
  distinct() %>%
  unnest(.custom) %>%
  filter(.metric2 %in% c("kap", "PPV", "NPV", "sensitivity", "specificity",  "accuracy", "balance_accuracy", "f_meas"))
p9

p9a <- p9 %>%
  filter(.class == "anthropogenic terrestrial biome") %>%
  ggplot(aes(x = .method, y = .value, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.metric2, scales = "free") +
  labs(title = "Anthropogenic terrestrial specific")
p9a
ggsave(p9a, filename = "outputs/p9a.svg", height = 6, width = 8)

p9b <- p9 %>%
  filter(.class != "anthropogenic terrestrial biome") %>%
  ggplot(aes(x = .method, y = .value, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.metric2~.class, scales = "free") +
  labs(title = "the other biomes")
p9b
ggsave(p9b, filename = "outputs/p9b.svg", height = 6, width = 8)
###########
###########
###########
#METHOD COMPARISON 50 COMPONENTS BALANCED ACCURACY
x1 <- results_validation %>%
  #filter(.model == "xgboost" & .metric == "kap" & .taxa == "species" & .no_components == 50) %>%
  filter(.metric == "kap" & .taxa == "species" & .no_components == 50) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.model) + 
  labs(title = "Methods Comparison", 
       subtitle = "50 components, KAP",
       caption = "95% CI")
x1
ggsave(p1, filename =  "outputs/p1kap.svg", height = 6, width = 10)

x1b <- results_validation %>%
  filter(.metric == "kap" & .taxa == "species" & .no_components == 50) %>%
  filter(.method != "pca" & .method != "word") %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~.model) + 
  labs(title = "Methods Comparison", 
       subtitle = "50 components, KAP",
       caption = "95% CI")
x1b
ggsave(p1b, filename =  "outputs/p1bkap.svg", height = 6, width = 10)

x2 <- results_validation %>%
  filter(.metric == "kap" & .taxa == "species") %>%
  filter(.method != "pca" & .method != "word") %>%
  distinct() %>%
  ggplot(aes(x = .no_components, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  facet_wrap(~.method) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~.model) + 
  labs(title = "Methods Comparison by # of components", 
       subtitle = "kap",
       caption = "95% CI")
x2
ggsave(p2, filename = "outputs/p2kap.svg", height = 6, width = 8)

x2b <- results_validation %>%
  filter(.metric == "kap" & .taxa == "species") %>%
  filter(.method == "pca" | .method == "word") %>%
  distinct() %>%
  ggplot(aes(x = .no_components, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  facet_grid(.model~.method) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison by # of components", 
       subtitle = "kap",
       caption = "95% CI")
x2b
ggsave(p2b, filename = "outputs/p2bkap.svg", height = 6, width = 8)

x3 <- results_validation %>%
  filter(.metric == "bal_accuracy" & .taxa == "species") %>%
  filter(.method != "pca" & .method != "word") %>%
  distinct() %>%
  unnest(.oob_error) %>%
  ggplot(aes(x = .no_components, y = .oob_error, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  facet_wrap(.model~.method) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison", 
       subtitle = "OOB Error",
       caption = "95% CI")
x3
ggsave(p3, filename = "outputs/p3.svg", height = 6, width = 8)

x4 <- results_validation %>%
  filter(.model == "randomforest" & .metric == "kap" & .taxa == "species" & .no_components == 200) %>%
  filter(.method != "word") %>%
  unnest(.oob_error) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Methods Comparison", 
       subtitle = "200 components, kap",
       caption = "95% CI")
p4
ggsave(p4, filename = "outputs/p4.svg", height = 6, width = 8)

x6 <- results_validation %>%
  filter(.taxa == "species" & .no_components == 200) %>%
  filter(.method != "word") %>%
  filter(.metric %in% c("kap", "accuracy", "bal_accuracy", "f_meas", "npv", "ppv", "sens", "spec")) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .model)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(~.metric) +
  labs(title = "Methods Comparison", 
       subtitle = "200 components, all metrics",
       caption = "95% CI")
x6
ggsave(p6, filename = "outputs/p6.svg", height = 6, width = 8)

p6b <- results_validation %>%
  filter(.model == "randomforest" & .taxa == "species" & .no_components == 200) %>%
  filter(.metric %in% c("kap", "accuracy", "bal_accuracy", "f_meas", "npv", "ppv", "sens", "spec")) %>%
  distinct() %>%
  ggplot(aes(x = .method, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~.metric) +
  labs(title = "Methods Comparison", 
       subtitle = "200 components, all metrics",
       caption = "95% CI")
p6b
ggsave(p6b, filename = "outputs/p6b.svg", height = 6, width = 8)

p7 <- results_validation %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.metric == "kap") %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison",
       subtitle = "kap")
p7
ggsave(p7, filename = "outputs/p7kap.svg", height = 6, width = 8)

p8 <- results_validation %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.metric == "kap") %>%
  unnest(.oob_error) %>%
  distinct() %>%
  ggplot(aes(x = .taxa, y = .oob_error, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison",
       subtitle = "OOB_error")
p8
ggsave(p8, filename = "outputs/p8.svg", height = 6, width = 8)
p9 <- results_validation %>%
  filter(.method == "base" | .method == "abundance") %>%
  filter(.metric %in% c("accuracy", "bal_accuracy", "f_meas", "npv", "ppv", "sens", "spec", "kap")) %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison") +
  facet_wrap(~.metric) 
p9
ggsave(p9, filename = "outputs/p9kap.svg", height = 6, width = 8)

p10 <- results_validation %>%
  filter(.metric == "kap") %>%
  filter(.taxa != "species") %>%
  distinct() %>%
  ggplot(aes(x = .taxa, y = .estimate, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") + 
  facet_wrap(.~.no_components) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Taxa Comparison with all methods")
p10
ggsave(p10, filename = "outputs/p10kap.svg", height = 6, width = 8)
###########


p9 <- results_validation %>%
  filter(.model == "randomforest" & .taxa == "species" & .no_components == 200) %>%
  filter(.method != "pca" & .method != "word") %>%
  filter(.metric == "accuracy") %>%
  distinct() %>%
  unnest(.custom) %>%
  filter(.metric2 %in% c("kap", "PPV", "NPV", "sensitivity", "specificity",  "accuracy", "balance_accuracy", "f_meas"))
p9

p9a <- p9 %>%
  filter(.class == "anthropogenic terrestrial biome") %>%
  ggplot(aes(x = .method, y = .value, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~.metric2, scales = "free") +
  labs(title = "Anthropogenic terrestrial specific")
p9a
ggsave(p9a, filename = "outputs/p9a.svg", height = 6, width = 8)

p9b <- p9 %>%
  filter(.class != "anthropogenic terrestrial biome") %>%
  ggplot(aes(x = .method, y = .value, fill = .method)) + 
  stat_summary(fun.data = quantiles_95, geom="boxplot", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.metric2~.class, scales = "free") +
  labs(title = "the other biomes")
p9b
ggsave(p9b, filename = "outputs/p9b.svg", height = 6, width = 8)