list.files()
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(tidyverse)
library(tidymodels)
res <- readRDS("envo_results.RDS")


metrics <- res %>%
  unnest(.metrics) %>%
  select(.model, .embed, .metric, .estimate) %>%
  ggplot(aes(x = .metric, y = .estimate, fill = .metric)) + facet_grid(.model~.embed) +
           stat_summary(fun.data = quantiles_95, geom="boxplot")+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
            
metrics
head(metrics)
metrics2 <- res %>%
  unnest(.metrics) %>%
  select(.model, .embed, .metric, .estimate) %>%
  ggplot(aes(x = .embed, y = .estimate, fill = .metric)) + facet_grid(.model~.metric) +
  stat_summary(fun.data = quantiles_95, geom="boxplot") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "metrics2.svg", width = 10, height = 6)
metrics2


metrics3 <- res %>%
  select(.model, .embed, .metrics2) %>%
  mutate_if(is.list, simplify_all) %>%    # flatten each list element internally 
  unnest(.metrics2) %>%
  reshape2::melt(id.vars = c(".model", ".embed", ".variable")) %>%
  #filter(variable != "TP" & variable != "FP" & variable != "FN" & variable != "TN") %>%
  filter(variable == "f_meas") %>%
  ggplot(aes(x = .embed, y = value, fill = .embed)) + facet_grid(.model~.variable) +
    stat_summary(fun.data = quantiles_95, geom="boxplot") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "metrics3.svg", width = 10, height = 6)
metrics3


metrics4 <- res %>%
  select(.model, .embed, pr_curves) %>%
  mutate_if(is.list, simplify_all) %>%    # flatten each list element internally 
  unnest(pr_curves) %>%
  ggplot(aes(x = recall, y = precision, color = .level)) + facet_grid(.model~.embed) +
  geom_path() +
  coord_equal() +
  theme_bw()
metrics4

metrics5 <- res %>%
  select(.model, .embed,roc_curves) %>%
  mutate_if(is.list, simplify_all) %>%    # flatten each list element internally 
  unnest(roc_curves) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = .level)) + facet_grid(.model~.embed) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw()
metrics5
