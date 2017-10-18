library(tidyverse)

t10 <- read_csv("./data/t10.csv", col_names = c("id", "period", "event"), skip = 1, n_max = 350)
inp <- read_csv("./data/inp.csv")

t10 %>%
cluster_knn()


t10 %>%
  pre_aggregated() %>%
  cluster_knn() %>%
    get_weighted_sequence() %>%
      filter_pattern(pattern_name = "variation", threshold = 0.5, noise = 5) %>%
       filter_pattern(pattern_name = "consensus", threshold = 0.75)


inp %>%
  aggregate_sequences(format = "%Y-%m-%d") %>%
  cluster_knn() %>%
   filter_pattern(threshold = 0.1) %>%
     format_sequence()

x <-
t10 %>%
  pre_aggregated() %>%
  cluster_knn() %>%
  get_weighted_sequence() %>%
  filter_pattern(pattern_name = "variation", threshold = 0.25, noise = 2) %>%
  filter_pattern(pattern_name = "consensus", threshold = 0.25, noise = 5)


# %>%
  # pull(consensus)
