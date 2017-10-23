library(tidyverse)
library(TraMineR)
library(stringr)
library(lubridate)

data("mvad")
df <- as_tibble(mvad)

df <-
df %>%
  select(-(2:14)) %>%
  gather(time, event, -id) %>%
  mutate(period =  str_c("01 ", time) %>%  str_replace("\\."," ")) %>%
  select(id, period, event) %>%
  mutate(period = readr::parse_date(period, format = "%d %b %y") %>%
    as.character())


df %>%
  aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 1) %>%
  cluster_knn(k = 50, use_cache = T)

df %>%
  aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 1) %>%
  cluster_knn(k = 15, use_cache = T)

df %>%
  aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 2) %>%
  cluster_knn(k = 25, use_cache = T)


y <- x %>%
  get_weighted_sequence()

y %>%
  filter_pattern(threshold = 0.3, pattern_name = "signal") %>%
  format_sequence() %>%
    View()

y <- x %>%
  filter_pattern(threshold = 0.6) %>%
    format_sequence()

x %>%
  get_weighted_sequence()

y$consensus_pattern


res <-
pre_agg_full %>%
  pre_aggregated() %>%
    cluster_knn(k = 4) %>%
      get_weighted_sequence()

res <-
res %>%
  mutate(weighted_sequence =
           map(weighted_sequence, function(x){
             class_it(x,"W_Sequence")
           }))

class(res$weighted_sequence) <- c("W_Sequence_List", class(res$weighted_sequence))
res <- class_it(res, "W_Sequence_Dataframe")
res %>% filter_pattern() %>%
  format_sequence(compare = T, n_as_percent = F) %>%
    View()

x <- pre_agg_test %>%
  pre_aggregated() %>%
  cluster_knn(k = 2) %>%
    .$df_sequences %>%
      .[[1]]
x <- x[c(1,3,2,5,4,6,7),]
y <- x$sequence
names(y) <- x$id
class(y) <- c("Sequence_List",class(y))
get_weighted_sequence(y)





pre_agg_full %>%
  pre_aggregated() %>%
  cluster_knn(k = 4, use_cache = T)



















# library(tidyverse)
#
# # pre_agg_full <- read_csv("./data/pre_agg.csv", col_names = c("id", "period", "event"), skip = 1)
# # date_agg <- read_csv("./data/date_agg.csv")
#
# pre_agg %>%
# cluster_knn()
#
#
# pre_agg %>%
#   pre_aggregated() %>%
#   cluster_knn(k = 2) %>%
#     get_weighted_sequence() %>%
#       filter_pattern(pattern_name = "variation", threshold = 0.5, noise = 5) %>%
#        filter_pattern(pattern_name = "consensus", threshold = 0.75)
#
#
# date_agg %>%
#   aggregate_sequences(format = "%Y-%m-%d") %>%
#   cluster_knn() %>%
#    filter_pattern(threshold = 0.1) %>%
#      format_sequence()
#
# x <-
# pre_agg %>%
#   pre_aggregated() %>%
#   cluster_knn(k = 2) %>%
#   get_weighted_sequence() %>%
#   filter_pattern(pattern_name = "variation", threshold = 0.25, noise = 2) %>%
#   filter_pattern(pattern_name = "consensus", threshold = 0.25, noise = 5)
#
#
# x
# # %>%
#   # pull(consensus)
#
#
#
#
# date_agg %>%
#   aggregate_sequences(format = "%Y-%m-%d") %>%
#     convert_to_sequence() %>%
#      .$sequence %>%
#       inter_sequence_distance()
