library(tidyverse)

pre_agg_test %>%
  pre_aggregated() %>%
    cluster_knn(k = 2) %>%
      filter_pattern(threshold = 0.4) %>%
        format_sequence()


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
