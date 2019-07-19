library(tidyverse)
library(glue)
library(approxmapR)
data("mvad")

#http://web.cs.ucla.edu/~weiwang/paper/SDM03_2.pdf

agg <-
  mvad %>%
  arrange(id, period) %>%
  mutate(event = str_to_lower(event)) %>%
  aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 24)

cluster <-
  agg %>%
  cluster_knn(k=3)

cluster$n %>% sum()



patterns <-
  cluster %>%
  filter_pattern(threshold = 0.5, pattern_name = "consensus") %>%
  filter_pattern(threshold = 0.3, pattern_name = "variation")


patterns %>%
  generate_reports()



  patterns %>%
  print_alignments()

patterns %>%
  save_alignment() %>%
  write_file("alignments.csv")




k <- 1:50

agg_df <-
  df %>%
  mutate(event = str_to_lower(event)) %>%
    aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 2)

read.csv("data/pre_agg_demo.csv", stringsAsFactors = F) %>%
  mutate_if(is.numeric, as.integer) %>%
  pre_aggregated() %>%
  cluster_knn(k=1) %>%
  filter_pattern(threshold = 0.4) %>%
  format_sequence()

k_vals <-
tibble(k = k) %>%
  mutate(cluster_ws = map(k, function(k){
    agg_df %>%
      cluster_knn(k = k) %>%
        filter_pattern(threshold = 0.4)
  }))

agg_seqs <-
mvad %>%
  mutate(event = str_to_lower(event)) %>%
  aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 6) %>%
    cluster_knn(k = 15)

aggs_aligned <-
agg_seqs %>%
  filter_pattern(threshold = 0.3, pattern_name = "variation") %>%
  filter_pattern(threshold = 0.5, pattern_name = "consensus")


formatted <-
  aggs_aligned %>%
    format_sequence(compare = TRUE,
                    html_format = TRUE,
                    truncate_patterns = T)

rmarkdown::render(input = "inst/rmd_w_sequence.Rmd", params = list(input = formatted))

aggs_aligned %>%
  generate_reports(truncate_patterns = T)

aggs_aligned %>%
  generate_reports()


itst <- aggs_aligned %>%
  pull(weighted_sequence) %>%
  .[[1]] %>%
  unclass() %>% .[[1]]
itst$itemset_weight <- NULL
flatten_dfc(itst)


k_vals <-
k_vals %>%
  mutate(num_clusters = map_int(cluster_ws, ~nrow(.)))

k_vals %>%
  filter(between(num_clusters, 5, 10))

k_vals %>%
  filter(k == 28) %>%
    pull(cluster_ws) %>% .[[1]] %>%
      format_sequence() %>%
  View()

k_vals %>%
    ggplot(aes(k, num_clusters)) +
    geom_line()





all_same <- function(x){
  l <- x %>%
  unique() %>%
  length()
  ifelse(l == 1, TRUE, FALSE)
}

x <- agg_seqs %>%
  filter(row_number() == 4)
attr(x,"class") = attr(agg_seqs, "class")
y <- x %>%
  get_weighted_sequence()
y2 <- y %>%
  as.data.frame() %>%
  .[1,4] %>% .[[1]] %>%
  attr("alignments")

all <- agg_seqs %>%
  get_weighted_sequence()

all <- all %>% mutate(cfw = weighted_sequence %>% map2_int(cluster, function(x,y){
  sink(paste0(y,".txt"))
  attr(x, "alignments") %>% print()
  sink()
  attr(x, "alignments") %>%
    map_int(length) %>%
    unique()
}))

map_dbl(y2, length) %>%
  all_same()

dem_aligned <-
  pre_agg_demo %>%
  pre_aggregated() %>%
  cluster_knn(k = 2) %>%
  get_weighted_sequence()

dem_aligned %>%
pull(weighted_sequence)%>%
  map(~attr(.,"alignments"))


# df %>%
#   aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 1) %>%
#   cluster_knn(k = 50, use_cache = T) %>%
#   filter_pattern(threshold = 0.5) %>%
#     format_sequence()
#
# df %>%
#   aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 1) %>%
#   cluster_knn(k = 15, use_cache = T)
#
# df %>%
#   aggregate_sequences(format = "%Y-%m-%d", unit = "month", n_units = 2) %>%
#   cluster_knn(k = 25, use_cache = T)
#
#
# y <- x %>%
#   get_weighted_sequence()
#
# y %>%
#   filter_pattern(threshold = 0.3, pattern_name = "signal") %>%
#   format_sequence() %>%
#     View()
#
# y <- x %>%
#   filter_pattern(threshold = 0.6) %>%
#     format_sequence()
#
# x %>%
#   get_weighted_sequence()
#
# y$consensus_pattern
#
#
# res <-
# pre_agg_full %>%
#   pre_aggregated() %>%
#     cluster_knn(k = 4) %>%
#       get_weighted_sequence()
#
# res <-
# res %>%
#   mutate(weighted_sequence =
#            map(weighted_sequence, function(x){
#              class_it(x,"W_Sequence")
#            }))
#
# class(res$weighted_sequence) <- c("W_Sequence_List", class(res$weighted_sequence))
# res <- class_it(res, "W_Sequence_Dataframe")
# res %>% filter_pattern() %>%
#   format_sequence(compare = T, n_as_percent = F) %>%
#     View()
#
# x <- pre_agg_test %>%
#   pre_aggregated() %>%
#   cluster_knn(k = 2) %>%
#     .$df_sequences %>%
#       .[[1]]
# x <- x[c(1,3,2,5,4,6,7),]
# y <- x$sequence
# names(y) <- x$id
# class(y) <- c("Sequence_List",class(y))
# get_weighted_sequence(y)
#
#
#
#
#
# pre_agg_full %>%
#   pre_aggregated() %>%
#   cluster_knn(k = 4, use_cache = T)
#
#
#
#
#
#
#
#
#
#
#








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
