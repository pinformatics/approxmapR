cluster_knn <- function(df){
  set.seed(1)
  ids <- unique(df$id)
  clust <- sample(1:3, length(ids), replace = T)
  density <- abs(rnorm(length(ids), 1, 3))
  cluster_info <- tibble(id = ids, cluster = clust, density = density)
  df_clustered <-
    df %>%
    left_join(cluster_info, by = "id") %>%
    group_by(cluster) %>%
    arrange(density,id) %>%
    nest(-density, .key = df_list) %>%
    mutate(n = map_int(df_list, nrow),
           df_list = map(df_list,
                         function(df_group){
                           if("Aggregated_Dataframe" %in% class(df)){
                             class(df_group) <- class(df)
                           }
                           df_group
                         })
    ) %>%
    arrange(desc(n)) %>%
    mutate(cluster = row_number())
  class(df_clustered) <- c("Clustered_Dataframe", class(df_clustered))

  df_clustered
}
