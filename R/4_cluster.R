# #dummy clustering
#
# cluster_knn_api <- function(df){
#   set.seed(1)
#   ids <- unique(df$id)
#   clust <- sample(1:3, length(ids), replace = T)
#   density <- abs(rnorm(length(ids), 1, 3))
#   cluster_info <- tibble(id = ids, cluster = clust, density = density)
#   df_clustered <-
#     df %>%
#     left_join(cluster_info, by = "id") %>%
#     group_by(cluster) %>%
#     arrange(density,id) %>%
#     nest(-density, .key = df_list) %>%
#     mutate(n = map_int(df_list, nrow),
#            df_list = map(df_list,
#                          function(df_group){
#                            if("Aggregated_Dataframe" %in% class(df)){
#                              class(df_group) <- class(df)
#                            }
#                            df_group
#                          })
#     ) %>%
#     arrange(desc(n)) %>%
#     mutate(cluster = row_number())
#   class(df_clustered) <- c("Clustered_Dataframe", class(df_clustered))
#
#   df_clustered
# }

calculate_density_info <- function(distance_matrix, k) {
  apply(distance_matrix, 2, function(distances){
    k_smallest_dist <- sort(distances, partial = k)[k]
    nearest_neighbours <- distances <= k_smallest_dist
    n = sum(nearest_neighbours, na.rm = T)
    density = n/k_smallest_dist
    list(density = density,
         nearest_neighbours = nearest_neighbours,
         distances = distances)
  })
}

cluster_lookup <- function(cluster_tbl){
  tb <- cluster_tbl %>%
    select(cluster_id, cluster_merge) %>%
      distinct()
  cluster_lookup_vec <- tb$cluster_merge
  names(cluster_lookup_vec) <- tb$cluster_id
  tb_new <-
    tb %>%
    mutate(cluster_merge = cluster_lookup_vec[as.character(cluster_merge)])

  names(tb_new$cluster_merge) <- NULL
  names(tb$cluster_merge) <- NULL

  if(identical(tb, tb_new)){
    tb_new
  } else{
    cluster_lookup(tb_new)
  }
}

merge_clusters <- function(df_cluster){
  cluster_lookup <-
    df_cluster %>%
    select(cluster_id, cluster_merge) %>%
    cluster_lookup()

  df_cluster <-
    df_cluster %>%
    select(-cluster_merge) %>%
    left_join(cluster_lookup, by = "cluster_id") %>%
    mutate(cluster_id = cluster_merge) %>%
    select(-cluster_merge) %>%
    group_by(cluster_id) %>%
    mutate(cluster_density = max(cluster_density)) %>%
    ungroup()
  names(df_cluster$sequence) <- df_cluster$id

  df_cluster
}

# x <- function(){
#   $xyz <- 2
# }



cluster_knn <- function(df_aggregated, k, use_cache = TRUE) {
  # message("------------Clustering------------")
  stopifnot("Aggregated_Dataframe" %in% class(df_aggregated))

  # .GlobalEnv$env_report$k <- k

  message("Clustering...")

  if(exists("env_dm", envir = globalenv()) &&
     identical(.GlobalEnv$env_dm$df_aggregated, df_aggregated)) {
    if(use_cache){
        message("Using cached distance matrix...")
        df_sequence <-  .GlobalEnv$env_dm$df_sequence
        distance_matrix <- .GlobalEnv$env_dm$distance_matrix
    } else{
      rm(env_dm)
      df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()
      message("Calculating distance matrix...")
      distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))
    }
  } else{
    if(use_cache){
      .GlobalEnv$env_dm <- new.env()

      df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()
      message("Calculating distance matrix...")
      distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))

      .GlobalEnv$env_dm$df_aggregated <- df_aggregated
      .GlobalEnv$env_dm$df_sequence <- df_sequence
      .GlobalEnv$env_dm$distance_matrix <- distance_matrix

      message("Caching distance matrix...")
    } else{
      df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()
      message("Calculating distance matrix...")
      distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))
    }
  }


  #step 1 - initialize every *unique* sequence as a cluster
  message("Initializing clusters...")
  df_cluster <-
    df_sequence %>%
    select(-sequence_formatted) %>%
    mutate(density_info = calculate_density_info(distance_matrix,k),
           sequence_density = map_dbl(density_info, "density"),
           cluster_id = row_number(),
           cluster_density = sequence_density)

  df_cluster <-
    df_cluster %>%
      mutate(same_sequence_clusters =
               map_int(density_info,
                       function(sequence_density_info){
                        distances <- sequence_density_info$distances
                        distances[is.na(distances)] <- 0
                        min(df_cluster$cluster_id[distances == 0])
                     }),
             cluster_id = same_sequence_clusters
             ) %>%
      select(-same_sequence_clusters)


  # step 2 - clustering based on criteria
  message("Clustering based on density...")
  df_cluster <-
  df_cluster %>%
    mutate(cluster_merge = cluster_id,
           cluster_merge =
             map2_int(density_info, cluster_id,
                      function(sequence_density_info, current_cluster){
                        density <- sequence_density_info$density
                        distances <- sequence_density_info$distances

                        checks <- (sequence_density_info$nearest_neighbours) &
                          (df_cluster$sequence_density > density) &
                          (df_cluster$cluster_id != current_cluster)

                        if(sum(checks) != 0){
                          candidate_clusters <- df_cluster$cluster_id[checks]
                          candidate_clusters[order(distances[checks])][1]
                        } else {
                          current_cluster
                        }
                      })
           )

  df_cluster <- merge_clusters(df_cluster)

  #step 3
  message("Resolving ties...")
  df_cluster <-
    df_cluster %>%
    mutate(cluster_merge = cluster_id,
           cluster_merge =
             pmap_int(list(density_info, cluster_id, cluster_density),
                      function(sequence_density_info, current_cluster, current_cluster_density){
                        density <- sequence_density_info$density
                        distances <- sequence_density_info$distances

                        checks <- (sequence_density_info$nearest_neighbours) &
                          (df_cluster$sequence_density == density) &
                          (df_cluster$cluster_id != current_cluster) &
                          (df_cluster$cluster_density > current_cluster_density)

                        if(sum(checks) != 0){
                          candidate_clusters <- df_cluster$cluster_id[checks]
                          candidate_clusters[order(distances[checks])][1]
                        } else {
                          current_cluster
                        }

                      })
    )


  df_cluster <- merge_clusters(df_cluster)

  df_cluster <-
    df_cluster %>%
      group_by(cluster_id) %>%
      arrange(desc(sequence_density)) %>%
      select(-sequence_density,
             -density_info,
             -cluster_density,
             cluster = cluster_id) %>%
      nest(.key = df_sequences) %>%
      mutate(n = map_int(df_sequences, nrow),
             df_sequences =
               map(df_sequences,
                   function(df_sequence){
                     class(df_sequence$sequence) <-
                       c("Sequence_List",class(df_sequence$sequence))
                     df_sequence
                   })) %>%
      arrange(desc(n)) %>%
      mutate(cluster = row_number())

  class(df_cluster) <- c("Clustered_Dataframe", class(df_cluster))

  message("----------Done Clustering----------")

  df_cluster
}
