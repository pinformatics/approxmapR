#dummy clustering

cluster_knn_api <- function(df){
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



cluster_knn <- function(df_aggregated, k) {
  # message("------------Clustering------------")
  stopifnot(class(df_aggregated) == "Aggregated_Dataframe")

  df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()

  message("Calculating distance matrix...")
  distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))


  #step 1 - initialize every *unique* sequence as a cluster
  message("Initializing clusters...")
  df_cluster <-
    df_sequence %>%
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


  #step 2 - clustering based on criteria
  message("Clustering based on density...")
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





  for(i in 1:length(sequences)) {
    (current_cluster <- cluster_info$cluster[i])
    (cluster_to_merge <- as.integer(cluster_info %>%
                                      dplyr::rename(id_new = id,density_array_new = density_array) %>%
                                      mutate(distances = distance_matrix[i,]) %>%
                                      filter(density_info[[i]]$NearestSequences) %>%
                                      filter(density_array_new > density_array[i]) %>%
                                      filter(cluster != current_cluster) %>%
                                      arrange(distances) %>%
                                      slice(1) %>% select(cluster)))

    if(!is.na(cluster_to_merge)) {
      cluster_info[cluster_info$cluster == current_cluster,]$cluster <- cluster_to_merge
      cluster_info[cluster_info$cluster == cluster_to_merge,]$cluster_density <- max(cluster_info[cluster_info$cluster == cluster_to_merge,]$cluster_density)
    }
  }



  # print(length(unique(cluster_info$cluster)))
  # print(table(cluster_info$cluster))
  #
  #step 3
  message("Resolving ties...")
  for(i in 1:length(sequences)) {
    (current_cluster <- cluster_info$cluster[i])
    (current_cluster_density <- cluster_info$cluster_density[i])
    current_sequence_density <- cluster_info$density_array[i]
    nearest_neighbours <- cluster_info %>%
      dplyr::rename(id_new = id,density_array_new = density_array) %>%
      filter(density_info[[i]]$NearestSequences) %>%
      filter(cluster != current_cluster)
    neighbour_density_check <- nearest_neighbours %>% filter(density_array_new > current_sequence_density) %>% nrow()

    if(!neighbour_density_check) {
      cluster_to_merge <- as.integer(nearest_neighbours %>%
                                       filter(near(density_array_new,current_sequence_density,0.1)) %>%
                                       filter(cluster_density > current_cluster_density) %>%
                                       arrange(desc(cluster_density)) %>% slice(1) %>%
                                       select(cluster))
      #filter(near(density_array_new,density_array[i],0.001)) %>%
      if(!is.na(cluster_to_merge)) {
        cluster_info[cluster_info$cluster == current_cluster,]$cluster <- cluster_to_merge
        cluster_info[cluster_info$cluster == cluster_to_merge,]$cluster_density <- max(cluster_info[cluster_info$cluster == cluster_to_merge,]$cluster_density)
      }
    }
  }

  cluster2 <- cluster <- cluster_info$cluster
  cluster_unique = unique(cluster)
  for(cl in 1:length(cluster_unique)) {
    cluster2[cluster==cluster_unique[cl]] = cl
  }
  cluster_info$cluster <- cluster2

  res = data.frame("ID" = cluster_info$id , "Density" = cluster_info$density_array, "Cluster" = cluster_info$cluster, "ClusterDensity" = round(cluster_info$cluster_density,2))
  message("----------Done Clustering----------")
  return(res)
}
