#' @export
calculate_density_info <- function(distance_matrix, k) {
    apply(distance_matrix, 2, function(distances) {
        k_smallest_dist <- sort(distances, partial = k)[k]
        nearest_neighbours <- (distances <= k_smallest_dist)
        n = sum(nearest_neighbours, na.rm = T)
        density = n/k_smallest_dist
        list(density = density, nearest_neighbours = nearest_neighbours, distances = distances)
    })
}

#' @export
cluster_lookup <- function(cluster_tbl) {
    # browser()
    cluster_tbl %>%
        select(cluster_id, cluster_merge) %>%
        distinct() %>%
        igraph::graph.data.frame(.,
                                 directed=TRUE,
                                 vertices=NULL) %>%
        igraph::clusters() %>%
        pluck("membership") %>%
        enframe() %>%
        rename(cluster_id = name,
               cluster_merge = value) %>%
        mutate_all(as.integer)


    # tb <-
    #     cluster_tbl %>%
    #     select(cluster_id, cluster_merge) %>%
    #     distinct() %>%
    #     add_count(cluster_id) %>%
    #     filter(!((n > 1) & (cluster_id == cluster_merge))) %>%
    #     select(-n)
    #
    # cluster_lookup_vec <- tb$cluster_merge
    # names(cluster_lookup_vec) <- tb$cluster_id
    # tb_new <- tb %>% mutate(cluster_merge = cluster_lookup_vec[as.character(cluster_merge)])
    #
    # names(tb_new$cluster_merge) <- NULL
    # names(tb$cluster_merge) <- NULL
    #
    # if (identical(tb, tb_new)) {
    #     tb_new
    # } else {
    #     cluster_lookup(tb_new)
    # }
}

#' @export
merge_clusters <- function(df_cluster) {
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




#' @export
cluster_knn <- function(df_aggregated, k, use_cache = TRUE) {
    # message('------------Clustering------------')
    stopifnot("Aggregated_Dataframe" %in% class(df_aggregated))

    # .GlobalEnv$env_report$k <- k

    message("Clustering... \n")

    if (exists("env_dm", envir = globalenv()) && identical(.GlobalEnv$env_dm$df_aggregated, df_aggregated)) {
        if (use_cache) {
            message("Using cached distance matrix... \n")
            df_sequence <- .GlobalEnv$env_dm$df_sequence
            distance_matrix <- .GlobalEnv$env_dm$distance_matrix
        } else {
            rm(env_dm)
            df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()
            message("Calculating distance matrix... \n")
            distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))
        }
    } else {
        if (use_cache) {
            .GlobalEnv$env_dm <- new.env()

            df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()
            message("Calculating distance matrix... \n")
            distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))

            .GlobalEnv$env_dm$df_aggregated <- df_aggregated
            .GlobalEnv$env_dm$df_sequence <- df_sequence
            .GlobalEnv$env_dm$distance_matrix <- distance_matrix

            message("Caching distance matrix... \n")
        } else {
            df_sequence <- df_aggregated %>% convert_to_sequence() %>% ungroup()
            message("Calculating distance matrix... \n")
            distance_matrix <- inter_sequence_distance(df_sequence %>% pull(sequence))
        }
    }



    # step 1 - initialize every *unique* sequence as a cluster
    message("Initializing clusters... \n")
    df_cluster <- df_sequence %>% select(-sequence_formatted) %>% mutate(density_info = calculate_density_info(distance_matrix,
        k), sequence_density = map_dbl(density_info, "density"), cluster_id = row_number(), cluster_density = sequence_density)

    df_cluster <- df_cluster %>% mutate(same_sequence_clusters = map_int(density_info, function(sequence_density_info) {
        distances <- sequence_density_info$distances
        distances[is.na(distances)] <- 0
        min(df_cluster$cluster_id[distances == 0])
    }), cluster_id = same_sequence_clusters) %>% select(-same_sequence_clusters)


    # step 2 - clustering based on criteria
    message("Clustering based on density... \n")
    df_cluster <- df_cluster %>% mutate(cluster_merge = cluster_id, cluster_merge = map2_int(density_info,
        cluster_id, function(sequence_density_info, current_cluster) {
            # browser()
            density <- sequence_density_info$density
            distances <- sequence_density_info$distances

            checks <- (sequence_density_info$nearest_neighbours) & (df_cluster$sequence_density > density) &
                (df_cluster$cluster_id != current_cluster)

            if (any(checks)) {
                candidate_clusters <- df_cluster$cluster_id[checks]
                candidate_clusters[order(distances[checks])][1]
            } else {
                current_cluster
            }
        }))

    df_cluster <- merge_clusters(df_cluster)

    # browser() step 3
    message("Resolving ties... \n")
    df_cluster <- df_cluster %>% mutate(cluster_merge = cluster_id, cluster_merge = pmap_int(list(density_info,
        cluster_id, cluster_density, cluster_id), function(sequence_density_info, current_cluster, current_cluster_density,
        id) {
        density <- sequence_density_info$density
        distances <- sequence_density_info$distances

        checks <- (sequence_density_info$nearest_neighbours) & (df_cluster$sequence_density == density) &
            (df_cluster$cluster_id != current_cluster) & (df_cluster$cluster_density > current_cluster_density)

        if (any(checks)) {
            candidate_clusters <- df_cluster$cluster_id[checks]
            candidate_clusters[order(distances[checks])][1]
        } else {
            current_cluster
        }

    }))

    df_cluster <- merge_clusters(df_cluster)


    df_cluster <- df_cluster %>% group_by(cluster_id) %>% arrange(desc(sequence_density)) %>% select(-sequence_density,
        -density_info, -cluster_density, cluster = cluster_id) %>% nest(df_sequences=c("id","sequence")) %>% mutate(n = map_int(df_sequences,
        nrow), df_sequences = map(df_sequences, function(df_sequence) {
        class(df_sequence$sequence) <- c("Sequence_List", class(df_sequence$sequence))
        names(df_sequence$sequence) <- df_sequence$id
        df_sequence
    })) %>% arrange(desc(n)) %>% ungroup() %>% mutate(cluster = row_number())

    class(df_cluster) <- c("Clustered_Dataframe", class(df_cluster))

    message("----------Done Clustering---------- \n")

    df_cluster
}




#' @export
cluster_kmedoids <- function(df_aggregated, k, use_cache = TRUE, estimate_k = FALSE, k.max = 10) {

  #df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()

  #message("Calculating distance matrix...")
  #distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))
  #distance_matrix[is.na(distance_matrix)] = 0

  #####
  # message('------------Clustering------------')
  stopifnot("Aggregated_Dataframe" %in% class(df_aggregated))

  # .GlobalEnv$env_report$k <- k

  message("Clustering... \n")

  if (exists("env_dm", envir = globalenv()) && identical(.GlobalEnv$env_dm$df_aggregated, df_aggregated)) {
      if (use_cache) {
          message("Using cached distance matrix... \n")
          df_cluster <- .GlobalEnv$env_dm$df_sequence
          distance_matrix <- .GlobalEnv$env_dm$distance_matrix
          distance_matrix[is.na(distance_matrix)] = 0

      } else {
          rm(env_dm)
          df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
          message("Calculating distance matrix... \n")
          distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))
          distance_matrix[is.na(distance_matrix)] = 0
      }
  } else {
      if (use_cache) {
          .GlobalEnv$env_dm <- new.env()

          df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
          message("Calculating distance matrix... \n")
          distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))

          .GlobalEnv$env_dm$df_aggregated <- df_aggregated
          .GlobalEnv$env_dm$df_sequence <- df_cluster
          .GlobalEnv$env_dm$distance_matrix <- distance_matrix

          distance_matrix[is.na(distance_matrix)] = 0

          message("Caching distance matrix... \n")
      } else {
          df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
          message("Calculating distance matrix... \n")
          distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))
          distance_matrix[is.na(distance_matrix)] = 0
      }
  }



  #######

  #######
  # Estimating the optimal number of clusters
  if (estimate_k) {

    message("Estimating Optimal K using the silhouette approach \n")
    res <- fviz_nbclust(distance_matrix, cluster::pam, method = "silhouette", diss = distance_matrix, k.max = k.max)
    best_k <- which.max(res$data$y)

    message("Optimal K = ", best_k, "\n")
    print(res)

    message("Clustering Based on PAM Algorithm \n")
    res <- pam(distance_matrix, k = best_k, diss = TRUE)

  } else {

    message("Clustering Based on PAM Algorithm \n")
    res <- pam(distance_matrix, k = k, diss = TRUE)

  }



  #######

  #message("Clustering Based on PAM Algorithm")
  #res <- pam(distance_matrix, k = k, diss = TRUE)


  df_cluster$cluster_id <- res$cluster



  df_cluster <- df_cluster %>%
                  group_by(cluster_id) %>%
                  select(-sequence_formatted, cluster = cluster_id) %>%
                  nest(df_sequences=c("id","sequence")) %>%
                  mutate(n = map_int(df_sequences, nrow),
                         df_sequences = map(df_sequences, function(df_sequence) {
                           class(df_sequence$sequence) <- c("Sequence_List", class(df_sequence$sequence))
                           names(df_sequence$sequence) <- df_sequence$id
                           df_sequence
                  })) %>%
                  arrange(desc(n)) %>%
                  ungroup() #%>%                   mutate(cluster = row_number())

  class(df_cluster) <- c("Clustered_Dataframe", class(df_cluster))

  message("----------Done Clustering----------")

  df_cluster

}








#' @export
find_optimal_k <- function(df_aggregated, clustering = 'k-nn', max_k = 10, use_cache = TRUE) {


  # Using cached information if present, otherwise calculating it and storing
  #   it as an environment variable
  stopifnot("Aggregated_Dataframe" %in% class(df_aggregated))

  message("Clustering... \n")

  if (exists("env_dm", envir = globalenv()) && identical(.GlobalEnv$env_dm$df_aggregated, df_aggregated)) {
    if (use_cache) {
      message("Using cached distance matrix... \n")
      df_cluster <- .GlobalEnv$env_dm$df_sequence
      distance_matrix <- .GlobalEnv$env_dm$distance_matrix
      distance_matrix[is.na(distance_matrix)] = 0

    } else {
      rm(env_dm)
      df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
      message("Calculating distance matrix... \n")
      distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))
      distance_matrix[is.na(distance_matrix)] = 0
    }
  } else {
    if (use_cache) {
      .GlobalEnv$env_dm <- new.env()

      df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
      message("Calculating distance matrix... \n")
      distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))

      .GlobalEnv$env_dm$df_aggregated <- df_aggregated
      .GlobalEnv$env_dm$df_sequence <- df_cluster
      .GlobalEnv$env_dm$distance_matrix <- distance_matrix

      distance_matrix[is.na(distance_matrix)] = 0

      message("Caching distance matrix... \n")
    } else {
      df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
      message("Calculating distance matrix... \n")
      distance_matrix <- inter_sequence_distance(df_cluster %>% pull(sequence))
      distance_matrix[is.na(distance_matrix)] = 0
    }
  }





  if (clustering == "k-nn") {

    algo = "K-Nearest Neighbors"

    sil <- sapply(2:max_k, function(k) {

              clustered <- df_aggregated %>% cluster_knn(k = k)
              cluster_id <- clustered %>% unnest(df_sequences) %>% select(id, cluster) %>% arrange(id)

              mean(silhouette(cluster_id$cluster, dmatrix = distance_matrix)[,3])

    })

  } else if (clustering == "k-medoids") {

    algo = "K-Medoids"

    sil <- sapply(2:max_k, function(k) {

      clustered <- df_aggregated %>% cluster_kmedoids(k = k)
      cluster_id <- clustered %>% unnest(df_sequences) %>% select(id, cluster) %>% arrange(id)

      mean(silhouette(cluster_id$cluster, dmatrix = distance_matrix)[,3])

    })

  } else {

    stop("Only clustering values of 'k-nn' and 'k-medoids' are supported right now.")

  }


  plot(c(sil), type = "o", col = "#20B2AA", bty = "l", oma = c(2, 3, 4, 5),
       xlab = "k Value",
       ylab = "Average silhouette width") +

    mtext(paste0("Optimal K Plot for ", algo), line = 2, adj = 0, cex = 1.5) +

    mtext(paste0("k =", which.max(sil), "; Max average silhouette width = ", round(max(sil), digits = 3)),
          line = .75, adj = 0) +

    abline(v = c(which.max(sil), max(sil)),
           lty = "dashed", lwd = .5, col = "#20B2AA")


}
