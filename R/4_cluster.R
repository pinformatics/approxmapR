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
cluster_kmedoids <- function(df_aggregated, k, use_cache = TRUE) {

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



  #######################
  # Clustering the data #
  #######################
  message("Clustering Based on PAM Algorithm \n")
  res <- pam(distance_matrix, k = k, diss = TRUE)


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
                  ungroup()

  class(df_cluster) <- c("Clustered_Dataframe", class(df_cluster))

  message("----------Done Clustering----------")

  df_cluster

}








#' @export
find_optimal_k <- function(df_aggregated, clustering = 'k-nn', min_k = 2, max_k = 10, validation_measure = 'silhouette',
                           use_cache = TRUE,
                           save_graph = TRUE, graph_file_name = "Optimal K Graph.png", size_width = 855, size_height = 317, output_directory = "~") {

  # Testing valid parameter entries
  if (!clustering %in% list("k-nn", "k-medoids")) {
    stop("Error: Not a valid clustering algorithm. Must be 'k-nn' or 'k-medoids'.")
  }

  if (save_graph) {
    if (!endsWith(graph_file_name, ".png")) {
      stop("Error: The graph file name must end with '.png'. Only PNG images are supported at this time.")
    }
  }

  # Using cached information if present, otherwise calculating it and storing
  #   it as an environment variable
  stopifnot("Aggregated_Dataframe" %in% class(df_aggregated))

  message("Clustering... \n \n")

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

  if (max_k == dim(distance_matrix)[1]) {

    message("\n \n Specified max_k value is too large for the size of the data, i.e. greater than the number of observations - 1 ... \n")
    message(cat("Set max_k to ", max_k - 1, " and run again ... \n \n"))

    stop()

  }





  calc <- function(n) {

    # Preallocating vectors
    k <- numeric(n)
    n_clusters <- numeric(n)
    cluster_size <- character(n)
    average_between <- numeric(n)
    average_within <- numeric(n)
    within_cluster_ss <- numeric(n)
    dunn <- numeric(n)
    wb_ratio <- numeric(n)
    #sum_diss <- numeric(n)

    average_silhouette_width_lower_ci <- numeric(n)
    average_silhouette_width <- numeric(n)
    average_silhouette_width_upper_ci <- numeric(n)
    

    for (i in min_k:n) {


      if (clustering == "k-nn") {

        algo = "K-Nearest Neighbors"
        clustered <- df_aggregated %>% cluster_knn(k = i)

      } else if (clustering == "k-medoids") {

        algo = "K-Medoids"
        clustered <- df_aggregated %>% cluster_kmedoids(k = i)

      }

      cluster_id <- clustered %>% unnest(df_sequences) %>% select(id, cluster) %>% arrange(id)

      if (max(clustered$cluster) == 1) {

        k[i] <- i
        n_clusters[i] <- 1
        cluster_size[i] <- toString(clustered$n)
        average_between[i] <- NA
        average_within[i] <- NA
        within_cluster_ss[i] <- NA
        dunn[i] <- NA
        wb_ratio[i] <- NA

        # Calculating information for silhouette
        sil <- NA

        nobs <- NA
        x_bar <- NA
        error <- NA

        average_silhouette_width_lower_ci[i] <- NA
        average_silhouette_width[i] <- NA
        average_silhouette_width_upper_ci[i] <- NA

      } else {

        clust_stats <- cluster.stats(d = distance_matrix, clustering = cluster_id$cluster)

        k[i] <- i
        n_clusters[i] <- clust_stats$cluster.number
        cluster_size[i] <- toString(clust_stats$cluster.size)
        average_between[i] <- clust_stats$average.between
        average_within[i] <- clust_stats$average.within
        within_cluster_ss[i] <- clust_stats$within.cluster.ss
        dunn[i] <- clust_stats$dunn
        wb_ratio[i] <- clust_stats$wb.ratio

        # Calculating information for silhouette
        sil <- silhouette(cluster_id$cluster, dmatrix = distance_matrix)

        nobs <- dim(sil)[1]
        x_bar <- mean(sil[,3])
        error <- qt(0.975, df = nobs - 1) * sd(sil[, 3])/sqrt(nobs)

        average_silhouette_width_lower_ci[i] <- x_bar - error
        average_silhouette_width[i] <- x_bar
        average_silhouette_width_upper_ci[i] <- x_bar + error

      }

    }

    slice(data.frame(k, n_clusters, cluster_size,
                     average_silhouette_width_lower_ci,
                     average_silhouette_width,
                     average_silhouette_width_upper_ci,
                     dunn,
                     average_within, average_between, wb_ratio,
                     within_cluster_ss), 2:n())


  }

  k_info <- calc(max_k)





  # Plotting Information
  if (validation_measure == 'silhouette') {

      measure = 'Average Silhouette Width'
      measure_values = k_info$average_silhouette_width

  } else if (validation_measure == 'dunn') {

      measure = 'Dunn Index'
      measure_values = k_info$dunn

  } else if (validation_measure == 'wb_ratio') {

      measure = 'Average Distance Within Cluster / Average Distance Between Clusters'
      measure_values = k_info$wb_ratio

  } else if (validation_measure == 'average_between') {

      measure = 'Average Distance Between Clusters'
      measure_values = k_info$average_between

  } else if (validation_measure == 'average_within') {

      measure = 'Average Distance Within Cluster'
      measure_values = k_info$average_within

  } else if (validation_measure == 'within_cluster_ss') {

      measure = 'Sum of Within Cluster / Cluster Size'
      measure_values = k_info$within_cluster_ss

  } else {

    stop("Only validation measures of silhouette, dunn, wb_ratio, average_between, average_within, and within_cluster_ss are supported.")

  }


  if (validation_measure == 'silhouette') {

   k_plot <- ggplot(k_info, aes(k, average_silhouette_width)) +

      geom_line(color = "#20B2AA") +

      geom_errorbar(aes(ymax = average_silhouette_width_upper_ci,
                        ymin = average_silhouette_width_lower_ci),
                    width = .25,
                    color = "#20B2AA") +

      geom_vline(xintercept = k_info$k[which.max(k_info$average_silhouette_width)],
                 color = "#20B2AA", linetype = 'dashed') +

      labs(title = "Optimal K Plot",
           subtitle = paste0("k =", k_info$k[which.max(k_info$average_silhouette_width)], "; Max average silhouette width = ", round(k_info$average_silhouette_width[which.max(k_info$average_silhouette_width)], digits = 3)),
           x = "K Value",
           y = measure)  +


     coord_cartesian(xlim = c(min(k_info$k), max(k_info$k))) +

     scale_x_continuous(labels = as.character(k_info$k), breaks = k_info$k)



  } else {


    k_plot <- ggplot(k_info, aes(k, measure_values)) +

      geom_line(color = "#20B2AA") +

      labs(title = "Optimal K Plot",
           x = "K Value",
           y = measure) +

      coord_cartesian(xlim = c(min(k_info$k), max(k_info$k))) +

      scale_x_continuous(labels = as.character(k_info$k), breaks = k_info$k)



  }





  # This portion saves the graph if the option is selected
  if (save_graph) {
    output_directory <- create_folder(output_directory, "approxmap_results")
    output_directory_graphs <- create_folder(output_directory, "graphs")

    png(file = paste0(output_directory_graphs, "/", file_check(output_directory_graphs,
               graph_file_name)), width = size_width, height = size_height)
    print(k_plot)
    dev.off()


  }


  print(k_plot)
  return(k_info)

}
