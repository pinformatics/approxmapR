#' @export
create_folder <- function(directory, folder) {
  output_directory = paste0(directory, "/", folder)

  if (!dir.exists(output_directory)) {
    dir.create(output_directory)
  }
  output_directory
}

#' @export
file_check <- function(dir = ".", file_name) {
  files_exist <- list.files(path = dir)
  file_part <- str_replace(file_name, "\\.[A-Za-z]*", "")
  extension <- str_replace(file_name, file_part, "")

  files_exist <- files_exist[str_detect(files_exist, file_part)]
  if (length(files_exist) == 0) {
    file_name
  } else {
    int <- str_replace(files_exist, file_part, "") %>%
      str_replace("_", "") %>%
      str_replace(extension, "") %>%
      as.integer()

    int <- int[!is.na(int)]
    if (length(int) == 0) {
      int <- 0
    }

    if (sum(is.na(int)) > 0) {
      paste0(file_part, "_1", extension)
    } else {
      paste0(file_part, "_", max(int) + 1, extension)
    }
  }
}

#' @export
generate_reports <- function(w_sequence_dataframe,
                             sil_table = NULL,
                             html_format = TRUE,
                             # truncate_patterns = FALSE,
                             output_directory = "~",
                             end_filename_with = "",
                             sequence_analysis_details = NULL,
                             sequence_analysis_details_definitions = NULL,
                             algorithm_comparison = FALSE) {
  stopifnot("W_Sequence_Dataframe" %in% class(w_sequence_dataframe))

  if (!is.null(sil_table)) {

    if (!identical(names(sil_table), c("id", "cluster", "neighbor", "sil_width"))) {

      stop("Error: Columns must be in the following order 'id', 'cluster', 'neighbor', 'sil_width'")

    }
  }

  if (!is.null(sequence_analysis_details)) {
    if (!"list" %in% class(sequence_analysis_details)) {
      stop("Error: sequence_analysis_details must be a list class.")
    }

    if (!"algorithm" %in% names(sequence_analysis_details))  {
      stop("Error: Missing the algorithm parameter in sequence_analysis_details. Must have algorithm, k_value, cluster_n, time_period, consensus_threshold, and notes parameters. See ??generate_reports() for help.")
    }
    if (!"k_value" %in% names(sequence_analysis_details)) {
      stop("Error: Missing the k_value parameter in sequence_analysis_details. Must have algorithm, k_value, cluster_n, time_period, consensus_threshold, and notes parameters. See ??generate_reports() for help.")
    }
    if (!"cluster_n" %in% names(sequence_analysis_details)) {
      stop("Error: Missing the cluster_n parameter in sequence_analysis_details. Must have algorithm, k_value, cluster_n, time_period, consensus_threshold, and notes parameters. See ??generate_reports() for help.")
    }
    if (!"time_period" %in% names(sequence_analysis_details)) {
      stop("Error: Missing the time_period parameter in sequence_analysis_details.Must have algorithm, k_value, cluster_n, time_period, consensus_threshold, and notes parameters. See ??generate_reports() for help.")
    }
    if (!"consensus_threshold" %in% names(sequence_analysis_details)) {
      stop("Error: Missing the consensus_threshold parameter in sequence_analysis_details. Must have algorithm, k_value, cluster_n, time_period, consensus_threshold, and notes parameters. See ??generate_reports() for help.")
    }
    if (!"notes" %in% names(sequence_analysis_details)) {
      stop("Error: Missing the notes parameter in sequence_analysis_details. Must have algorithm, k_value, cluster_n, time_period, consensus_threshold, and notes parameters. See ??generate_reports() for help.")
    }
  }

  folder = "approxmap_results"
  output_directory <- create_folder(output_directory, folder)
  output_directory_public <-
    create_folder(output_directory, "public")
  output_directory_private <-
    create_folder(output_directory, "private")


  # Checking which report file structure to be used
  if (!is.null(sequence_analysis_details) & is.null(sequence_analysis_details_definitions)) {

    report_rmd <- system.file("rmd_w_sequence_analysis_details.Rmd", package = "approxmapR")

  } else if (!is.null(sequence_analysis_details) & !is.null(sequence_analysis_details_definitions)) {

    report_rmd <- system.file("rmd_w_sequence_analysis_details_definitions.Rmd", package = "approxmapR")

  } else {

    report_rmd <- system.file("rmd_w_sequence.Rmd", package = "approxmapR")

  }



  formatted <-
    w_sequence_dataframe %>%
    format_sequence(
      compare = TRUE,
      html_format = html_format,
      truncate_patterns = FALSE
    )

  formatted_trunc <-
    w_sequence_dataframe %>%
    select(-weighted_sequence) %>%
    format_sequence(
      compare = TRUE,
      html_format = html_format,
      truncate_patterns = TRUE
    ) %>%
    # filter(pattern != "weighted_sequence") %>%
    mutate(pattern = pattern %>% str_c("_truncated"))

  df_unique_items <-
    w_sequence_dataframe %>%
    mutate(
      pattern = "unique_items",
      sequence =
        weighted_sequence %>%
        map_chr(function(pattern) {
          pattern %>%
            map("elements") %>%
            unlist() %>%
            unique() %>%
            str_c(collapse = ", ")
        }),
      n = as.double(n),
      n_percent = str_c(round(n / sum(n) * 100, digits = 2), "%")
    ) %>%
    select(-ends_with("_pattern"), -weighted_sequence, -df_sequences)


  if (!is.null(sequence_analysis_details_definitions)) {
    event_defs <- subset(sequence_analysis_details_definitions, event %in% unique(str_split(df_unique_items$sequence, ", ") %>% unlist()))
  }


  formatted <-
    formatted %>%
    bind_rows(formatted_trunc) %>%
    arrange(cluster, pattern) %>%
    bind_rows(df_unique_items) %>%
    arrange(cluster)

  if (!is.null(sequence_analysis_details) & is.null(sequence_analysis_details_definitions)) {
    rmarkdown::render(
      report_rmd,
      params = append(list(input = formatted,
                           title = "All Sequences"), sequence_analysis_details),
                    #output_file = file_check(output_directory_private, "all_sequences.html"),
                    output_file = file_check(output_directory_private, paste0("all_sequences", end_filename_with, ".html")),
                    output_dir = output_directory_private
                  )

  } else if (!is.null(sequence_analysis_details) & !is.null(sequence_analysis_details_definitions)) {

    rmarkdown::render(
      report_rmd,
      params = append(list(input = formatted,
                           title = "All Sequences",
                           event_definitions = event_defs), sequence_analysis_details),
                    #output_file = file_check(output_directory_private, "all_sequences.html"),
                    output_file = file_check(output_directory_private, paste0("all_sequences", end_filename_with, ".html")),
                    output_dir = output_directory_private
                  )
  } else {
    rmarkdown::render(
      report_rmd,
      params = list(input = formatted,
                    title = "All Sequences"),
                    #output_file = file_check(output_directory_private, "all_sequences.html"),
                    output_file = file_check(output_directory_private, paste0("all_sequences", end_filename_with, ".html")),
                    output_dir = output_directory_private
                  )
  }


  patterns <-
    formatted %>%
    filter(pattern != "weighted_sequence")

  if (!is.null(sequence_analysis_details)) {
    rmarkdown::render(
      report_rmd,
      params = append(list(input = patterns,
                           title = "Patterns"), sequence_analysis_details),
                    #output_file = file_check(output_directory_public, "patterns.html"),
                    output_file = file_check(output_directory_public, paste0("patterns", end_filename_with, ".html")),
                    output_dir = output_directory_public
                  )
 } else {
   rmarkdown::render(
     report_rmd,
     params = list(input = patterns,
                   title = "Patterns"),
                   #output_file = file_check(output_directory_public, "patterns.html"),
                   output_file = file_check(output_directory_public, paste0("patterns", end_filename_with, ".html")),
                   output_dir = output_directory_public
                 )
 }


  w_sequences <-
    formatted %>%
    filter(pattern == "weighted_sequence") %>%
    select(-pattern)

  if (!is.null(sequence_analysis_details)) {
    rmarkdown::render(
      report_rmd,
      params = append(list(input = w_sequences,
                           title = "Weighted Sequences"), sequence_analysis_details),
                    #output_file = file_check(output_directory_private, "weighted_sequences.html"),
                    output_file = file_check(output_directory_private, paste0("weighted_sequences", end_filename_with, ".html")),
                    output_dir = output_directory_private
                  )
 } else {
   rmarkdown::render(
     report_rmd,
     params = list(input = w_sequences,
                   title = "Weighted Sequences"),
                   #output_file = file_check(output_directory_private, "weighted_sequences.html"),
                   output_file = file_check(output_directory_private, paste0("weighted_sequences", end_filename_with, ".html")),
                   output_dir = output_directory_private
                 )
 }

  # output_directory_alignments <- create_folder(output_directory_private, "alignments")

  message("saving alignments...")

  if (!is.null(sil_table)) {

    sil_table <- sil_table %>% group_by(cluster) %>% mutate(cluster_average = mean(sil_width))
    sil_table$grand_avg_sil <- mean(sil_table$sil_width)
    names(sil_table) <- c("id2", "cluster", "neighbor", "sil_width",  "cluster_average", "grand_avg_sil")

  }


  if (!is.null(sil_table)) {

    paste0("Grand Silhouette Width", ", ", sil_table$grand_avg_sil[[1]], "\n", "\n",
           "Cluster", ", ", "Cluster Silhoutte Width", ", ", "Neighbor Cluster", ", ", "id", ", ", "id Silhouette Width", ", ", "Sequence", "\n",
           w_sequence_dataframe %>% save_alignment(save_date = FALSE, sil_table = sil_table, algorithm_comparison = algorithm_comparison)
          ) %>%
      write_file(paste0(
        output_directory_private,
        "/",
        #file_check(output_directory_private, "alignments_with_date.csv")
        file_check(output_directory_private, paste0("alignments", end_filename_with, ".csv"))
      ))

  } else {

    w_sequence_dataframe %>%
      save_alignment(save_date = FALSE, sil_table =sil_table, algorithm_comparison = algorithm_comparison) %>%
      write_file(paste0(
        output_directory_private,
        "/",
        #file_check(output_directory_private, "alignments.csv")
        file_check(output_directory_private, paste0("alignments", end_filename_with, ".csv"))
      ))

  }


  if (!is.null(sil_table)) {

    paste0("Grand Silhouette Width", ", ", sil_table$grand_avg_sil[[1]], "\n", "\n",
           "Cluster", ", ", "Cluster Silhoutte Width", ", ", "Neighbor Cluster", ", ", "id", ", ", "id Silhouette Width", ", ", "Sequence", "\n",
           w_sequence_dataframe %>% save_alignment(save_date = TRUE, sil_table = sil_table, algorithm_comparison = algorithm_comparison)
          ) %>%
      write_file(paste0(
        output_directory_private,
        "/",
        #file_check(output_directory_private, "alignments_with_date.csv")
        file_check(output_directory_private, paste0("alignments_with_date", end_filename_with, ".csv"))
      ))

  } else {

    w_sequence_dataframe %>%
      save_alignment(save_date = TRUE, sil_table = sil_table, algorithm_comparison = algorithm_comparison) %>%
      write_file(paste0(
        output_directory_private,
        "/",
        #file_check(output_directory_private, "alignments_with_date.csv")
        file_check(output_directory_private, paste0("alignments_with_date", end_filename_with, ".csv"))
      ))

  }






  message("Reports generated.")

}

#' @export
generate_summary_stats <- function(input_data,
                                   results_directory = "~",
                                   noise_threshold = 0,
                                   write_files = FALSE) {
  input_data <- as_tibble(input_data) %>% ungroup()

  n_seq <- input_data %>%
    pull(id) %>%
    unique() %>%
    length()

  cat(noquote(paste0("The number of sequences is ", n_seq)))

  n_unique_items <- input_data %>%
    pull(event) %>% unique() %>%
    length()

  cat(noquote(
    sprintf("\n\nThe number of unique items is %i\n", n_unique_items)
  ))

  items <-
    input_data %>%
    count(event, sort = T) %>%
    top_n(20, n) %>%
    filter(n > 5) %>%
    mutate(n = n / max(n)) %>%
    rename(relative_freq = n, item = event)
  # pull(event)

  # cat(noquote("The unique items are \n"))
  # cat(str_c(items, collapse = ", "))
  print(items)

  cat(noquote("\nStatistics for the number of sets per sequence:\n"))
  n_sets <-
    input_data %>%
    select(id, period) %>%
    group_by(id) %>% unique() %>%
    summarise(len = n()) %>%
    pull(len)
  print(summary(n_sets))


  cat(noquote("\nStatistics for the number of items in a set:\n"))
  n_items_in_set <-
    input_data %>%
    group_by(id, period) %>%
    unique() %>%
    summarise(len = n()) %>%
    pull(len)
  print(summary(n_items_in_set))

  cat(noquote("\nFrequencies of items:\n"))
  count_items <-
    input_data %>%
    count(event) %>%
    arrange(desc(n)) %>%
    filter(n > noise_threshold)

  print(summary(count_items$n))
  cat("\n")

  if (write_files) {
    results_directory <- paste0(results_directory, "/public")
    if (!dir.exists(results_directory))
      dir.create(results_directory)
    write_csv(count_items,
              #paste0(results_directory, "/count_items.csv")
              paste0(results_directory, "/count_items", end_filename_with, ".csv")
            )
  }
}

#' @export
generate_summary_stats_dup <- function(input_data,
                                       results_directory = "~",
                                       noise_threshold = 0,
                                       write_files = FALSE) {
  input_data <- as_tibble(input_data) %>% ungroup()

  n_unique_items <- input_data %>%
    pull(event) %>% unique() %>%
    length()

  n_sets <-
    input_data %>%
    select(id, period) %>%
    group_by(id) %>% unique() %>%
    summarise(len = n()) %>%
    pull(len)

  n_sets_info <-
    summary(n_sets) %>%
    broom::tidy() %>%
    mutate(type = "set")


  n_items_in_set <-
    input_data %>%
    group_by(id, period) %>%
    unique() %>%
    summarise(len = n()) %>%
    pull(len)

  n_itemset_info <-
    summary(n_items_in_set) %>%
    broom::tidy() %>%
    mutate(type = "set_items")


  count_items <-
    input_data %>%
    count(event) %>%
    arrange(desc(n)) %>%
    filter(n > noise_threshold)

  count_items_info <-
    summary(count_items$n) %>%
    broom::tidy() %>%
    mutate(type = "count_items")

  bind_rows(n_sets_info, n_itemset_info, count_items_info)
}


#' @export
print_alignments <- function(df) {
  wseqs <- df$weighted_sequence
  # browser()
  wseqs %>%
    walk2(df$cluster, function(x, i) {
      # browser()
      cat("\n\n")
      cat(glue("Cluster {i}"))
      cat("\n")
      print(attr(x, "alignments"))
    })
}

align_date_to_seq <- function(current_id, seq){

  #browser()

  id_elements <- unlist(seq)

  id_dates <-
    env_dm$df_aggregated %>%
    filter(id == current_id) %>%
    pull(date)

  id_period <-
    env_dm$df_aggregated %>%
    filter(id == current_id) %>%
    pull(period)



  new_dates <- rep("_", length(id_elements))
  new_date_index <- 1

  new_period <- rep("_", length(id_elements))
  new_period_index <- 1


  for (i in seq(1, length(id_elements))) {
    if (id_elements[i] != "_"){
      new_dates[i] <- id_dates[new_date_index] %>% as.character()
      new_date_index = new_date_index + 1

      new_period[i] <- id_period[new_period_index] %>% as.character()
      new_period_index = new_period_index + 1
    }
  }



  newer_dates <- NULL
  newer_dates_index <- 1

  for (i in seq(1, length(new_period))) {

    if (i == 1) {

      if (new_period[i] == "_") {

        val = "_"

      } else {

        val <- new_dates[i]

      }

    } else {

      if (new_period[i] == new_period[i - 1]) {

        if (new_period[i] == "_") {

          newer_dates <- c(newer_dates, val)
          newer_dates_index = newer_dates_index + 1

          val <- new_dates[i]

        } else {

          val <- paste0(val, "; ", new_dates[i])

        }

      } else {

        newer_dates <- c(newer_dates, val)
        newer_dates_index = newer_dates_index + 1

        val <- new_dates[i]
      }

    }

  }
  newer_dates <- c(newer_dates, val)

  newer_dates %>% str_c(collapse = ", ")

}



#' @export
save_alignment <- function(x, ...) {
  UseMethod("save_alignment")
}


save_alignment.Sequence <- function(sequence) {
  sequence %>%
    map_chr(function(x) {
      str_c('"', paste0(x, collapse = "; "), '"')
    }) %>%
    paste0(collapse = ", ")
}

save_alignment.Sequence_List <- function(alignment, save_date=TRUE, sil_table = NULL, algorithm_comparison = algorithm_comparison) {
  map2_chr(alignment, names(alignment), function(seq, id) {

    if (algorithm_comparison) {

      temp <- id %>% str_split("_", simplify = TRUE)

      seqs <- str_c(temp[1], ", ", temp[2], ", ", temp[3], ", ",save_alignment(seq), "\n")

    } else {

        if (!is.null(sil_table)) {

          seqs <- str_c((filter(sil_table, id2 == id))$cluster, ", ",
                        (filter(sil_table, id2 == id))$cluster_average, ", ",
                        (filter(sil_table, id2 == id))$neighbor, ", ",
                        id, ", ",
                        (filter(sil_table, id2 == id))$sil_width, ", ",
                        save_alignment(seq), "\n")

        } else {

          seqs <- str_c(id, ", ", save_alignment(seq), "\n")

        }

    }

    if (save_date == TRUE) {

      if (!is.null(sil_table)) {

        dates <- str_c((filter(sil_table, id2 == id))$cluster, ", ",
                      (filter(sil_table, id2 == id))$cluster_average, ", ",
                      (filter(sil_table, id2 == id))$neighbor, ", ",
                      id, ", ",
                      (filter(sil_table, id2 == id))$sil_width, ", ",
                      align_date_to_seq(id, seq), "\n")
        paste0(seqs, dates)

      } else {

        dates <- str_c(id, ", ", align_date_to_seq(id, seq), "\n")
        paste0(seqs, dates)

      }

    } else {

      seqs
    }

  }) %>%

    str_c(collapse = "")
}




save_alignment.W_Sequence_Dataframe <- function(df, sil_table, ...) {

  if (!is.null(sil_table)) {

    map2_chr(df$cluster, df$weighted_sequence, function(c, seq) {

        str_c(save_alignment(attr(seq, "alignments"), sil_table = sil_table, ...))
            }) %>% str_c(collapse = "")

  } else {

    map2_chr(df$cluster, df$weighted_sequence, function(c, seq) {
        str_c("Cluster ", c, "\n", save_alignment(attr(seq, "alignments"), ...))
        }) %>% str_c(collapse = "\n")
  }
}



save_alignment.W_Sequence_Dataframe2 <- function(df, sil_table, ...) {

  if (!is.null(sil_table)) {

    map2_chr(df$cluster, df$weighted_sequence, function(c, seq) {


        str_c("Grand Silhouette Width", ", ", sil_table$grand_avg_sil[[1]], "\n", "\n",
              "Cluster", ", ", "Cluster Silhoutte Width", ", ", "Neighbor Cluster", ", ", "id", ", ", "id Silhouette Width", ", ", "Sequence", "\n",
              save_alignment(attr(seq, "alignments"), sil_table = sil_table, ...))
            }) %>% str_c(collapse = "")

  } else {

    map2_chr(df$cluster, df$weighted_sequence, function(c, seq) {
        str_c("Cluster ", c, "\n", save_alignment(attr(seq, "alignments"), ...))
        }) %>% str_c(collapse = "\n")
  }
}



#' @export
algorithm_comparison <- function(formatted1, formatted1_pars = "No Pars1",
                                 formatted2, formatted2_pars = "No Pars2",
                                 formatted3, formatted3_pars = "No Pars2",
                                 approxmapR_results = "~",
                                 output_directory = "~",
                                 file_name = "algorithm-comparison") {

  report_rmd <- system.file("algo_compar.Rmd", package = "approxmapR")

  rmarkdown::render(

    report_rmd,

    params = list(formatted1 = formatted1,
                  formatted1_pars = formatted1_pars,

                  formatted2 = formatted2,
                  formatted2_pars = formatted2_pars,

                  formatted3 = formatted3,
                  formatted3_pars = formatted3_pars,

                  approxmapR_results = approxmapR_results),

    output_file = file_check(output_directory, paste0(file_name, ".html")),
    output_dir = output_directory


  )
}
