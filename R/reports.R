create_folder <- function(directory, folder){
  output_directory = paste0(directory, "/", folder)

  if(!dir.exists(output_directory)){
    dir.create(output_directory)
  }
  output_directory
}

file_check <- function(dir = ".", file_name){

  files_exist <- list.files(path = dir)
  file_part <- str_replace(file_name, "\\.[A-Za-z]*","")
  extension <- str_replace(file_name, file_part,"")

  files_exist <- files_exist[str_detect(files_exist, file_part)]
  if(length(files_exist) == 0){
    file_name
  } else {
    int <- str_replace(files_exist, file_part, "") %>%
    str_replace("_","") %>%
    str_replace(extension, "") %>%
    as.integer()

    int <- int[!is.na(int)]

    if(sum(is.na(int)) > 0) {
      paste0(file_part, "_1", extension)
    } else {
      paste0(file_part, "_", max(int) + 1, extension)
    }
  }
}


generate_reports <- function(w_sequence_dataframe,
                             html_format = TRUE,
                             truncate_patterns = FALSE,
                             output_directory = "~",
                             folder = "approxmap_results"){
  stopifnot("W_Sequence_Dataframe" %in% class(w_sequence_dataframe))


  output_directory <- create_folder(output_directory, folder)
  output_directory_public <- create_folder(output_directory, "public")
  output_directory_private <- create_folder(output_directory, "private")

  report_rmd <- system.file("rmd_w_sequence.Rmd", package="approxmapR")



  formatted <-
    w_sequence_dataframe %>%
    format_sequence(compare = TRUE,
                    html_format = html_format,
                    truncate_patterns = truncate_patterns)

  rmarkdown::render(report_rmd,
                    params = list(input = formatted,
                                  title = "All Sequences"),
                    output_file = file_check(output_directory_private,"all_sequences.html"),
                    output_dir = output_directory_private)
  patterns <-
    formatted %>%
    filter(pattern != "weighted_sequence")

  rmarkdown::render(report_rmd,
                    params = list(input = patterns,
                                  title = "Patterns"),
                    output_file = file_check(output_directory_public, "patterns.html"),
                    output_dir = output_directory_public)


  w_sequences <-
    formatted %>%
    filter(pattern == "weighted_sequence") %>%
    select(-pattern)

  rmarkdown::render(report_rmd,
                    params = list(input = w_sequences,
                                  title = "Weighted Sequences"),
                    output_file = file_check(output_directory_private,"weighted_sequences.html"),
                    output_dir = output_directory_private)

}


generate_summary_stats <- function(input_data,
                                   results_directory = "~",
                                   noise_threshold = 0,
                                   write_files = FALSE){

  input_data <- as_tibble(input_data) %>% ungroup()

  n_unique_items <- input_data %>%
    pull(event) %>% unique() %>%
    length()
  cat(noquote(sprintf("The number of unique items is %i\n", n_unique_items)))


  cat(noquote("\nStatistics for the number of sets per sequence:\n"))
  n_sets <-
    input_data %>%
    select(id,period) %>%
    group_by(id) %>% unique() %>%
    summarise(len = n()) %>%
    pull(len)
  print(summary(n_sets))


  cat(noquote("\nStatistics for the number of items in a set\n"))
  n_items_in_set <-
    input_data %>%
    group_by(id,period) %>%
    unique() %>%
    summarise(len = n()) %>%
    pull(len)
  print(summary(n_items_in_set))

  cat(noquote("\nFrequencies of items\n"))
  count_items <-
    input_data %>%
    count(event) %>%
    arrange(desc(n)) %>%
    filter(n > noise_threshold)

  print(summary(count_items$n))
  cat("\n")

  if(write_files){
    results_directory <- paste0(results_directory,"/public")
    if(!dir.exists(results_directory)) dir.create(results_directory)
    write_csv(count_items,paste0(results_directory,"/count_items.csv"))
  }
}


generate_summary_stats_dup <- function(input_data,
                                   results_directory = "~",
                                   noise_threshold = 0,
                                   write_files = FALSE){

  input_data <- as_tibble(input_data) %>% ungroup()

  n_unique_items <- input_data %>%
    pull(event) %>% unique() %>%
    length()

  n_sets <-
    input_data %>%
    select(id,period) %>%
    group_by(id) %>% unique() %>%
    summarise(len = n()) %>%
    pull(len)

  n_sets_info <-
       summary(n_sets) %>%
       broom::tidy() %>%
       mutate(type = "set")


  n_items_in_set <-
    input_data %>%
    group_by(id,period) %>%
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

