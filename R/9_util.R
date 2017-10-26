format_sequence.W_Sequence_Dataframe <- function(df,
                                                 compare = FALSE,
                                                 truncate_patterns = FALSE,
                                                 html_format = FALSE){

  column_patterns <- names(df)[str_detect(names(df),"_pattern")]
  columns <- c(column_patterns, "weighted_sequence")

  if(truncate_patterns){
    df <-
      df %>%
      mutate_at(column_patterns, truncate_pattern)
  }



  df <-
    df %>%
      select(one_of("cluster", "n", columns)) %>%
      mutate_at(columns, function(x) format_sequence(x, html_format = html_format)) %>%
      mutate(n = as.double(n),
             n_percent = str_c(round(n/sum(n) * 100, digits = 2),"%")) %>%
    select(one_of("cluster", "n", "n_percent", columns))

   if(compare){
    compare_sequences(df)
  } else {
    df
  }

}

compare_sequences <- function(df){
  df %>%
    gather(-cluster, -n, -n_percent, key = "pattern", value = "sequence") %>%
      arrange(cluster) %>%
      mutate(pattern = stringr::str_replace(pattern, "_pattern", ""))
}



generate_summary_stats <- function(input_data, results_directory = "~", noise_threshold = 0, write_files = FALSE){

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

class_it <- function(obj, class_name){
  class(obj) <- c(class_name, class(obj)) %>%
    unique()
  obj
}

truncate_pattern <- function(x, ...){
  UseMethod("truncate_pattern")
}

truncate_pattern.W_Sequence_Pattern_List <- function(w_sequence_list){
  class_it(map(w_sequence_list, truncate_pattern), "W_Sequence_List")
}

truncate_pattern.W_Sequence_Pattern <- function(w_sequence){
  truncate_index <- rep(FALSE, length(w_sequence))
  for(i in seq_along(w_sequence)){
    if(i == length(w_sequence)) break()
    e_1 <- sort(w_sequence[[i]]$elements)
    e_2 <- sort(w_sequence[[i + 1]]$elements)
    if(identical(e_1, e_2)){
      truncate_index[i] <- TRUE
    }
  }
  w_sequence[truncate_index] <- NULL
  w_sequence
}

