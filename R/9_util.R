format_sequence.W_Sequence_Dataframe <- function(df, compare = FALSE, n_as_percent = TRUE){
  columns <- c(names(df)[str_detect(names(df),"_pattern")], "weighted_sequence")

  df <-
    df %>%
      select(one_of("cluster", "n", columns)) %>%
      mutate_at(columns, format_sequence) %>%
      mutate(n = as.double(n),
             n = if_else(rep(n_as_percent, nrow(df)),
                         round(n/sum(n) * 100, digits = 2),
                         n)
             )


   if(compare){
    compare_sequences(df)
  } else {
    df
  }

}

compare_sequences <- function(df){
  df %>%
    gather(-cluster,-n, key = "pattern", value = "w_sequence") %>%
      arrange(cluster)
}



generate_summary_stats <- function(input_data, results_directory = "~", noise_threshold = 0, write_files = FALSE){

  input_data <- as_tibble(input_data) %>% ungroup()
  # message("Getting summary statistics...")

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

