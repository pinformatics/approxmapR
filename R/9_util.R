format_sequence.W_Sequence_Dataframe <- function(df, compare = FALSE){
  columns <- c(names(df)[str_detect(names(df),"_pattern")], "weighted_sequence")

  df <-
    df %>%
      select(one_of("cluster", "n", columns)) %>%
      mutate_at(columns, format_sequence)

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
