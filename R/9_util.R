# pull.list <- function(df, column){
#   column <- quo_name(column)
#   df[,column] %>%
#     unlist()
# }

format_sequence.W_Sequence_Dataframe <- function(df){
  columns <- c(names(df)[str_detect(names(df),"_pattern")], "weighted_sequence")

  df %>%
    select(one_of("cluster", "n", columns)) %>%
    mutate_at(columns, format_sequence)
}

compare_sequences <- function(df){
  df %>%
    gather(-cluster,-n, key = "pattern", value = "w_sequence") %>%
      arrange(cluster)
}
