format_sequence.W_Sequence_Dataframe <- function(df,
                                                 compare = FALSE,
                                                 truncate_patterns = FALSE,
                                                 html_format = FALSE){

  column_patterns <- names(df)[str_detect(names(df),"_pattern")]
  if("weighted_sequence" %in% names(df)){
    columns <- c(column_patterns, "weighted_sequence")
  } else {
    columns <- column_patterns
  }

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

w_sequence_to_tibble <- function(w_sequence){
  tibble(element = map(w_sequence, "elements") %>% unlist(),
         element_weight = map(w_sequence, "element_weights") %>% unlist(),
         itemset = map2(1:length(w_sequence), w_sequence, ~rep(.x ,length(.y$elements))) %>% unlist()) %>%
  mutate(element_no = row_number())
}

plot_weighted_sequence <- function(w_sequence){
  df_sequence <-
    w_sequence %>%
    w_sequence_to_tibble()

  df_itemset <-
  df_sequence %>%
      group_by(itemset) %>%
      filter(element_no == max(element_no))

  df_sequence %>%
    ggplot(aes(element_no, element_weight)) +
    geom_point() +
    geom_label(aes(y = element_weight + 0.02*element_weight, label = element)) +
    geom_vline(data = df_itemset, aes(xintercept = element_no))
}

