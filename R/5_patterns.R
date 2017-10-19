filter_pattern <- function(x, ...){
  UseMethod("filter_pattern")
}

filter_pattern.Clustered_Dataframe <- function(df_clustered, ...){
  df_clustered %>%
    get_weighted_sequence() %>%
      filter_pattern(...)
}


filter_pattern.W_Sequence_Dataframe <- function(df_w_sequences, pattern_name = "consensus", ...){

  if(!("weighted_sequence" %in% names(df_w_sequences))){
      stop("Please use a dataframe where weighted sequences can be calculated")
  }


  # expr <- enquo(pattern_name)
  pattern <- paste0(quo_name(pattern_name),"_pattern")

  df_pattern <-
    df_w_sequences %>%
    mutate(!!pattern := structure(map(weighted_sequence, filter_pattern, ... = ...),
                                  class = "W_Sequence_List")) %>%
    select(cluster, n, !!pattern, everything())

  class(df_pattern) <- class(df_w_sequences)

  df_pattern
}


filter_pattern.W_Sequence <- function(w_sequence,
                                      threshold = 0.5,
                                      noise_threshold = 0,
                                      blank_if_absent = F,
                                      pattern_name = NULL){
  n <- attr(w_sequence, "n")

  elements <- unlist(map(w_sequence,
                         function(x) {
                          if(length(x$elements)==0) return(NULL)
                          x$elements
                         }
                         )
                     )
  if(length(elements) == 0) return(weighted_seq)
  min_occurences <- n * threshold
  min_occurences <- max(min_occurences, noise_threshold)



  pattern <-
    map(w_sequence, function(w_sequence_itemset){

      threshold_check <- w_sequence_itemset$element_weights >= min_occurences

      if(blank_if_absent && (sum(threshold_check) == 0)) {
        w_sequence_itemset$elements <- ""
        w_sequence_itemset$element_weights <- 0
      } else if(sum(threshold_check) == 0){
        w_sequence_itemset$elements <- NULL
        w_sequence_itemset$element_weights <- NULL
        w_sequence_itemset$itemset_weight <- NULL
      } else {
         w_sequence_itemset$elements <- w_sequence_itemset$elements[threshold_check]
         w_sequence_itemset$element_weights <- w_sequence_itemset$element_weights[threshold_check]
      }

      w_sequence_itemset
    })

  class(pattern) <- c("W_Sequence",  class(pattern))

  if(!is.null(pattern_name)){
    class(pattern) <- c(paste0(pattern_name,"_Pattern"),  class(pattern))
  }

  attr(pattern, "n") <- n

  pattern

}



