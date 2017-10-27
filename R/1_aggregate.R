get_n_days <- function(unit,n_units) {
  if(unit == "day") {
    n_days <- n_units
  } else if(unit == "week"){
    n_days <- n_units * 7
  } else if(unit == "month"){
    n_days <- n_units * 30
  } else if(unit == "6 months"){
    n_days <- n_units * 180
  } else if(unit == "year"){
    n_days <- n_units * 365
  }
  lubridate::ddays(n_days)
}

pre_aggregated <- function(df, summary_stats = TRUE){
  if(!(all(names(df) == c("id", "period", "event")))) {
    stop("There should be 3 columns named id, period and event (all lower case).")
  }

  stopifnot(is.integer(df$period))
  # if(!){
  #   stop("period should be numeric")
  # }

  if(!is.character(df$event)) df$event = as.character(df$event)


  class(df) <- c("Aggregated_Dataframe", "tbl_df", "tbl", "data.frame")
  if(summary_stats){
    message("Generating summary statistics of aggregated data...")
    generate_summary_stats(aggregated_data)
  }

  # .GlobalEnv$env_report <- new.env()
  # .GlobalEnv$env_report$aggregated_data <- aggregated_data

  df
}

aggregate_sequences <- function(unaggregated_data,
                     format = "%m-%d-%Y",
                     calendar = FALSE,
                     unit = "week",
                     n_units = 1,
                     anchor_table = NA,
                     anchor_vector = NA,
                     base_date = NA,
                     occurence = min,
                     multiset = FALSE,
                     include_date = FALSE,
                     summary_stats = TRUE) {

  n_days <- get_n_days(unit,n_units)

  names(unaggregated_data) <- c("id","date","event")

  unaggregated_data2 <-
    unaggregated_data %>%
    mutate(date = readr::parse_date(date,format))

  if(!is.na(anchor_table) || !is.na(anchor_vector)){

    if(is.data.frame(anchor_table)){

      names(anchor_table) <- c("id","anchor_date")

      anchor_table <-
        anchor_table %>%
        mutate(anchor_date = readr::parse_date(anchor_date,format))

    } else if(!is.na(anchor_vector)) {
      anchor_table <-
        unaggregated_data2 %>%
        filter(event %in% anchor_vector) %>%
        group_by(id) %>%
        summarise(anchor_date = min(date))
    }

    aggregated_data <-
      unaggregated_data2 %>%
      left_join(anchor_table, by = "id") %>%
      group_by(id) %>%
      mutate(n_ndays = (date - anchor_date)/n_days,
             agg_n_ndays = if_else(n_ndays < 0, floor(n_ndays), ceiling(n_ndays)),
             agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)) %>%
      arrange(id, agg_n_ndays) %>%
      mutate(agg_period =  dense_rank(agg_n_ndays))

  } else if(calendar){

    aggregated_data <-
      unaggregated_data2 %>%
      mutate(agg_date = lubridate::floor_date(date, paste0(n_units, " ", unit, "s"))) %>%
      group_by(id) %>%
      mutate(agg_period =  dense_rank(agg_date))

  } else if(!is.na(base_date)){

    if(typeof(base_date) == "closure") {
      base_date <- base_date(unaggregated_data2$date)
    }

    aggregated_data <-
      unaggregated_data2 %>%
      mutate(n_ndays = (date - base_date)/n_days,
             agg_n_ndays = if_else(n_ndays < 0, floor(n_ndays), ceiling(n_ndays)),
             agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)) %>%
      group_by(id) %>%
      arrange(id, agg_n_ndays) %>%
      mutate(agg_period =  dense_rank(agg_n_ndays))

  } else {
    aggregated_data <-
      unaggregated_data2 %>%
      group_by(id) %>%
      mutate(n_ndays = (date - occurence(date))/n_days,
             agg_n_ndays = if_else(n_ndays < 0, floor(n_ndays), ceiling(n_ndays)),
             agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)) %>%
      arrange(id, agg_n_ndays) %>%
      mutate(agg_period =  dense_rank(agg_n_ndays))
  }



  if(!multiset){
    aggregated_data <-
      aggregated_data %>%
      group_by(id,agg_period,event) %>%
      slice(1)
  }



  if(include_date){
    aggregated_data <-
      aggregated_data %>%
      select(id, date, period = agg_period, event) %>%
      arrange(id, period, date ,event)
  } else{
    aggregated_data <-
      aggregated_data %>%
      select(id, period = agg_period, event) %>%
      arrange(id, period, event)
  }

  class(aggregated_data) <- c("Aggregated_Dataframe", class(aggregated_data))
  if(summary_stats){
    message("Generating summary statistics of aggregated data...")
    generate_summary_stats(aggregated_data)
  }

  # .GlobalEnv$env_report <- new.env()
  # .GlobalEnv$env_report$aggregated_data <- aggregated_data

  aggregated_data
}




convert_to_sequence <- function(df_seq){
  if(!"Aggregated_Dataframe" %in% class(df_seq)){
    warning("Are you sure the sequence dataframe you passed is already aggregated?")
  }
  df_seq <-
    df_seq %>%
    group_by(id) %>%
    nest(.key = "nested_id") %>%
    mutate(sequence = map(nested_id,
                 function(df_id){
                   seqs <-
                     df_id %>%
                     group_by(period) %>%
                     nest(.key = "list_data") %>%
                     mutate(seqs = map(list_data,
                              function(itemset){
                                events <- itemset$event
                                class(events) <- c("Sequence_Itemset", class(events))
                                events
                              })) %>%
                     pull(seqs)
                   class(seqs) <- c("Sequence", class(seqs))
                   seqs
                 }),
           sequence_formatted = map_chr(sequence, format_sequence)
          ) %>%
    select(id, sequence, sequence_formatted)
  names(df_seq$sequence) <- df_seq$id
  class(df_seq$sequence) <- c("Sequence_List", class(df_seq$sequence))

  df_seq
}

format_sequence <- function(x, ...){
  UseMethod("format_sequence")
}

format_sequence.Sequence <- function(sequence){
  sequence <-
    sequence %>%
    map_chr(function(itemset){
      itemset <- str_c(itemset, collapse = ", ")
      paste0("(", itemset, ")")
    }) %>%
    str_c(collapse = " ")
  paste0("<", sequence, ">")
}

print.Sequence <- function(sequence){
  sequence %>%
  format_sequence() %>%
  cat()
}

print.Sequence_List <- function(sequences){
  # print(sequences)
  if(is.null(names(sequences))){
    warning("id for the sequences not present")
    walk(sequences, function(sequence_obj){
      cat(format_sequence(sequence_obj), "\n")
    })
  } else{
    walk2(sequences, names(sequences), function(sequence_obj, id){
      cat(id,": ", format_sequence(sequence_obj), "\n")
    })
  }
}

print_raw <- function(obj){
  obj %>%
    unclass() %>%
    print()
}

#' Pipe
#'
#' @importFrom magrittr %>%
#' @name
#' @rdname pipe
#' @export
NULL


# "select","filter","mutate","summarise","group_by","dense_rank"
