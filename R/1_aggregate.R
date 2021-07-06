



#' Aggregation functions
#'
#' A function used by \code{\link{aggregate_sequences}} to get
#' the unit of aggregation
#'
#' @param unit The unit of aggregation. Takes one of c('day', 'week',
#' 'month', '6 months', 'year')
#' @param n_units The number of units to aggregate
#'
#' @return Number of days to aggregate in the input data
#'
#' @examples get_n_days('week', 2)
get_n_days <- function(unit, n_units) {
    if (unit == "day") {
        n_days <- n_units
    } else if (unit == "week") {
        n_days <- n_units * 7
    } else if (unit == "month") {
        n_days <- n_units * 30
    } else if (unit == "6 months") {
        n_days <- n_units * 180
    } else if (unit == "year") {
        n_days <- n_units * 365
    }
    lubridate::ddays(n_days)
}


#' Aggregation functions
#'
#' Used whenever the df to be analyzed in preaggregated i.e.
#' instead of dates, we have periods (corresponding to itemsets)
#'
#' @param df The preaggregated dataframe
#' @param summary_stats Logical controlling printing of summary
#' statistics regarding aggregation. Defaults to TRUE
#'
#' @return Returns a dataframe that has the properly classes dataframe
#' @export
#'
#' @examples pre_agg_demo %>% pre_aggregated()
pre_aggregated <- function(df, summary_stats = TRUE) {
    if (!(all(names(df) == c("id", "period", "event")))) {
        stop("There should be 3 columns named id, period and event (all lower case).")
    }

    stopifnot(is.integer(df$period))
    # if(!){ stop('period should be numeric') }

    if (!is.character(df$event))
        df$event = as.character(df$event)


    class(df) <-
        c("Aggregated_Dataframe", "tbl_df", "tbl", "data.frame")
    if (summary_stats) {
        message("Generating summary statistics of aggregated data...")
        generate_summary_stats(df)
    }

    df
}

#' Aggregation functions
#'
#' A dataframe having id, date and event column
#' can be aggregated using this function
#'
#' @param unaggregated_data A dataframe that has exactly 3 columns
#' in this order - id, date, event
#' @param format String specifying format of the date
#' @param calendar Boolean indicating whether or not to use calendar aggregation.
#' Defaults to false
#' @param unit String specifying unit of aggregation. Takes one of c('day', 'week',
#' 'month', '6 months', 'year')
#' @param n_units Integer specifying number of units.
#' @param anchor_table Beta
#' @param anchor_vector Beta
#' @param base_date Beta
#' @param occurence Beta
#' @param multiset Beta
#' @param include_date Beta
#' @param summary_stats Logical controlling printing of summary
#' statistics regarding aggregation. Defaults to TRUE
#'
#' @return Aggregated dataframe that has sequence id, itemset (period) and event (item)
#' @export
#'
#' @examples mvad %>%
#' aggregate_sequences(format = '%Y-%m-%d', unit = 'month', n_units = 6)
aggregate_sequences <-
    function(unaggregated_data,
             format = "%m-%d-%Y",
             calendar = FALSE,
             index_date = FALSE,
             unit = "week",
             n_units = 1,
             anchor_table = NA,
             anchor_vector = NA,
             base_date = NA,
             occurence = min,
             multiset = FALSE,
             include_date = FALSE,
             summary_stats = TRUE,
             output_directory = "~") {

        n_days <- get_n_days(unit, n_units)

        names(unaggregated_data) <- c("id", "date", "event")

        if (index_date) {

          unaggregated_data2 <- unaggregated_data

        } else {

          unaggregated_data2 <- unaggregated_data %>% mutate(date = readr::parse_date(date, format))
       }

        .GlobalEnv$env_dates <- new.env()
        .GlobalEnv$env_dates$df_unaggregated <- unaggregated_data2 %>% arrange(id,date)

        if (!is.na(anchor_table) || !is.na(anchor_vector)) {
            if (is.data.frame(anchor_table)) {
                names(anchor_table) <- c("id", "anchor_date")

                anchor_table <-
                    anchor_table %>% mutate(anchor_date =
                                                readr::parse_date(anchor_date, format))

            } else if (!is.na(anchor_vector)) {
                anchor_table <-
                    unaggregated_data2 %>% filter(event %in% anchor_vector) %>%
                    group_by(id) %>%
                    summarise(anchor_date = min(date))
            }

            aggregated_data <-
                unaggregated_data2 %>% left_join(anchor_table, by = "id") %>%
                group_by(id) %>%
                mutate(
                    n_ndays = (date - anchor_date) / n_days,
                    agg_n_ndays = if_else(n_ndays < 0, floor(n_ndays),
                                          ceiling(n_ndays)),
                    agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)
                ) %>% arrange(id,
                              agg_n_ndays) %>% mutate(agg_period = dense_rank(agg_n_ndays))

        } else if (calendar) {
            aggregated_data <-
                unaggregated_data2 %>%
                mutate(agg_date = lubridate::floor_date(date,
                                                        paste0(n_units,
                                                               " ", unit, "s"))) %>%
                group_by(id) %>%
                mutate(agg_period = dense_rank(agg_date))

        } else if (!is.na(base_date)) {
            if (typeof(base_date) == "closure") {
                base_date <- base_date(unaggregated_data2$date)
            }

            aggregated_data <-
                unaggregated_data2 %>% mutate(
                    n_ndays = (date - base_date) / n_days,
                    agg_n_ndays = if_else(n_ndays <
                                              0, floor(n_ndays), ceiling(n_ndays)),
                    agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)
                ) %>%
                group_by(id) %>% arrange(id, agg_n_ndays) %>%
                mutate(agg_period = dense_rank(agg_n_ndays))

        } else if (index_date) {

            aggregated_data <-
                unaggregated_data2 %>%
                group_by(id) %>%
                mutate(
                    n_ndays = date / (as.numeric(n_days) / (60 * 60 * 24)),
                    agg_n_ndays = if_else(n_ndays < 0, floor(n_ndays), ceiling(n_ndays)),
                    agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)
                ) %>%
                print() %>%
                arrange(id, agg_n_ndays) %>%
                mutate(agg_period = dense_rank(agg_n_ndays))

        }  else {
            aggregated_data <-
                unaggregated_data2 %>%
                group_by(id) %>%
                mutate(
                    n_ndays = (date - occurence(date)) / n_days,
                    agg_n_ndays = if_else(n_ndays < 0, floor(n_ndays), ceiling(n_ndays)),
                    agg_n_ndays = if_else(agg_n_ndays == 0, 1, agg_n_ndays)
                ) %>%
                print() %>%
                arrange(id, agg_n_ndays) %>%
                mutate(agg_period = dense_rank(agg_n_ndays))
        }



        if (!multiset) {
            aggregated_data <-
                aggregated_data %>% group_by(id, agg_period, event) %>% slice(1)
        }



        if (include_date) {
            aggregated_data <-
                aggregated_data %>% select(id, date, period = agg_period, event) %>%
                arrange(id,period, date, event)
        } else {
            aggregated_data <-
                aggregated_data %>% select(id, period = agg_period, event) %>%
                arrange(id, period, event)
        }

        class(aggregated_data) <-
            c("Aggregated_Dataframe", class(aggregated_data))
        if (summary_stats) {
            message("Generating summary statistics of aggregated data...")

            output_directory <-
                create_folder(output_directory, "approxmap_results")
            output_directory_public <-
                create_folder(output_directory, "public")
            sink(paste0(
                output_directory_public,
                "/",
                file_check(output_directory_public, "summary_text.txt")
            ),
            split = TRUE)
            generate_summary_stats(aggregated_data)
            sink()
        }


        aggregated_data
    }




#' Dataframe to sequence
#'
#' @param df_seq The aggregated dataframe that has the sequence ids, itemset ids and events.
#'
#' @return Returns a dataframe with a sequence object
#' (and a formatted sequence for readability)
#' corresponding to each id
#' @export
#'
#' @examples pre_agg_demo %>%
#' pre_aggregated %>%
#' convert_to_sequence()
convert_to_sequence <- function(df_seq) {
    if (!"Aggregated_Dataframe" %in% class(df_seq)) {
        warning("Are you sure the sequence dataframe you passed is already aggregated?")
    }

    # since the items have to be aggregated into itemsets and itemsets, into sequences,
    # this method has 2
    # mutate calls to do that all while ensuring the classes are appropriately maintained
    df_seq <-
        df_seq %>% group_by(id) %>% nest() %>% mutate(
            sequence = map(nested_id,
                           function(df_id) {
                               seqs <-
                                   df_id %>% group_by(period) %>%
                                   nest() %>%
                                   mutate(seqs = map(list_data,
                                                     function(itemset) {
                                                         events <- itemset$event
                                                         class(events) <-
                                                             c("Sequence_Itemset",
                                                               class(events))
                                                         events
                                                     })) %>% pull(seqs)
                               class(seqs) <-
                                   c("Sequence", class(seqs))
                               seqs
                           }),
            sequence_formatted = map_chr(sequence, format_sequence)
        ) %>% select(id, sequence, sequence_formatted)
    names(df_seq$sequence) <- df_seq$id
    class(df_seq$sequence) <-
        c("Sequence_List", class(df_seq$sequence))

    df_seq
}

#' Print methods
#'
#' A generic function that print methods use to display formatted results.
#'
#' @param x A sequence,w_sequence, w_sequence_list,
#' w_sequence_pattern or a w_sequence_dataframe object
#' @param ... Additional parameters to the function that is invoked
#'
#' @export
format_sequence <- function(x, ...) {
    UseMethod("format_sequence")
}

#' Print methods
#'
#' @param sequence A Sequence object
#'
#' @return Returns the sequence as a formatted string
#' @export
#'
#' @examples
#' sequences <- pre_agg_demo %>% pre_aggregated() %>%
#'     convert_to_sequence() %>% pull(sequence)
#' format_sequence(sequences[[1]])
format_sequence.Sequence <- function(sequence) {
    sequence <- sequence %>% map_chr(function(itemset) {
        itemset <- str_c(itemset, collapse = ", ")
        paste0("(", itemset, ")")
    }) %>% str_c(collapse = " ")
    paste0("<", sequence, ">")
}

#' Print methods
#'
#' @param sequence A Sequence object
#'
#' @return Null. Just prints the sequence as a formatted string
#' using the format_sequence method
#'
#' @export
#'
#' @examples
#' sequences <- pre_agg_demo %>% pre_aggregated() %>%
#'     convert_to_sequence() %>% pull(sequence)
#' print(sequences[[1]])
print.Sequence <- function(sequence) {
    sequence %>% format_sequence() %>% cat()
}

#' Print methods
#'
#' @param sequence A Sequence object
#'
#' @return Null. Just prints the sequence as a formatted string
#' using the format_sequence method
#'
#' @export
#'
#' @examples
#' sequences <- pre_agg_demo %>% pre_aggregated() %>%
#'     convert_to_sequence()
#' sequences <- sequences$sequence
#' print(sequences)
print.Sequence_List <- function(sequences) {
    if (is.null(names(sequences))) {
        warning("id for the sequences not present")
        walk(sequences, function(sequence_obj) {
            cat(format_sequence(sequence_obj), "\n")
        })
    } else {
        walk2(sequences, names(sequences), function(sequence_obj, id) {
            cat(id, ": ", format_sequence(sequence_obj), "\n")
        })
    }
}


#' Print methods
#'
#' Just a method to override default print and print the raw data structure.
#' Can be used to examine what it looks like internally
#'
#' @param obj Any object (useful for s3 classes)
#'
#' @export
#'
#' @examples
#' sequences <- pre_agg_demo %>% pre_aggregated() %>%
#'     convert_to_sequence() %>% pull(sequence)
#' print_raw(sequences[[1]])
print_raw <- function(obj) {
    obj %>% unclass() %>% print()
}
