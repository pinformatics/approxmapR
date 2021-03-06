#' @export
format_sequence.W_Sequence_Dataframe <-
    function(df,
             compare = FALSE,
             truncate_patterns = FALSE,
             html_format = FALSE) {
        column_patterns <- names(df)[str_detect(names(df), "_pattern")]
        if ("weighted_sequence" %in% names(df)) {
            columns <- c(column_patterns, "weighted_sequence")
        } else {
            columns <- column_patterns
        }

        if (truncate_patterns) {
            df <- df %>% mutate_at(column_patterns, truncate_pattern)
        }



        df <-
            df %>% select(one_of("cluster", "n", columns)) %>% mutate_at(columns, function(x)
                format_sequence(x,
                                html_format = html_format)) %>% mutate(n = as.double(n), n_percent = str_c(round(n /
                                                                                                                     sum(n) * 100,
                                                                                                                 digits = 2), "%")) %>% select(one_of("cluster", "n", "n_percent", columns))

        if (compare) {
            compare_sequences(df)
        } else {
            df
        }

    }

#' @export
view_formatted_sequence <- function(seq) {
    format_sequence(seq, html = TRUE) %>% stringr:::str_view_widget()
}

#' @export
compare_sequences <- function(df) {
    df %>% gather(-cluster,-n,-n_percent, key = "pattern", value = "sequence") %>% arrange(cluster) %>%
        mutate(pattern = stringr::str_replace(pattern, "_pattern", ""))
}




#' @export
class_it <- function(obj, class_name) {
    class(obj) <- c(class_name, class(obj)) %>% unique()
    obj
}

#' @export
truncate_pattern <- function(x, ...) {
    UseMethod("truncate_pattern")
}

#' @export
truncate_pattern.W_Sequence_Pattern_List <-
    function(w_sequence_list) {

        class_it(map(w_sequence_list, truncate_pattern),
                 "W_Sequence_List")
    }

## [ISSUE HERE]
#' @export
truncate_pattern.W_Sequence_Pattern <- function(w_sequence) {

    #browser()


    truncate_index <- rep(FALSE, length(w_sequence))
    for (i in seq_along(w_sequence)) {
        if (i == length(w_sequence))
            (break)()
        e_1 <- sort(w_sequence[[i]]$elements)
        e_2 <- sort(w_sequence[[i + 1]]$elements)
        if (identical(e_1, e_2)) {
            truncate_index[i] <- TRUE
        }
    }
    w_sequence[truncate_index] <- NULL

    ## May need to uncomment this out
    #compressed_n <-
    #    (truncate_index %>%
    #         as.integer() %>%
    #         as.character() %>%
    #         str_c(collapse = "") %>%
    #         str_split("0") %>%
    #         pluck(1) %>%
    #         str_subset(".") %>%
    #         str_count("1")) + 1

    #for(i in seq(1,length(w_sequence))){
    #    w_sequence[[i]]$itemset_weight <- compressed_n[i]
    #}


    w_sequence %>% class_it("W_Sequence_Pattern_Compressed")
}

format_sequence.W_Sequence_Pattern_Compressed <-
    function(w_sequence, html_format = FALSE) {
        n <- attr(w_sequence, "n")
        if (html_format) {
            if(n > 1){
                colors <-
                    rev(colormap::colormap(colormap = "viridis", nshades = n) %>%
                            stringr::str_sub(1, -3))
            } else {
                colors <- colormap::colormap(nshades = 2)[1]
            }


            # cuts <- floor(n*seq(0,1,0.2))[2:5]
            w_sequence %>%
                map_chr(function(w_itemset) {
                    tibble(
                        element = as.character(w_itemset$elements),
                        weight = as.integer(w_itemset$element_weights)
                    ) %>%
                        mutate(
                            ratio = weight / n,
                            # i = ceiling(ratio),
                            # color = map_chr(i, ~colors[.]),
                            color = colors[weight],
                            font_size = paste0(floor((1 + ratio * .6) * 100), "%"),
                            font_weight = signif(460 + ratio * 340, 1),
                            otag = str_c(
                                '<span style="',
                                "color: ",
                                color,
                                "; ",
                                "font-size: ",
                                font_size,
                                "; ",
                                "font-weight: ",
                                font_weight,
                                ";",
                                '">'
                            ),
                            ctag = "</span>",
                            element_html = str_c(otag, element, ":", weight, ctag)
                        ) %>%
                        pull(element_html) %>%
                        str_c(collapse = ", ") %>%
                        paste0("(", ., ")", ":", w_itemset$itemset_weight)
                }) %>%
                str_c(collapse = " ") %>%
                paste0("<", ., ">", " : ", n) %>%
                stringr::str_replace("<\\(", " < ( ")

        } else{
            w_sequence %>%
                map_chr(function(w_itemset) {
                    if (length(w_itemset$elements) > 0) {
                        str_c(w_itemset$elements, ":", w_itemset$element_weights) %>%
                            str_c(collapse = ", ") %>%
                            paste0("(", ., ")", ":", w_itemset$itemset_weight)
                    }
                    else{
                        NA
                    }

                }) %>%
                .[!is.na(.)] %>%
                str_c(collapse = " ") %>%
                paste0("<", ., ">", " : ", n)
        }

    }


#' @export
w_sequence_to_tibble <- function(w_sequence) {
    tibble(
        element = map(w_sequence, "elements") %>% unlist(),
        element_weight = map(w_sequence, "element_weights") %>%
            unlist(),
        itemset = map2(1:length(w_sequence), w_sequence, ~ rep(.x, length(.y$elements))) %>% unlist()
    ) %>%
        mutate(element_no = row_number())
}

#' @export
plot_weighted_sequence <- function(w_sequence) {
    df_sequence <- w_sequence %>% w_sequence_to_tibble()

    df_itemset <-
        df_sequence %>% group_by(itemset) %>% filter(element_no == max(element_no))

    df_sequence %>% ggplot(aes(element_no, element_weight)) + geom_point() + geom_label(aes(y = element_weight +
                                                                                                0.02 * element_weight, label = element)) + geom_vline(data = df_itemset, aes(xintercept = element_no))
}





#' @export
convert_to_events <- function(data, id_column, sequence_column) {
  
  data %>%
    mutate(event_set = str_split(data[[sequence_column]], "[\\(\\)]")) %>%
    unnest(cols = c(event_set)) %>%
    filter(event_set != "") %>% filter(event_set != " ") %>%
    group_by({{ sequence_column }}) %>%
    mutate(period = row_number()) %>%
    mutate(event = str_split(event_set, "[, ]")) %>%
    unnest(cols = c(event)) %>%
    filter(event != "") %>% filter(event != " ") %>% ungroup() %>%
    select(id_column, period, event)

}
