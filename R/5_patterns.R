#' @export
filter_pattern <- function(x, ...) {
  UseMethod("filter_pattern")
}

#' @export
filter_pattern.Clustered_Dataframe <- function(df_clustered, ...) {
  df_clustered %>%
    get_weighted_sequence() %>%
    filter_pattern(...)
}

#' @export
filter_pattern.W_Sequence_Dataframe <-
  function(df_w_sequences, pattern_name = "consensus", ...) {
    if (!("weighted_sequence" %in% names(df_w_sequences))) {
      stop("Please use a dataframe where weighted sequences can be calculated")
    }


    pattern <- paste0(quo_name(pattern_name), "_pattern")

    df_pattern <-
      df_w_sequences %>%
      mutate(pat_name = structure(
                          map(weighted_sequence, filter_pattern, ... = ...)
                          # Commented out the two lines below on 11.2.2020
                          #,
                          #class = c("W_Sequence_List")
                          #
                        )
    )

    df_pattern$pat_name = class_it(df_pattern$pat_name, "W_Sequence_Pattern_List")

    df_pattern <-
      df_pattern %>%
      mutate(!!pattern := pat_name) %>%
      select(cluster, n, !!pattern, everything(), -pat_name)

    class(df_pattern) <- class(df_w_sequences)

    # print(attr(df_w_sequences$weighted_sequence[[1]], "alignments"))


    df_pattern
  }

#' @export
filter_pattern.W_Sequence <- function(w_sequence,
                                      threshold = 0.5,
                                      noise_threshold = 0,
                                      blank_if_absent = F,
                                      pattern_name = NULL,
                                      output_w_sequence = FALSE) {
  n <- attr(w_sequence, "n")

  elements <- unlist(map(w_sequence,
                         function(x) {
                           if (length(x$elements) == 0)
                             return(NULL)
                           x$elements
                         }))
  if (length(elements) == 0)
    return(weighted_seq)
  min_occurences <- n * threshold
  min_occurences <- max(min_occurences, noise_threshold)



  pattern <-
    map(w_sequence, function(w_sequence_itemset) {
      threshold_check <-
        w_sequence_itemset$element_weights >= min_occurences

      if (blank_if_absent && (sum(threshold_check) == 0)) {
        w_sequence_itemset$elements <- ""
        w_sequence_itemset$element_weights <- ""
      } else if (sum(threshold_check) == 0) {
        # w_sequence_itemset$elements <- NULL
        # w_sequence_itemset$element_weights <- NULL
        # w_sequence_itemset$itemset_weight <- NULL

        w_sequence_itemset <- NULL
      } else {
        w_sequence_itemset$elements <-
          w_sequence_itemset$elements[threshold_check]
        w_sequence_itemset$element_weights <-
          w_sequence_itemset$element_weights[threshold_check]
      }

      w_sequence_itemset
    })

  pattern <- pattern[map_lgl(pattern,  ~ !is.null(.))]


  attr(pattern, "n") <- n

  pattern <- class_it(pattern, "W_Sequence")
  if (!output_w_sequence)
    pattern <- class_it(pattern, "W_Sequence_Pattern")

  pattern
}

#' @export
format_sequence.W_Sequence_Pattern <-
  function(w_sequence_pattern, html_format = FALSE) {
    n <- attr(w_sequence_pattern, "n")

    if (html_format) {
      if(n > 1){
        colors <-
          rev(colormap::colormap(colormap = "bluered", nshades = n) %>%
                stringr::str_sub(1, -3))
      } else {
        colors <- colormap::colormap(nshades = 2)[1]
      }

      # cuts <- floor(n*seq(0,1,0.2))[2:5]
      w_sequence_pattern %>%
        map_chr(function(w_itemset) {
          tibble(
            element = as.character(w_itemset$elements),
            weight = as.integer(w_itemset$element_weights)
          ) %>%
            mutate(
              ratio = weight / n,
              # color = colors[floor(ratio*n)],
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
              element_html = str_c(otag, element, ctag)
            ) %>%
            pull(element_html) %>%
            str_c(collapse = ", ") %>%
            paste0("(", ., ")", "<br>")
        }) %>%
        str_c(collapse = " ")

    } else{
      w_sequence_pattern %>%
        map_chr(function(w_itemset) {
          if (length(w_itemset$elements) > 0) {
            str_c(w_itemset$elements, collapse = ", ") %>%
              paste0("(", ., ")")
          }
          else {
            NA
          }

        }) %>%
        .[!is.na(.)] %>%
        str_c(collapse = " ")
    }
  }



#' @export
format_sequence.W_Sequence_Pattern_List <-
  function(w_sequence_pattern_list, ...) {
    map_chr(w_sequence_pattern_list, function(w_sequence_pattern) {
      format_sequence(w_sequence_pattern, ...)
    })
  }

#' @export
print.W_Sequence_Pattern <- function(w_sequence_pattern, ...) {
  format_sequence(w_sequence_pattern) %>%
    cat()
}

#' @export
print.W_Sequence_Pattern_List <-
  function(w_sequence_pattern_list, ...) {
    walk(w_sequence_pattern_list, function(w_sequence_pattern) {
      cat(format_sequence(w_sequence_pattern), "\n")
    })
  }





#' @export
format_sequence.percentage_bar <-
  function(w_sequence_pattern, html_format = FALSE) {
    n <- attr(w_sequence_pattern, "n")

    if (html_format) {
      if(n > 1){
        colors <-
          rev(colormap::colormap(colormap = "bluered", nshades = n) %>%
                stringr::str_sub(1, -3))
      } else {
        colors <- colormap::colormap(nshades = 2)[1]
      }

      # cuts <- floor(n*seq(0,1,0.2))[2:5]
      w_sequence_pattern %>%
        map_chr(function(w_itemset) {
          tibble(
            element = as.character(w_itemset$elements),
            weight = as.integer(w_itemset$element_weights)
          ) %>%
            mutate(
              ratio = weight / n,
              # color = colors[floor(ratio*n)],
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
              element_html = str_c(otag, element, ctag)
            ) %>%
            pull(element_html) %>%
            str_c(collapse = ", ") %>%
            paste0("(", ., ")")
        }) %>%
        str_c(collapse = " ")

    } else {
      w_sequence_pattern %>%
        map_chr(function(w_itemset) {
          if (length(w_itemset$elements) > 0) {
            str_c(w_itemset$elements, collapse = ", ") %>%
              paste0("(", ., ")")
          }
          else {
            NA
          }

        }) %>%
        .[!is.na(.)] %>%
        str_c(collapse = " ")
    }
  }
