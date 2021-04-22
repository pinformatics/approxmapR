#' @export
get_weighted_sequence <- function(x, ...) {
  UseMethod("get_weighted_sequence")
}

#' @export
get_weighted_sequence.Sequence <- function(sequence_1, sequence_2) {
  w_sequence <-
    map2(sequence_1, sequence_2,
         function(sequence_itemset_1,
                  sequence_itemset_2) {
           w_sequence_itemset = list()
           elements <- c(sequence_itemset_1, sequence_itemset_2)
           elements <- elements[elements != "_"]
           freq_tbl <-
             elements %>%
             table()
           w_sequence_itemset$elements <- names(freq_tbl)
           w_sequence_itemset$element_weights <- unname(freq_tbl)
           if (("_" %in% sequence_itemset_1) |
               ("_" %in% sequence_itemset_2)) {
             w_sequence_itemset$itemset_weight <- 1
           } else {
             w_sequence_itemset$itemset_weight <- 2
           }
           class_it(w_sequence_itemset, "W_Sequence_Itemset")
         })
  attr(w_sequence, "n") <- 2
  class_it(w_sequence, "W_Sequence")
}

#' @export
get_weighted_sequence.W_Sequence <- function(w_sequence, sequence) {
  n <- attr(w_sequence, "n")
  x <- F
  w_sequence_new <-
    map2(w_sequence, sequence,
         function(w_sequence_itemset, sequence_itemset) {
           if ("_" %in% sequence_itemset) {
             x <- T
             class_it(w_sequence_itemset, "W_Sequence_Itemset")

           } else if ("_" %in% w_sequence_itemset$elements) {
             x <- T
             freq_tb <- table(sequence_itemset)
             w_sequence_itemset$elements <- names(freq_tb)
             w_sequence_itemset$element_weights <- unname(freq_tb)
             w_sequence_itemset$itemset_weight <- 1
             class_it(w_sequence_itemset, "W_Sequence_Itemset")
           } else {
             freq_tb <-
               rep(w_sequence_itemset$elements,
                   w_sequence_itemset$element_weights) %>%
               c(sequence_itemset) %>%
               table()

             w_sequence_itemset$elements <- names(freq_tb)
             w_sequence_itemset$element_weights <- unname(freq_tb)
             w_sequence_itemset$itemset_weight <-
               w_sequence_itemset$itemset_weight + 1
             class_it(w_sequence_itemset, "W_Sequence_Itemset")
           }
         })


  attr(w_sequence_new, "n") <- n + 1
  class_it(w_sequence_new, "W_Sequence")
}

#' @export
get_weighted_sequence.Sequence_List <- function(sequence_list,
                                                fun = sorenson_distance) {
  if (length(sequence_list) == 1) {
    # browser()
    sequence <- sequence_list[[1]]
    w_sequence <- align_sequences(sequence,
                                  sequence,
                                  fun)

    for (i in 1:length(sequence)) {
      w_sequence[[i]]$itemset_weight <- 1
      w_sequence[[i]]$element_weights <-
        (w_sequence[[i]]$element_weights) / 2
    }
    attr(w_sequence, "n") <- 1
    alignments <-
      attr(w_sequence, "alignments")[1]

  } else {
    # browser()
    w_sequence <-
      align_sequences(sequence_list[[1]], sequence_list[[2]], fun)

    if (length(sequence_list) > 2) {
      for (i in 3:length(sequence_list)) {
        w_sequence <-  align_sequences(w_sequence, sequence_list[[i]], fun)
      }
    }

    # browser()
    alignments <- attr(w_sequence, "alignments")
  }

  alignments <- class_it(alignments, "Sequence_List")
  names(alignments) <- names(sequence_list)
  attr(w_sequence, "alignments") <- alignments

  w_sequence
}

#' @export
get_weighted_sequence.Aggregated_Dataframe <-
  function(aggregated_dataframe,
           fun = sorenson_distance) {
    aggregated_dataframe %>%
      convert_to_sequence() %>%
      .$sequence %>%
      get_weighted_sequence(fun)
  }

#' @export
get_weighted_sequence.data.frame <- function(dataframe,
                                             fun = sorenson_distance) {
  warning(
    "Assuming that the dataframe you are using is pre-aggregated. If not, please use use the aggregate_sequences function.\n"
  )
  dataframe %>%
    pre_aggregated() %>%
    get_weighted_sequence(fun)
}

#' @export
get_weighted_sequence.Clustered_Dataframe <- function(df_clusters,
                                                      fun = sorenson_distance) {
  df_clusters <-
    df_clusters %>%
    mutate(weighted_sequence =
             map(df_sequences, function(df_sequence) {
               get_weighted_sequence(df_sequence$sequence, fun)
             }))

  class(df_clusters$weighted_sequence) <-
    c("W_Sequence_List",
      class(df_clusters$weighted_sequence))

  class(df_clusters) <-
    c("W_Sequence_Dataframe", class(df_clusters))

  df_clusters
}

#' @export
align_sequences <- function(x, ...) {
  UseMethod("align_sequences")
}

#' @export
align_sequences.Sequence <-
  function(sequence_1, sequence_2, fun = sorenson_distance) {
    distance_matrix <-
      inter_sequence_distance(sequence_1, sequence_2, fun)$distance_matrix
    aligned_sequence_1 <- structure(list(), class = "Sequence")
    aligned_sequence_2 <- structure(list(), class = "Sequence")

    i <- length(sequence_1) + 1
    j <- length(sequence_2) + 1



    backtrack <-
      function(i,
               j,
               aligned_sequence_1,
               aligned_sequence_2,
               operations = "") {
        if ((i == 1) & (j == 1)) {
          return(
            list(
              aligned_sequence_1 = structure(aligned_sequence_1, class = "Sequence"),
              aligned_sequence_2 = structure(aligned_sequence_2, class = "Sequence"),
              operations = operations
            )
          )

        }

        if ((i > 1) & (j > 1)) {
          #is left plus cost?
          if (distance_matrix[i, j] == distance_matrix[i, j - 1] + 1)  {
            #left
            aligned_sequence_1 <-
              append(aligned_sequence_1,
                     class_it("_", "Sequence_Itemset"),
                     0)
            aligned_sequence_2 <-
              append(aligned_sequence_2, sequence_2[j - 1][1], 0)
            backtrack(i,
                      j - 1,
                      aligned_sequence_1,
                      aligned_sequence_2,
                      operations)

          } else if (distance_matrix[i, j] == distance_matrix[i - 1, j] + 1) {
            #is up plus cost
            aligned_sequence_2 <-
              append(aligned_sequence_2,
                     class_it("_", "Sequence_Itemset"),
                     0)
            aligned_sequence_1 <-
              append(aligned_sequence_1, sequence_1[i - 1][1], 0)
            backtrack(i - 1,
                      j,
                      aligned_sequence_1,
                      aligned_sequence_2,
                      operations)

          } else {
            #diag
            aligned_sequence_1 <-
              append(aligned_sequence_1, sequence_1[i - 1][1], 0)
            aligned_sequence_2 <-
              append(aligned_sequence_2, sequence_2[j - 1][1], 0)
            backtrack(i - 1,
                      j - 1,
                      aligned_sequence_1,
                      aligned_sequence_2,
                      operations)

          }
        } else if ((i == 1) & (j > 1)) {
          aligned_sequence_1 <-
            append(aligned_sequence_1,
                   class_it("_", "Sequence_Itemset"),
                   0)
          aligned_sequence_2 <-
            append(aligned_sequence_2, sequence_2[j - 1][1], 0)
          backtrack(i,
                    j - 1,
                    aligned_sequence_1,
                    aligned_sequence_2,
                    operations)

        } else if ((j == 1) & (i > 1)) {
          aligned_sequence_2 <-
            append(aligned_sequence_2,
                   class_it("_", "Sequence_Itemset"),
                   0)
          aligned_sequence_1 <-
            append(aligned_sequence_1, sequence_1[i - 1][1], 0)
          backtrack(i - 1,
                    j,
                    aligned_sequence_1,
                    aligned_sequence_2,
                    operations)

        }
      }
    # debug(backtrack)
    aligned_sequences <-
      backtrack(i, j, aligned_sequence_1, aligned_sequence_2)

    w_sequence <-
      get_weighted_sequence(aligned_sequences$aligned_sequence_1,
                            aligned_sequences$aligned_sequence_2)

    attr(w_sequence, "alignments") <-
      list(aligned_sequences$aligned_sequence_1,
           aligned_sequences$aligned_sequence_2)

    w_sequence

  }


#' @export
insert_blank_w_itemset <- function(w_sequence) {
  blank_w_itemset <- class_it(list(list(
    elements = "_",
    element_weights = 0,
    itemset_weight = 0
  )), "W_Sequence_Itemset")

  append(w_sequence, blank_w_itemset, 0)
}

#' @export
align_sequences.W_Sequence <- function(w_sequence,
                                       sequence,
                                       fun = sorenson_distance) {
  distance_matrix <-
    inter_sequence_distance(w_sequence, sequence, fun)$distance_matrix
  n <- attr(w_sequence, "n")
  pre_alignments <- attr(w_sequence, "alignments")

  aligned_sequence <- list()
  aligned_w_sequence <- list()
  aligned_alignments <- list()

  i <- length(sequence) + 1
  j <- length(w_sequence) + 1

  emp <- class_it("_", "Sequence_Itemset")

  backtrack <- function(i,
                        j,
                        aligned_sequence,
                        aligned_w_sequence,
                        aligned_alignments = list()) {
    if ((i == 1) & (j == 1)) {
      result = list(
        aligned_sequence = structure(aligned_sequence,
                                     class = "Sequence"),
        aligned_w_sequence = structure(
          aligned_w_sequence,
          class = "W_Sequence",
          n = attr(w_sequence, "n")
        ),
        aligned_alignments = aligned_alignments
      )
      return(result)
    }

    if ((i > 1) & (j > 1)) {
      h =  sorenson_distance(w_sequence[[j - 1]], sequence[i - 1], n)
      if (distance_matrix[i, j] == (distance_matrix[i - 1, j - 1] + h)) {
        aligned_sequence <- append(aligned_sequence, sequence[i - 1][1], 0)
        aligned_w_sequence <-
          append(aligned_w_sequence, w_sequence[j - 1][1], 0)
        if (length(aligned_alignments) == 0) {
          aligned_alignments <- map(pre_alignments, j - 1) %>% map( ~ list(.))
        } else{
          alignment_elements <- map(pre_alignments, j - 1)
          aligned_alignments <-
            map2(aligned_alignments,
                 alignment_elements,
                 ~ append(.x, list(.y), 0))
        }

        backtrack(i - 1,
                  j - 1,
                  aligned_sequence,
                  aligned_w_sequence,
                  aligned_alignments)

      } else if (distance_matrix[i, j] == distance_matrix[i, j - 1] + 1)  {
        #is left plus cost?
        aligned_sequence <- append(aligned_sequence, emp, 0)
        aligned_w_sequence <-
          append(aligned_w_sequence, w_sequence[j - 1][1], 0)
        if (length(aligned_alignments) == 0) {
          aligned_alignments <- map(pre_alignments, j - 1) %>% map( ~ list(.))
        } else{
          alignment_elements <- map(pre_alignments, j - 1)
          aligned_alignments <-
            map2(aligned_alignments,
                 alignment_elements,
                 ~ append(.x, list(.y), 0))
        }

        backtrack(i,
                  j - 1,
                  aligned_sequence,
                  aligned_w_sequence,
                  aligned_alignments)

      } else if (distance_matrix[i, j] == distance_matrix[i - 1, j] + 1) {
        # browser()

        #is up plus cost
        aligned_w_sequence <-
          insert_blank_w_itemset(aligned_w_sequence)
        aligned_sequence <-
          append(aligned_sequence, sequence[i - 1][1], 0)
        if (length(aligned_alignments) == 0) {
          aligned_alignments <- rep(list(emp), length(pre_alignments))
        } else{
          aligned_alignments <- map(aligned_alignments,
                                    function(l) {
                                      l <- append(l, emp, 0)
                                      class(l[[1]]) <- class(emp)
                                      l
                                    })
        }

        backtrack(i - 1,
                  j,
                  aligned_sequence,
                  aligned_w_sequence,
                  aligned_alignments)

      } else {
        #diag
        aligned_sequence <-
          append(aligned_sequence, sequence[i - 1][1], 0)
        aligned_w_sequence <-
          append(aligned_w_sequence, w_sequence[j - 1][1], 0)
        if (length(aligned_alignments) == 0) {
          aligned_alignments <- map(pre_alignments, j - 1) %>% map( ~ list(.))
        } else{
          alignment_elements <- map(pre_alignments, j - 1)
          aligned_alignments <-
            map2(aligned_alignments,
                 alignment_elements,
                 ~ append(.x, list(.y), 0))
        }

        backtrack(i - 1,
                  j - 1,
                  aligned_sequence,
                  aligned_w_sequence,
                  aligned_alignments)

      }
    } else if ((i == 1) & (j > 1)) {
      aligned_sequence <- append(aligned_sequence, emp, 0)
      aligned_w_sequence <-
        append(aligned_w_sequence, w_sequence[j - 1][1], 0)
      if (length(aligned_alignments) == 0) {
        aligned_alignments <- map(pre_alignments, j - 1) %>% map( ~ list(.))
      } else{
        alignment_elements <- map(pre_alignments, j - 1)
        aligned_alignments <-
          map2(aligned_alignments,
               alignment_elements,
               ~ append(.x, list(.y), 0))
      }

      backtrack(i,
                j - 1,
                aligned_sequence,
                aligned_w_sequence,
                aligned_alignments)

    } else if ((j == 1) & (i > 1)) {
      # browser()
      aligned_w_sequence <-
        insert_blank_w_itemset(aligned_w_sequence)
      aligned_sequence <-
        append(aligned_sequence, sequence[i - 1][1], 0)
      if (length(aligned_alignments) == 0) {
        aligned_alignments <- rep(list(emp), length(pre_alignments))
      } else{
        aligned_alignments <- map(aligned_alignments,
                                  function(l) {
                                    l <- append(l, emp, 0)
                                    class(l[[1]]) <- class(emp)
                                    l
                                  })
      }

      backtrack(i - 1,
                j,
                aligned_sequence,
                aligned_w_sequence,
                aligned_alignments)
    }
  }
  # debug(backtrack)
  aligned_sequences <-
    backtrack(i, j, aligned_sequence, aligned_w_sequence)

  # message(get_weighted_sequence(aligned_sequences$aligned_w_sequence,
  # aligned_sequences$aligned_sequence))
  w_sequence <-
    get_weighted_sequence(aligned_sequences$aligned_w_sequence,
                          aligned_sequences$aligned_sequence)
  aligned_sequences <-
    append(
      aligned_sequences$aligned_alignments,
      list(aligned_sequences$aligned_sequence)
    )

  aligned_sequences <-
    aligned_sequences %>%
    map(function(x) {
      # x <- map(function(y){
      #
      # })
      class_it(x, "Sequence")
    })


  attr(w_sequence, "alignments") <- aligned_sequences

  w_sequence
}

#' @export
format_sequence.W_Sequence <-
  function(w_sequence, html_format = FALSE) {
    n <- attr(w_sequence, "n")
    if (html_format) {
      if(n > 1){
        colors <-
          rev(colormap::colormap(colormap = "bluered", nshades = n) %>%
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
            paste0("(", ., ")", ":", w_itemset$itemset_weight, "<br>")
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
format_sequence.W_Sequence_List <- function(w_sequence_list, ...) {
  map_chr(w_sequence_list, function(w_sequence) {
    format_sequence(w_sequence, ...)
  })
}

#' @export
print.W_Sequence <- function(w_sequence, ...) {
  format_sequence(w_sequence) %>%
    cat()
}

#' @export
print.W_Sequence_List <- function(w_sequence_list, ...) {
  walk(w_sequence_list, function(w_sequence) {
    cat(format_sequence(w_sequence), "\n")
  })
}
