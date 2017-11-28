#' Distance Calculation
#'
#' A genric that helps calculate sorenson distance between two itemsets
#' (either both from a Sequence or 1 from a sequence and another from a
#' weighted sequence)
#'
#' @param x A sequence or weighted sequence itemset
#' @param ... Contains the other itemset and other parameters to be passed
#' to the invoked parameter
#'
#' @return Returns the distance between the 2 itemsets
#' @export
sorenson_distance <- function(x, ...){
  UseMethod("sorenson_distance")
}


#' Distance Calculation
#'
#' @param sequence_itemset_1 The first sequence itemset
#' @param sequence_itemset_2 The second sequence itemset
#'
#' @return Returns the distance between the 2 itemsets
#' @export
sorenson_distance.Sequence_Itemset <- function(sequence_itemset_1, sequence_itemset_2) {

  set_1 <- setdiff(sequence_itemset_1,sequence_itemset_2)
  set_2 <- setdiff(sequence_itemset_2,sequence_itemset_1)
  set_union <- union(set_1,set_2)

  length(set_union) / (length(sequence_itemset_1) + length(sequence_itemset_2))
}

#' Distance Calculation
#'
#' @param w_sequence_itemset The weighted sequence itemset
#' @param sequence_itemset The sequence itemset
#'
#' @return Returns the distance between the 2 itemsets
#' @export
sorenson_distance.W_Sequence_Itemset <- function(w_sequence_itemset,sequence_itemset,n) {
  v <- w_sequence_itemset$itemset_weight
  w <- w_sequence_itemset$element_weights
  x <- w_sequence_itemset$elements
  y <- sequence_itemset
  y_no <- length(sequence_itemset)
  eR <- (sum(w) + (y_no * v) - (2 * sum(w[x %in% y]))) / (sum(w) + (y_no * v))

  ((eR * v) + n - v) / n
}


#' Distance Calculation
#'
#' Used to calculate replacement costs.
#' Passes the supplied function to the subsequent functions.
#'
#' @param x Either a sequence or weighted_sequence itemset
#' @param ... Has the other sequence itemset and function to be passed to calcaute costs
#'
#' @return Returns the replacement costs between the 2 itemsets
#' @export
repl <- function(x, ...){
  UseMethod("repl")
}


#' Distance Calculation
#'
#' @param sequence_itemset_1 The first sequence itemset
#' @param sequence_itemset_2 The second sequence itemset
#' @param fun Uses the sorenson distance by default.
#' Can pass other custom distance functions
#'
#' @return Returns the replacement costs between the 2 itemsets
#' @export
repl.Sequence_Itemset <- function(sequence_itemset_1, sequence_itemset_2, fun = sorenson_distance) {
  fun(sequence_itemset_1, sequence_itemset_2)
}

#' Title
#'
#' @param w_sequence_itemset The weighted sequence itemset
#' @param sequence_itemset The sequence itemset
#' @param n The number of sequences the weighted sequence from which the
#' weighted sequence itemset, is a summary of.
#' @param fun The distance measure to be to be used
#'

#' @export
repl.W_Sequence_Itemset <- function(w_sequence_itemset,
                                       sequence_itemset,
                                       n,
                                       fun = sorenson_distance) {
  fun(w_sequence_itemset,sequence_itemset, n)
}

#' @export
indel <- function(x, ...){
  UseMethod("indel")
}

#' @export
indel.Sequence_Itemset <- function(sequence_itemset,fun) {
  repl(sequence_itemset,"",fun)
}

#' @export
indel.W_Sequence_Itemset <- function(w_sequence_itemset, n, fun) {
  repl(w_sequence_itemset, "", n, fun)
}

#' @export
inter_sequence_distance <- function(x, ...){
  UseMethod("inter_sequence_distance")
}

#' @export
inter_sequence_distance.Sequence <- function(sequence_1,
                                             sequence_2,
                                             fun = sorenson_distance) {

  distance_matrix <- matrix(nrow = length(sequence_1) + 1,
                           ncol = length(sequence_2) + 1)

  distance_matrix[1,] <- 0:length(sequence_2)
  distance_matrix[,1] <- 0:length(sequence_1)

  for(i in 2:nrow(distance_matrix)) {
    for(j in 2:ncol(distance_matrix)) {
      replace <- distance_matrix[i-1,j-1] + repl(sequence_1[[i-1]],
                                                 sequence_2[[j-1]], fun)
      indel_r <- distance_matrix[i,j-1] + indel(sequence_2[[j-1]], fun)
      indel_d <- distance_matrix[i-1,j] + indel(sequence_1[[i-1]], fun)
      distance_matrix[i,j] <- min(replace, indel_d, indel_r)
    }
  }

  list(distance_matrix = distance_matrix,
       distance = distance_matrix[nrow(distance_matrix),ncol(distance_matrix)])
}

#' @export
inter_sequence_distance.W_Sequence <- function(w_sequence,
                                               sequence,
                                               fun = sorenson_distance) {

  n <- attr(w_sequence, "n")
  distance_matrix <- matrix(nrow = length(sequence) + 1,
                            ncol = length(w_sequence) + 1)

  distance_matrix[1,] <- 0:length(w_sequence)
  distance_matrix[,1] <- 0:length(sequence)

  for(i in 2:nrow(distance_matrix)) {
    for(j in 2:ncol(distance_matrix)) {
      sequence_itemset <- sequence[[i-1]]
      w_sequence_itemset <- w_sequence[[j-1]]
      replace <- distance_matrix[i-1,j-1] + repl(w_sequence_itemset, sequence_itemset, n, fun)
      indel_r <- distance_matrix[i,j-1] + 1 #indel(w_sequence_itemset, n, fun)
      indel_d <- distance_matrix[i-1,j] + 1 #indel(sequence_itemset, fun)
      distance_matrix[i,j] <- min(replace,indel_d,indel_r)
    }
  }

  list(distance_matrix = distance_matrix,
       distance = distance_matrix[nrow(distance_matrix), ncol(distance_matrix)])
}

#' @export
inter_sequence_distance.Sequence_List <- function(sequence_list){
  inter_sequence_distance_cpp(sequence_list)
}

