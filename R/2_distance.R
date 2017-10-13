sorenson_distance <- function(itemset1, itemset2) {

  set1 <- setdiff(itemset1,itemset2)
  set2 <- setdiff(itemset2,itemset1)
  set_union <- union(set1,set2)

  length(set_union) / (length(itemset1) + length(itemset2))

}


replace <- function(itemset1, itemset2, fun) {
  fun(itemset1,itemset2)
}



indel <- function(itemset,fun) {
  repl_btw_itemsets(itemset,"",fun)
}





inter_sequence_distance <- function(seq_1, seq_2,fun = sorenson_distance) {

  distance_matrix <- matrix(nrow = length(seq_1) + 1,
                           ncol = length(seq_2) + 1)

  distance_matrix[1,] <- 0:length(seq_2)
  distance_matrix[,1] <- 0:length(seq_1)

  for(i in 2:nrow(distance_matrix)) {
    for(j in 2:ncol(distance_matrix)) {
      repl <- distance_matrix[i-1,j-1] + repl_btw_itemsets(seq_1[[i-1]],seq_2[[j-1]],fun)
      indel_r <- distance_matrix[i,j-1] + indel(seq_2[[j-1]],fun)
      indel_d <- distance_matrix[i-1,j] + indel(seq_1[[i-1]],fun)
      distance_matrix[i,j] <- min(repl,indel_d,indel_r)
    }
  }

  list(distance_matrix = distance_matrix,
       distance = distance_matrix[nrow(distance_matrix),ncol(distance_matrix)])
}
