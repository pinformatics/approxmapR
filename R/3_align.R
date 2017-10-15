get_weighted_sequence <- function(x, ...){
  UseMethod("get_weighted_sequence")
}

get_weighted_sequence.Sequence <- function(sequence_1, sequence_2) {

  weighted_sequence <- list(length(sequence_1))

  for(i in 1:length(sequence_1)) {
    all_elements = c(sequence_1[[i]],sequence_2[[i]])
    elements = union(sequence_1[[i]],sequence_2[[i]])
    elements = elements[elements != "_"]

    element_weights = numeric(length(elements))
    for(j in 1:length(elements)) {
      element_weights[j] = sum(all_elements==elements[j])
    }

    if(("_" %in% sequence_1[[i]])|("_" %in% sequence_2[[i]])) {
      itemset_weight = 1
    } else {
      itemset_weight = 2
    }
    weighted_sequence[[i]] = list(elements = elements,
                                  element_weights = element_weights,
                                  itemset_weight = itemset_weight)
    class(weighted_sequence[[i]]) <- c("W_Sequence_Itemset", class(weighted_sequence[[i]]))
  }

  attr(weighted_sequence, "n") <- 2
  class(weighted_sequence) <- "W_Sequence"

  weighted_sequence
}

get_weighted_sequence.W_Sequence = function(w_sequence, sequence) {

  for(i in 1:length(w_sequence)) {
    w_sequence_itemset = w_sequence[[i]]
    itemset = sequence[[i]]
    if("_" %in% itemset) {
      next

    } else if("_" %in% w_sequence_itemset$elements) {

      w_sequence[[i]]$elements = unique(itemset)
      for(j in 1:length(w_sequence[[i]]$elements)) {
        w_sequence[[i]]$element_weights[j] = sum(w_sequence[[i]]$elements[j]==itemset)
      }
      w_sequence[[i]]$itemset_weight = 1

    } else {

      all_elements = NULL
      for(j in 1:length(w_sequence_itemset$elements)) {
        all_elements = c(all_elements,rep(w_sequence_itemset$elements[j],w_sequence_itemset$element_weights[j]))
      }

      all_elements = c(itemset,all_elements)
      unique_elements = unique(all_elements)

      w_sequence[[i]]$elements = unique_elements
      w_sequence[[i]]$element_weights = NULL

      for(j in 1:length(w_sequence[[i]]$elements)) {
        w_sequence[[i]]$element_weights[j] = sum(w_sequence[[i]]$elements[j]==all_elements)
      }

      w_sequence[[i]]$itemset_weight = w_sequence_itemset$itemset_weight + 1

    }
  }

  attr(w_sequence,"n") = attr(w_sequence,"n") + 1

  w_sequence
}


get_weighted_sequence.Sequence_List <- function(sequence_list,
                                                fun = sorenson_distance){

  if(length(sequence_list)==1) {
    sequence = sequence_list[[1]]
    w_sequence = vector(mode="list", length(sequence))
    for(i in 1:length(sequence)) {
      w_sequence[[i]]$elements = sequence[[i]]
      w_sequence[[i]]$itemset_weight = 1
      w_sequence[[i]]$element_weights = rep(1,length(sequence[[i]]))
    }
    attr(w_sequence, "n") <- 1

  } else {
    w_sequence = align_sequences(sequence_list[[1]],sequence_list[[2]], fun)

    if(length(sequence_list)>2) {
      for(i in 3:length(sequence_list)) {

        w_sequence =  align_sequences(w_sequence, sequence_list[[i]], fun)
      }
    }
  }

  w_sequence
}



align_sequences <- function(x, ...){
  UseMethod("align_sequences")
}

align_sequences.Sequence <- function(sequence_1, sequence_2, fun = sorenson_distance) {

  distance_matrix = inter_sequence_distance(sequence_1, sequence_2, fun)$distance_matrix
  aligned_sequence_1 = structure(list(), class = "Sequence")
  aligned_sequence_2 = structure(list(), class = "Sequence")

  i = length(sequence_1)+1
  j = length(sequence_2)+1



  backtrack <- function(i, j, aligned_sequence_1, aligned_sequence_2, operations = "") {
    if((i==1) & (j==1)) {
      return(list(aligned_sequence_1 = structure(aligned_sequence_1, class = "Sequence"),
                  aligned_sequence_2 = structure(aligned_sequence_2, class = "Sequence"),
                  operations = operations))

    }

    if((i>1) & (j>1)) {
      #is left plus cost?
      if (distance_matrix[i,j]==distance_matrix[i,j-1] + 1)  {
        #left
        aligned_sequence_1 = append(aligned_sequence_1, "_", 0)
        aligned_sequence_2 = append(aligned_sequence_2, sequence_2[j-1][1], 0)
        backtrack(i, j - 1, aligned_sequence_1, aligned_sequence_2, operations)

      } else if (distance_matrix[i,j]==distance_matrix[i-1,j]+1) {
        #is up plus cost
        aligned_sequence_2 = append(aligned_sequence_2, "_", 0)
        aligned_sequence_1 = append(aligned_sequence_1, sequence_1[i-1][1], 0)
        backtrack(i - 1, j, aligned_sequence_1, aligned_sequence_2, operations)

      } else {
        #diag
        aligned_sequence_1 = append(aligned_sequence_1, sequence_1[i-1][1],0)
        aligned_sequence_2 = append(aligned_sequence_2, sequence_2[j-1][1],0)
        backtrack(i - 1, j - 1, aligned_sequence_1, aligned_sequence_2, operations)

      }
    } else if ((i==1) & (j>1)) {

      aligned_sequence_1 = append(aligned_sequence_1, "_", 0)
      aligned_sequence_2 = append(aligned_sequence_2, sequence_2[j-1][1], 0)
      backtrack(i, j - 1, aligned_sequence_1, aligned_sequence_2, operations)

    } else if((j==1) & (i>1)) {

      aligned_sequence_2 = append(aligned_sequence_2, "_", 0)
      aligned_sequence_1 = append(aligned_sequence_1, sequence_1[i-1][1], 0)
      backtrack(i - 1, j, aligned_sequence_1, aligned_sequence_2, operations)

    }
  }
  # debug(backtrack)
  aligned_sequences <- backtrack(i, j, aligned_sequence_1, aligned_sequence_2)

  get_weighted_sequence(aligned_sequences$aligned_sequence_1,
                        aligned_sequences$aligned_sequence_2)

}

insert_blank_w_itemset <- function(w_sequence) {
  w_sequence_itemset = list(list(elements = "_",element_weights = 0, itemset_weight = 0))
  append(w_sequence, w_sequence_itemset,0)
}


align_sequences.W_Sequence = function(w_sequence,
                                      sequence,
                                      fun = sorenson_distance) {

  distance_matrix = inter_sequence_distance(w_sequence, sequence, fun)$distance_matrix
  n = attr(w_sequence, "n")

  aligned_sequence = list()
  aligned_w_sequence = list()

  i = length(sequence)+1
  j = length(w_sequence)+1

  backtrack = function(i, j, aligned_sequence, aligned_w_sequence, operations = "") {
    if((i==1) & (j==1))
    {
      result = list(aligned_sequence = structure(aligned_sequence,
                                                 class = "Sequence"),
                    aligned_w_sequence = structure(aligned_w_sequence,
                                                   class = "W_Sequence",
                                                   n = attr(w_sequence, "n")),
                    operations = operations)
      return(result)
    }

    if((i>1) & (j>1)) {
      h =  sorenson_distance(w_sequence[[j-1]], sequence[i-1], n)
      if(distance_matrix[i,j] == (distance_matrix[i-1,j-1] + h)) {

        aligned_sequence = append(aligned_sequence, sequence[i-1][1], 0)
        aligned_w_sequence = append(w_sequence[j-1][1], aligned_w_sequence, 0)
        backtrack(i - 1, j - 1, aligned_sequence, aligned_w_sequence, operations)

      } else if (distance_matrix[i,j] == distance_matrix[i, j-1] + 1)  {

        #is left plus cost?
        aligned_sequence = append(aligned_sequence, "_", 0)
        aligned_w_sequence = append(aligned_w_sequence, w_sequence[j-1][1], 0)
        backtrack(i, j - 1, aligned_sequence, aligned_w_sequence, operations)

      } else if (distance_matrix[i,j]==distance_matrix[i-1,j]+1) {

        #is up plus cost
        aligned_w_sequence = insert_blank_w_itemset(aligned_w_sequence)
        aligned_sequence = append(aligned_sequence, sequence[i-1][1], 0)
        backtrack(i - 1, j, aligned_sequence, aligned_w_sequence, operations)

      } else {

        #diag
        aligned_sequence = append(aligned_sequence, sequence[i-1][1], 0)
        aligned_w_sequence = append(aligned_w_sequence, w_sequence[j-1][1], 0)
        backtrack(i - 1, j - 1, aligned_sequence, aligned_w_sequence, operations)

      }
    } else if ((i==1) & (j>1)) {

      aligned_sequence = append(aligned_sequence, "_", 0)
      aligned_w_sequence = append(aligned_w_sequence, w_sequence[j-1][1], 0)
      backtrack(i, j - 1, aligned_sequence, aligned_w_sequence, operations)

    } else if((j==1) & (i>1)) {

      aligned_w_sequence = insert_blank_w_itemset(aligned_w_sequence)
      aligned_sequence = append(aligned_sequence, sequence[i-1][1], 0)
      backtrack(i - 1, j, aligned_sequence, aligned_w_sequence, operations)

    }
  }
  #debug(backtrack)
  aligned_sequences = backtrack(i, j, aligned_sequence, aligned_w_sequence)

  # message(get_weighted_sequence(aligned_sequences$aligned_w_sequence,
                        # aligned_sequences$aligned_sequence))
  get_weighted_sequence(aligned_sequences$aligned_w_sequence,
                        aligned_sequences$aligned_sequence)
}


format_w_sequence <- function(w_sequence){
  w_sequence %>%
    map_chr(function(w_itemset){
      str_c(w_itemset$elements, ":", w_itemset$element_weights) %>%
        str_c(collapse = ", ") %>%
        paste0("(", ., ")", ":", w_itemset$itemset_weight)
    }) %>%
    str_c(collapse = " ") %>%
    paste0("<", ., ">", " : ", attr(w_sequence, "n"))
}

print.W_Sequence <- function(w_sequence, ...){
  format_w_sequence(w_sequence) %>%
  cat()
}



