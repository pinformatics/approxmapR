#include <Rcpp.h>
// #include <algorithm>
// #include <string>
// #include <set>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace std;
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
float calculate_sorenson_distance_cpp(StringVector itemset_1, StringVector itemset_2)
{
  int seq_1_length = itemset_1.length();
  int seq_2_length = itemset_2.length();
  int seq_1_diff = 0, seq_2_diff = 0;

  for(int i = 0; i < seq_1_length; i++)
  {
    bool found  = FALSE;
    for(int j = 0; j < seq_2_length; j++)
    {
      if(itemset_1[i]==itemset_2[j])
      {
        found  = TRUE;
        break;
      }
    }
    if(!found)
    {
      seq_1_diff = seq_1_diff + 1;
    }
  }

  for(int i = 0; i < seq_2_length; i++)
  {
    bool found  = FALSE;
    for(int j = 0; j < seq_1_length; j++)
    {

      if(itemset_2[i]==itemset_1[j])
      {
        found  = TRUE;
        break;
      }
    }
    if(!found)
    {
      seq_2_diff = seq_2_diff + 1;
    }
  }

  float dist = ((float)(seq_1_diff + seq_2_diff)) / ((float)(seq_1_length + seq_2_length));
  return(dist);
}


//' @export
// [[Rcpp::export]]
float dist_bw_sequences_cpp(List sequence_1, List sequence_2)
{
  int seq_1_length = sequence_1.length();
  int seq_2_length = sequence_2.length();


  NumericMatrix distance_matrix(seq_1_length+1,seq_2_length+1);

  for(int i = 0; i<distance_matrix.nrow() ; i++)
  {
    distance_matrix(i,0) = i;
  }

  for(int i = 0; i<distance_matrix.ncol() ; i++)
  {
    distance_matrix(0,i) = i;
  }

  int max_length = max(seq_1_length,seq_2_length);
  //distance_matrix[(distance_matrix/max_length)>max_dist] = Inf


  for(int i = 1; i<(distance_matrix.nrow()); i++)
  {
    for(int j = 1; j<(distance_matrix.ncol()); j++)
    {

      float repl = ((float)(distance_matrix(i-1,j-1))) + calculate_sorenson_distance_cpp(sequence_1[i-1],sequence_2[j-1]);
      float indel_d = ((float)(distance_matrix(i-1,j))) + (float)1;
      float indel_r = ((float)(distance_matrix(i,j-1))) + (float)1;
      distance_matrix(i,j) = min(NumericVector::create(repl,indel_d,indel_r));
    }
  }

  float distance = distance_matrix(distance_matrix.nrow()-1,distance_matrix.ncol()-1) / max_length;
  return(distance);
}



//' @export
// [[Rcpp::export]]
NumericMatrix inter_sequence_distance_cpp(List sequences)
{
  int seq_list_length = sequences.length();
  Progress p(seq_list_length*seq_list_length, true);
  NumericMatrix distance_matrix(seq_list_length,seq_list_length);
  for(int i = 0; i < seq_list_length; i++)
  {
    for(int j = 0; j < seq_list_length; j++)
    {
      if (Progress::check_abort() )
        return -1.0;
      if(i==j){
        distance_matrix(i,j) = NumericVector::get_na();
      } else if(distance_matrix(j,i) != 0) {
        distance_matrix(i,j) = distance_matrix(j,i);
      } else {
        distance_matrix(i,j) = dist_bw_sequences_cpp(sequences[i],sequences[j]);
      }
      p.increment();
    }
  }
  return(distance_matrix);
}

