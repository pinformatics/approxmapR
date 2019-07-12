// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calculate_sorenson_distance_cpp
float calculate_sorenson_distance_cpp(StringVector itemset_1, StringVector itemset_2);
RcppExport SEXP _approxmapR_calculate_sorenson_distance_cpp(SEXP itemset_1SEXP, SEXP itemset_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type itemset_1(itemset_1SEXP);
    Rcpp::traits::input_parameter< StringVector >::type itemset_2(itemset_2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_sorenson_distance_cpp(itemset_1, itemset_2));
    return rcpp_result_gen;
END_RCPP
}
// dist_bw_sequences_cpp
float dist_bw_sequences_cpp(List sequence_1, List sequence_2);
RcppExport SEXP _approxmapR_dist_bw_sequences_cpp(SEXP sequence_1SEXP, SEXP sequence_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type sequence_1(sequence_1SEXP);
    Rcpp::traits::input_parameter< List >::type sequence_2(sequence_2SEXP);
    rcpp_result_gen = Rcpp::wrap(dist_bw_sequences_cpp(sequence_1, sequence_2));
    return rcpp_result_gen;
END_RCPP
}
// inter_sequence_distance_cpp
NumericMatrix inter_sequence_distance_cpp(List sequences);
RcppExport SEXP _approxmapR_inter_sequence_distance_cpp(SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(inter_sequence_distance_cpp(sequences));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_approxmapR_calculate_sorenson_distance_cpp", (DL_FUNC) &_approxmapR_calculate_sorenson_distance_cpp, 2},
    {"_approxmapR_dist_bw_sequences_cpp", (DL_FUNC) &_approxmapR_dist_bw_sequences_cpp, 2},
    {"_approxmapR_inter_sequence_distance_cpp", (DL_FUNC) &_approxmapR_inter_sequence_distance_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_approxmapR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}