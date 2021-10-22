\name{generate_reports}
\alias{generate_reports}
\title{Report generation}
\usage{
  generate_reports(w_sequence_dataframe,  sil_table = NULL, html_format = TRUE,
                  output_directory = "~", end_filename_with = "",
                  sequence_analysis_details = NULL,
                  sequence_analysis_details_definitions == NULL,
                  algorithm_comparison = FALSE)
}
\arguments{
\item{w_sequence_dataframe}{A dataframe with class "W_Sequence_Dataframe". This will be
the dataframe that resulted from extracting the patterns after clustering.}

\item{sil_table}{The silhouette_object which is produced by the -find_optimal_k- function that contains the silhouette Information
for the K value that was selected.}

\item{html_format}{A boolean value to indicate if the exports should have HTML formatting.}

\item{output_directory}{The path to where the exports should be placed. This creates a folder with the name of "approxmap_results".}

\item{end_filename_with}{The option of appending to the end of the default file names. This
is useful if running multiple algorithms that will be exported to the same output_directory.}

\item{sequence_analysis_details}{This will generate a report that includes details pertaining to the sequence analysis. This
must be a list with the following structure list("algorithm" = "a string", "k_value" = a number, "time_period" = "a string", "consensus_threshold" = a number, "notes" = "Any special notes as a string.")}

\item{sequence_analysis_details_definitions}{This needs to be an data frame object with column 1 being labelled "event" which contains the events in the data, while column 2 can be any label which
contains the definitions or descriptions of the event.}

\item{algorithm_comparison}{The option to indicate if the report being generated is
one that is comparing multiple algorithms, for example the outcome of using the K-NN
and K-Medoids algorithm. This function separates the *id* column using
id \%>\% str_split("_", simplify = TRUE).

If using this option the criteria is specific for the id column. The *id* column
must represent the algorithm used, cluster, and number of sequences within the cluster.
For example, an id should look like "kmed_cluster1_n288" where "kmed" represents the
clustering algorithm used, "_cluster1" indicates the pattern came from the first
cluster, and "_n228" indicates that 228 sequences were apart of cluster 1. An example
of how this can be created is:

    formatted_kmed \%>\% mutate(id = paste0("kmed", "_cluster", cluster, "_n", n))
}
}
\value{
Nothing is returned, only exports results.
}
\description{
A function that generates reports and exports the files to the default or specified location;
the default placement is in the user's document folder. This creates a main folder with
two subfolders, public and private, which exports report documents to the corresponding
folder based on data sensitivity.

**Public Folder**

    * Pattern file: Contains the clusters, number of sequences within the cluster (n),
                      the percentage of sequences within that cluster (n_percent), the
                      pattern type (ex. consensus), and the resulting sequence. This also
                      contains a truncated weighted sequence, and the unique events
                      within the sequence.

**Private Folder**

    * Alignments file: contains the sequences and their alignment within each cluster.

    * All Sequences file: Contains the clusters, number of sequences within the cluster (n),
                            the percentage of sequences within that cluster (n_percent), the
                            pattern type (ex. consensus), and the resulting sequence. This also
                            contains a truncated weighted sequence, the full weighted sequence, and
                            the unique events within the sequence.

    * Weighted Sequences file: Contains the clusters, number of sequences within the cluster (n),
                                the percentage of sequences within that cluster (n_percent),
                                and the full weighted sequence of the cluster.
}
\examples{
  data("mvad")

  clustered_kmed <- mvad \%>\%
                          aggregate_sequences(format = "\%Y-\%m-\%d",
                                              unit = "month",
                                              n_units = 1,
                                              summary_stats=FALSE) \%>\%
                          cluster_kmedoids(k = 5)

  patterns_kmed <- clustered_kmed \%>\%  filter_pattern(threshold = 0.5,
                                                        pattern_name = 'consensus')

  patterns_kmed \%>\% generate_reports(end_filename_with = "_kmed")
}
