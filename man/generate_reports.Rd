\name{generate_reports}
\alias{generate_reports}
\title{Report generation}
\usage{
  generate_reports(w_sequence_dataframe, html_format = TRUE,
                    output_directory = "~", end_filename_with = "")
}
\arguments{
\item{w_sequence_dataframe}{A dataframe with class "W_Sequence_Dataframe". This will be
the dataframe that resulted from extracting the patterns after clustering.}

\item{html_format}{A boolean value to indicate if the exports should have HTML formatting.}

\item{output_directory}{The path to where the exports should be placed. This creats a
folder with the name of "approxmap_results".}

\item{end_filename_with}{The option of appending to the end of the default file names. This
is useful if running multiple algorithms that will be exported to the same output_directory.}
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

  cclustered_kmed <- mvad \%>\%
                          aggregate_sequences(format = "\%Y-\%m-\%d",
                                              unit = "month",
                                              n_units = 1,
                                              summary_stats=FALSE) \%>\%
                          cluster_kmedoids(k = 5)

  patterns_kmed <- clustered_kmed \%>\%  filter_pattern(threshold = 0.5,
                                                        pattern_name = 'consensus')

  patterns_kmed \%>\% generate_reports(end_filename_with = "_kmed")
}
