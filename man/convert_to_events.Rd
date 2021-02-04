\name{convert_to_events}
\alias{convert_to_events}
\title{Sequence to dataframe}
\usage{
  convert_to_events(data, id_column, sequence_column)
}
\arguments{
\item{data}{A data frame that has a column which is an identifier variable (ex. id),
and another that represents a sequence with as a character class.}

\item{id_column}{The column which is the identifier variable.}

\item{sequence_column}{The column which contains the sequences as a character class.}
}
\value{
Returns a dataframe with 3 columns:

    * id
    * period
    * event
}
\description{
A function that converts the id and sequence into a dataframe which contains
the id, period, and event.
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

  formatted_kmed <- patterns_kmed \%>\% format_sequence()

  formatted_kmed \%>\% convert_to_events("cluster", "consensus_pattern")
}
