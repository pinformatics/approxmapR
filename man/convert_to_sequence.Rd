% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_aggregate.R
\name{convert_to_sequence}
\alias{convert_to_sequence}
\title{Dataframe to sequence}
\usage{
convert_to_sequence(df_seq)
}
\arguments{
\item{df_seq}{The aggregated dataframe that has the sequence ids, itemset ids and events.}
}
\value{
Returns a dataframe with a sequence object
(and a formatted sequence for readability)
corresponding to each id
}
\description{
Dataframe to sequence
}
\examples{
pre_agg_demo \%>\%
pre_aggregated \%>\%
convert_to_sequence()
}
