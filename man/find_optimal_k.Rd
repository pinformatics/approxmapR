\name{find_optimal_k}
\alias{find_optimal_k}
\title{Clustering Functions}
\usage{
  find_optimal_k(df_aggregated, clustering = 'k-nn', min_k = 2, max_k = 10, validation_measure = 'silhouette',
                              use_cache = TRUE, save_table = TRUE, file_name = NULL , output_directory = "~")
}
\arguments{
\item{df_aggregated}{A dataframe that has been aggregated using either aggregate_sequences(), or pre_aggregated()
  from approxmapR. This dataframe will have exactly 3 columns: id, period, and event.}

\item{clustering}{The type of clustering algorithm to be used. Currently, K-Nearest Neighbors (clustering = 'k-nn') and
  K-Medoids (clustering = 'k-medoids') are supported.}

\item{mink_k}{The starting K value.}

\item{max_k}{The ending K value.}

\item{validation_measure}{The clustering validation measure to be used. The following measures are
cuurently supported:

    * Average silhouette width (validation_measure = 'silhouette')

    * Dunn index (validation_measure = 'dunn')

    * Average distance within clusters (validation_measure = 'average_within')

    * Average distance between clusters (validation_measure = 'average_between')

    * Ratio of average distance within cluster / average distance between clusters (validation_measure = 'wb_ratio')

    * Half the sum of the within cluster squared dissimilarities divided by the cluster size (validation_measure = 'within_cluster_ss')}

\item{use_cache}{A boolean value to indicate weather or not to use the cached distance matrix.}

\item{save_table}{Default value is TRUE which will save the table as a CSV file.}

\item{file_name}{Allows user to specify the file name for the table that is being saved, if nothing is specified then a default file name is used.}

\item{output_directory}{The path to where the exports should be placed. This creates a folder with the name of "approxmap_results".}

}
\value{
Returns a dataframe with the K value, number of clusters, size of clusters, average silhouette
width and it's 95% confidence interval, Dunn index, average distance within clusters,
average distance between clusters, average distance within clusters / average distance between clusters,
and the sum of average distance within clusters. Additionally, a plot is returned which
show the validation measure of choice and it's corresponding K value.
}
\description{
A function which allows one to find the optimal K value to be used in a supported
clustering algorithm based on a variety of clustering validation measures. Returns a
dataframe with information regarding the clusters and a plot of the stated
validation measure and K values.
}
\examples{
  data("mvad")

  mvad \%>\%
    aggregate_sequences(format = "\%Y-\%m-\%d",
                        unit = "month",
                        n_units = 1,
                        summary_stats = FALSE) \%>\%
    find_optimal_k(clustering = 'k-nn', min_k = 2, max_k = 10,
                   validation_measure = 'silhouette', use_cache = TRUE)
}
