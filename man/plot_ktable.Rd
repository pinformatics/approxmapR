\name{plot_ktable}
\alias{plot_ktable}
\title{Plotting Functions}
\usage{
  plot_ktable <- (ktable, validation_measure = 'silhouette', save_graph = TRUE, graph_file_name = NULL,
                          size_width = 855, size_height = 317, output_directory = "~")
}
\arguments{
\item{ktable}{A dataframe that has been produced from find_optimal_k which contains various K values and their respective
              cluster validation measures values.}

\item{validation_measure}{The clustering validation measure to be used for plotting. The following measures are
cuurently supported:

    * Average silhouette width (validation_measure = 'silhouette')

    * Dunn index (validation_measure = 'dunn')

    * Average distance within clusters (validation_measure = 'average_within')

    * Average distance between clusters (validation_measure = 'average_between')

    * Ratio of average distance within cluster / average distance between clusters (validation_measure = 'wb_ratio')

    * Half the sum of the within cluster squared dissimilarities divided by the cluster size (validation_measure = 'within_cluster_ss')}

\item{save_graph}{Default value is TRUE which will save the plot as a PNG file.}

\item{graph_file_name}{Allows user to specify the file name for the graph that is being saved, if nothing is specified then a default file name is used.}

\item{size_width}{The width of the graph.}

\item{size_height}{The height of the graph.}

\item{output_directory}{The path to where the exports should be placed.}

}
\value{
Returns a plot with the K value on the x-axis and the cluster validation measure on the y-axis.
}
\description{
A function which allows one to plot a cluster validation measure against the corresponding K value
from the table that is produced from the find_optimal_k function.
}
\examples{
  data("demo1")
  demo1 <- data.frame(do.call("rbind", strsplit(as.character(demo1$id.date.item), ",")))
  names(demo1) <- c("id", "period", "event")


  agg <- demo1 %>% aggregate_sequences(format = "%m/%d/%Y",
                                     unit = "month",
                                     n_units = 1,
                                     include_date = TRUE,
                                     summary_stats = TRUE)

  ktable <- agg %>% find_optimal_k(clustering = "k-medoids", min_k = 2, max_k = 9,
                                   use_cache = TRUE, save_table = TRUE)

  plot_ktable(ktable)

}
