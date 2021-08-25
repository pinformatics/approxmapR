\name{plot_silhouette}
\alias{plot_silhouette}
\title{Plotting Functions}
\usage{
  plot_silhouette <- (ktable,
                      save_graph = TRUE, graph_file_name_individual = NULL, graph_file_name_cluster = NULL,
                      save_table = TRUE, table_file_name = NULL,
                      size_width = 855, size_height = 317,
                      output_directory = "~")
}
\arguments{
\item{ktable}{A silhouette_object from the ktable which is produced from find_optimal_k function. This object is found in the silhouette_object column.}


\item{save_graph}{Default value is TRUE which will save the plot as a PNG file.}

\item{graph_file_name_individual}{Allows user to specify the file name for the graph that is being saved, if nothing is specified then a default file name is used.}

\item{graph_file_name_cluster}{Allows user to specify the file name for the graph that is being saved, if nothing is specified then a default file name is used.}

\item{save_table}{Default value is TRUE which will save the table as a CSV file.}

\item{table_file_name}{Allows user to specify the file name for the table that is being saved, if nothing is specified then a default file name is used.}

\item{size_width}{The width of the graph.}

\item{size_height}{The height of the graph.}

\item{output_directory}{The path to where the exports should be placed.}

}
\value{
Returns two plots, one for the silhouette values at the individual level and another for the silhouette values at the cluster level.
Also has the ability to save the images as well as the silhouette table to a file.
}
\description{
A function which allows one to visualize the silhouette values for the individual as well as for the cluster.
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

  plot_silhouette(ktable$silhouette_object[[8]])

}
