\name{validation_measures}
\alias{validation_measures}
\title{Validation Measures}
\usage{
  Cluster validation measures supported in find_optimal_k()
}
\arguments{
  \item{Average Distance Between Clusters}{The average distance between clusters (FPC Documentation)}


  \item{Average Distance Within Clusters}{The average distance within clusters (reweighted so that every observation, rather than every distance, has the same weight) (FPC Documentation)}


  \item{Average Distance Within Clusters / Average Distance Between Clusters}{Average Distance Within Clusters /Average Distance Between (FPC Documentation)}


  \item{Within Cluster Sum of Squares}{A generalization of the within clusters sum of squares (k-means objective function), which is obtained if d is a Euclidean distance matrix. For general distance measures,
                                        this is half the sum of the within cluster squared dissimilarities by the cluster size (FPC Documentation)}


  \item{Dunn Index}{
                      Is equal to the Minimum Separation / Maximum Diameter. Where the minimum separation is the minimum distance of a point in the cluster to a point of another cluster
                      and the maximum diameter is the maximum within cluster distance (FPC Documentation). High values on the Dunn Index indicate the presence of compact a well-separated clusters.
                   }


  \item{Silhouette Index}{
                            “The silhouette width of a point measures the proximity of the point to its own cluster relative to the proximity to other clusters” (Halkidi, Batistaki, & Vazirgiannis, 2001).
                            The silhouette index is the average silhouette width of clusters and ranges from -1 to 1. A value closer to -1 indicates that the point is on average closer to another cluster
                            than the one to which it belongs, while a value closer to 1 indicates the point's average distance to its own cluster is smaller than any other cluster. A higher silhouette index
                            value indicates the presence of compact and separated clusters.
                         }


}
\references{
    1. Dalton, L., Ballarin, V., & Brun, M. (2009). Clustering algorithms: on learning, validation, performance, and applications to genomics. Current genomics, 10(6), 430–445. https://doi.org/10.2174/138920209789177601

    2. Halkidi, M., Batistakis, Y., and Vazirgiannis, M. (2001). On clustering validation techniques. Journal of Intelligent Information Systems. 17(2:3), 107 – 145.

    3. Rousseeuw, P. (1987), Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. Journal of Computational and Applied Mathematics, 20, 53-65.

    4. https://www.rdocumentation.org/packages/fpc/versions/2.2-9/topics/cluster.stats
}
