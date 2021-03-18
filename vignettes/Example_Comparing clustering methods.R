
library(approxmapR)
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)





data("mvad")


# Step 1. Aggregation
agg <- mvad %>% aggregate_sequences(format = "%Y-%m-%d",
                                    unit = "month",
                                    n_units = 1,
                                    summary_stats=FALSE)



# Step 2. Clustering - Using K-NN and K-Medoids; must specify K value
clustered_knn <- agg %>% cluster_knn(k = 5, use_cache = TRUE)

clustered_kmed <- agg %>% cluster_kmedoids(k = 5, use_cache = TRUE)



# Step 3. Pull the pattern
patterns_knn <- clustered_knn %>% filter_pattern(threshold = 0.3, pattern_name = "consensus")

patterns_kmed <- clustered_kmed %>% filter_pattern(threshold = 0.3, pattern_name = "consensus")



# Step 4. Generate Reports based on initial clustering approach
patterns_knn %>% generate_reports(end_filename_with = "_knn",
                                  output_directory = "~")

patterns_kmed %>% generate_reports(end_filename_with = "_kmed",
                                   output_directory = "~")
## -output_directory- does not need to be specified, the default location above



# Step 5. Format the patterns and create data frame which is the combination of the formatted patterns
formatted_knn <- patterns_knn %>% format_sequence()

formatted_kmed <- patterns_kmed %>% format_sequence()


# Creating a unique ID for each cluster and consensus patterns based on the
#   algorithm and cluster; id = algorithm name (knn or kmed) + cluster # + n
#  
#   Use "_" so the output for the algorithm comparison can parse the cluster and
#     n information into separate columns
knn_output <- formatted_knn %>% 
                            mutate(id = paste0("knn", "_cluster", cluster, "_n", n)) %>%
                            select(id, consensus_pattern)

kmed_output <- formatted_kmed %>% 
                              mutate(id = paste0("kmed", "_cluster", cluster, "_n", n)) %>%
                              select(id, consensus_pattern)



# This combines the seperate consensus patterns into a single data.frame and
#   unnests the consensus patterns into a long format (like original data.frame)
df_consensus <- rbind(kmed_output, knn_output) %>% 
                                               convert_to_events("id", "consensus_pattern")



##### Now, one repeats the above steps #####


# Step 6. Clustering - Use K-Medoids for algorithm comparison; must specify K value
clustered_consensus <- df_consensus %>% 
                                    pre_aggregated() %>% 
                                    cluster_kmedoids(4, use_cache = FALSE)



# Step 7. Pull the pattern
patterns_consensus <- clustered_consensus %>% 
                                          filter_pattern(threshold = 0.3,
                                                         pattern_name = "consensus")



# Step 8. Generate Reports based on initial clustering approach
patterns_consensus %>% generate_reports(end_filename_with = "_comparison", 
                                        output_directory = "~",
                                        algorithm_comparison = TRUE)
## -output_directory- does not need to be specified, the default location above
## Specifying -algorithm_comparison = TRUE- indicates the alignment file
##  should attempt to parse out the "_cluster" and "_n" to seperate columns
##  withint the CSV file.
