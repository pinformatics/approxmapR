
library(approxmapR)
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)





data("mvad")


# Step 1. Clustering - Must specify K value
clustered <- mvad %>%
                  aggregate_sequences(format = "%Y-%m-%d",
                                      unit = "month",
                                      n_units = 1,
                                      summary_stats=FALSE) %>%
                  cluster_knn(k = 5, use_cache = TRUE)



# Step 2. Filter consensus pattern based on threshold
patterns <- clustered %>% filter_pattern(threshold = 0.3, pattern_name = "consensus")



# Step 3. Export to specified folder
patterns %>% generate_reports(end_filename_with = "_knn",
                              output_directory = "C:\Users\...\Documents")
## -output_directory- does not need to be specified, the default location above
