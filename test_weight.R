library(tidyverse)
data("demo1")
demo1 <- data.frame(do.call("rbind", strsplit(as.character(demo1$id.date.item), ",")))
names(demo1) <- c("id", "period", "event")


a <- demo1 %>%
  aggregate_sequences(format = "%m/%d/%Y",
                      unit = "month",
                      n_units = 1,
                      summary_stats = FALSE)

b <- a %>%
  cluster_knn(k = 5, use_cache = TRUE)



###############################################
w <- read.csv("data/w.csv")
w <- unname(data.matrix(w))
w <- w[1:14, 2:15]
w
c <- a %>%
  cluster_knn(k = 5, use_cache = TRUE, weight = w)













