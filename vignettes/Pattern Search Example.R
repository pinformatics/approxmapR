library(approxmapR)
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)





# While there are example data sets available within approxmapR, for this example
#   more cases are being created for a larger demonstration
data("demo1")
demo1 <- data.frame(do.call("rbind", strsplit(as.character(demo1$id.date.item), ",")))
names(demo1) <- c("id", "period", "event")



##################################
#  CREATING DATA FOR TEST CASES  #
##################################

#### Adding sequence
## Copying ID = 8 to switch the last 2 events
to_add <- demo1 %>% filter(id == 8)
to_add$id <- 808
to_add[6, 3] <- "M"
to_add[7, 3] <- "L"


demo1 <- rbind(demo1, to_add)

(demo1 %>% filter(id == 808))[, 1:3]




#### Adding sequence
## Copying ID 8 again and removing some event sets
to_add <- demo1 %>% filter(id == 8)
to_add$id <- 8082
to_add <- to_add[-c(5), ]

demo1 <- rbind(demo1, to_add)




#### Adding sequence
## Copying ID = 7 to test specificity of regex; only keeping the first event set
to_add <- (demo1 %>% filter(id == 7))[1:3, 1:3]
to_add$id <- 707

demo1 <- rbind(demo1, to_add)




#### Adding sequence
to_add <- demo1 %>% filter(id == 8)
to_add$id <- "900"


to_add[1:4, 3] <- c("I", "K", "Z", "M")
to_add[3:4, 2] <- c("12/20/2015", "12/21/2015")

to_add[6, 2:3] <- c("7/23/2016", "B")
to_add[7, 2:3] <- c("7/24/2016", "C")

to_add <- rbind(to_add, c("900", "2/03/2017", "Z"))
to_add <- rbind(to_add, c("900", "2/03/2017", "I")) # If multiple events occur on same day then it choose the first one (alphabetically sorted)
to_add <- rbind(to_add, c("900", "2/04/2017", "Z")) # If multiple same events occur in event set, it keeps the first one (time sorted)
to_add <- rbind(to_add, c("900", "2/05/2017", "Z"))
to_add <- rbind(to_add, c("900", "2/05/2017", "J"))
to_add <- rbind(to_add, c("900", "2/06/2017", "Z"))
to_add <- rbind(to_add, c("900", "2/07/2017", "Z"))


demo1 <- rbind(demo1, to_add)




#### Adding multiple sequences
to_add <- data.frame(id = c("901", "901", "901", "901", "901", "901", "901", "901", "901", "901"),
                     period = c("01/01/2017", "01/02/2017", "01/03/2017", "01/04/2017", "01/05/2017", "02/04/2017", "02/05/2017", "02/06/2017", "02/07/2017", "02/08/2017"),
                     event = c("I", "K", "M", "Q", "Z", "B", "C", "D", "G", "M"))


to_add <- rbind(to_add, c("500", "01/01/2017", "A"))
to_add <- rbind(to_add, c("500", "01/02/2017", "B"))
to_add <- rbind(to_add, c("500", "01/03/2017", "C"))
to_add <- rbind(to_add, c("500", "03/04/2017", "B"))
to_add <- rbind(to_add, c("500", "03/05/2017", "I"))
to_add <- rbind(to_add, c("500", "03/06/2017", "J"))
to_add <- rbind(to_add, c("500", "04/07/2017", "X"))
to_add <- rbind(to_add, c("500", "06/01/2017", "B"))
to_add <- rbind(to_add, c("500", "06/02/2017", "K"))
to_add <- rbind(to_add, c("500", "06/03/2017", "L"))
to_add <- rbind(to_add, c("500", "06/03/2017", "M"))


to_add <- rbind(to_add, c("501", "01/01/2017", "A2"))
to_add <- rbind(to_add, c("501", "01/02/2017", "B4"))
to_add <- rbind(to_add, c("501", "01/03/2017", "C5"))
to_add <- rbind(to_add, c("501", "03/04/2017", "B5"))
to_add <- rbind(to_add, c("501", "03/05/2017", "I6"))
to_add <- rbind(to_add, c("501", "03/06/2017", "J6"))
to_add <- rbind(to_add, c("501", "04/07/2017", "X4"))
to_add <- rbind(to_add, c("501", "06/01/2017", "B4"))
to_add <- rbind(to_add, c("501", "06/02/2017", "K1"))
to_add <- rbind(to_add, c("501", "06/03/2017", "L2"))
to_add <- rbind(to_add, c("501", "06/03/2017", "M5"))


# Adding cases that are all numeric
to_add <- rbind(to_add, c("502", "01/01/2017", "2"))
to_add <- rbind(to_add, c("502", "01/02/2017", "4"))
to_add <- rbind(to_add, c("502", "01/03/2017", "5"))
to_add <- rbind(to_add, c("502", "03/04/2017", "5"))
to_add <- rbind(to_add, c("502", "03/05/2017", "6"))
to_add <- rbind(to_add, c("502", "03/06/2017", "6"))
to_add <- rbind(to_add, c("502", "04/07/2017", "4"))
to_add <- rbind(to_add, c("502", "06/01/2017", "4"))
to_add <- rbind(to_add, c("502", "06/02/2017", "1"))
to_add <- rbind(to_add, c("502", "06/03/2017", "2"))
to_add <- rbind(to_add, c("502", "06/03/2017", "5"))


# Adding cases that are a combination of numeric_alpha characters
to_add <- rbind(to_add, c("503", "01/01/2017", "2A"))
to_add <- rbind(to_add, c("503", "01/02/2017", "4B"))
to_add <- rbind(to_add, c("503", "01/03/2017", "5C"))
to_add <- rbind(to_add, c("503", "03/04/2017", "5B"))
to_add <- rbind(to_add, c("503", "03/05/2017", "6I"))
to_add <- rbind(to_add, c("503", "03/06/2017", "6J"))
to_add <- rbind(to_add, c("503", "04/07/2017", "4X"))
to_add <- rbind(to_add, c("503", "06/01/2017", "4B"))
to_add <- rbind(to_add, c("503", "06/02/2017", "1K"))
to_add <- rbind(to_add, c("503", "06/03/2017", "2L"))
to_add <- rbind(to_add, c("503", "06/03/2017", "5M"))


demo1 <- rbind(demo1, to_add) %>% arrange(id, period)





#######################
# CLUSTERING THE DATA #
#######################

agg <- demo1 %>% aggregate_sequences(format = "%m/%d/%Y",
                                     unit = "month",
                                     n_units = 1,
                                     summary_stats = FALSE)


clustered <- agg %>% cluster_kmedoids(k = 2, use_cache = FALSE)


pats <- clustered %>% filter_pattern(threshold = 0.3, pattern_name = "consensus")





####################
# Finding Patterns #
####################

# The documentation can be viewed using
??pattern_search

# To see the sequences use the following
agg %>% convert_to_sequence()



# The default approach allows one to specify the events of interest to occur within
#   an event set of a sequence. The method will then pull all sequences that
#   has those events occurring within the same event set.

# [NOTE] With the default approach the order of the events specified does not matter

# Here we are searching for sequences that have events B and D occurring within the same event set
pattern_search(pats, find_pattern = "(B, D)")
pattern_search(pats, find_pattern = "(D, B)")


# Here we are searching for sequences that have events B4 and C5 occurring within an event set,
#   and events K1 and L2 occurring with an event set - both must occur
#
# Each event set needs to be separated by a space
pattern_search(pats, find_pattern = "(C5, B4) (K1, L2)")





## Finding patterns with more specificity ##

# The default approach is flexible, if one wants to have more granular control
#   then set event_set = TRUE

# [NOTE] Here the order of the events specified matter and one must indicate if
#         events are allowed to be occur before, after, or between events specified

# For example, the following 2 searched will not return the same results
pattern_search(pats, find_pattern = "(I, Z) (B)")
pattern_search(pats, find_pattern = "(I, Z) (B)", event_set = TRUE)
# In order to duplicate the results returned with event_set = TRUE, the correct
#   pattern specification would be:
pattern_search(pats, find_pattern = "(I, event*, Z) (B, event*)", event_set = TRUE)
