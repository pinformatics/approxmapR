library(approxmapR)
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)





## Loading in data
data("demo2")





# This data is already in the same index format, i.e. the date represents
#   the number of days since the index date. It has not yet been
#   aggregated into time periods for the event set. 
#
# From here there are two options:
#   (1) Use the -aggregate_sequence()- function and pass -index_date = TRUE-
#   (2) Create a column labelled 'period' which represents the custom grouping
#         that is desired

## Route (1)
agg1 <- demo2 %>% aggregate_sequences(unit = "month",
                                      n_units = 1,
                                      include_date = TRUE,
                                      index_date = TRUE)

## Route (2)
custom_agg <- demo2 %>% 
                group_by(id) %>% 
                mutate(indexed_date = date,
                       n_ndays7 = date / 7,
                       period = as.integer(case_when(date == 0 ~ 1,
                                                     ceiling(n_ndays7) < 5 ~ ceiling(n_ndays7) + 1,
                                                     TRUE ~ floor(n_ndays7) + 2))) %>%
                select(id, date, period, item) %>% arrange(id, date)

## The custom grouping time frame demonstrated in -Route (2)- has the index date (0)
##  be an event set of it's own, the first 4 weeks calculated will each be 
##  their own event set, i.e. every 7 days will be grouped together, and every
##  30 days after the first 4 weeks will be their own event set.


# Have to rename the columns due to the package requirements
names(custom_agg) <- c("id", "date", "period", "event") 


agg2 <- custom_agg %>% pre_aggregated(include_date = TRUE, summary_stats = TRUE)





# Once the data has been based through either -aggregate_sequences()- or
#   -pre_aggregated()-, the data is ready for -pattern_search()-.
#
# If you want to visualize the sequences without running the rest of the algorithm, 
#   one case use -convert_to_sequence()- and each individual will have their
#   events converted to a sequence and displayed in the 'sequence_formatted' column.

agg1 %>% convert_to_sequence()
agg2 %>% convert_to_sequence()

## Example 1. Finding all people with events L and M within the same event set. Keeping true
##              to the methodology where the order of the events within the event set do no matter
##              but the order of the event set within the sequence does, the order in which the events
##              are entered in the event set does not matter while the order of the event sets of the
##              sequence does matter.

## For agg1 this will return 5 individuals
pattern_search(agg1, find_pattern = "(L, M)")
pattern_search(agg1, find_pattern = "(M, L)")

## For agg2 this will return 2 individuals
pattern_search(agg2, find_pattern = "(L, M)")
pattern_search(agg2, find_pattern = "(M, L)")


## Demonstrating the order of the event sets within the sequence matters
pattern_search(agg1, find_pattern = "(I) (L, M)")

pattern_search(agg1, find_pattern = "(L, M) (I)")



