\name{pattern_search}
\alias{pattern_search}
\title{pattern_search}
\description{
A function which finds patterns and has the ability to handle regular expression using the
R syntax. For ease, some shorthands have been created for common pattern structures.

    * "event*," will translate into "[[:alnum:], ]*" and will match an event zero or more times.
    * "event+," will translate into "[[:alnum:], ]+" and will match an event one or more times.
    * "**" will translate into "[[:print:]]*" and will match any printable characters which are: [:alnum:], [:punct:] and space.
    * "eventset*, " will translate into "[\\([[:alnum:], ]*[[:alnum:]]+\\), ]*" and will match an event set structure zero or more times.
    * "eventset+, " will translate into "[\\([[:alnum:], ]*[[:alnum:]]+\\), ]*" and will match an event set structure one or more times.
}
\usage{
  Pulls the observations that match an event, event set, sequence, or a list which contains
  a combination of those.
}
\arguments{
  \item{Clustered_Dataframe}{The data frame with class of "Clustered_Dataframe" which the operation will be performed on.}

  \item{find_pattern}{The event, event set, sequence, or a list which contains a combinations of those that are to be found.}

  \item{event_set}{A boolean indicator to indicate if the pattern contains the event set structure of '(event, ..., event)'}

  \item{exact}{a boolean indicator to indicate if the pattern should be searched for how it's exactly entered; this uses stringer 'fixed()' modifier.}
}
\value{
Returns a data frame with the columns of cluster, id, sequence, and sequences which has a
class of ( "tbl_df", "tbl", "data.frame").
}
\examples{
data("demo1")
demo1 <- data.frame(do.call("rbind", strsplit(as.character(demo1$id.date.item), ",")))
names(demo1) <- c("id", "period", "event")


agg <- demo1 %>% aggregate_sequences(format = "%m/%d/%Y",
                                     unit = "month",
                                     n_units = 1,
                                     summary_stats = FALSE)


clustered <- agg %>% cluster_kmedoids(k = 2, use_cache = TRUE)
pats <- clustered %>% filter_pattern(threshold = 0.3, pattern_name = "consensus")


# Find all sequences which contain an event set which has any number of events before the final event L
pattern_search(pats, find_pattern = "(event*, L)")

# Find all sequences which contain any event sets before and after the event set (K, M)
finds <- list("<(I, M), (K), (K, M), (L, M)>", "Z", "(I, J)")
pattern_search(pats, find_pattern = finds)
}
