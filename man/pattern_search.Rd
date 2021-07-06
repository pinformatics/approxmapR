\name{pattern_search}
\alias{pattern_search}
\title{pattern_search}
\description{
A function which finds patterns and has the ability to handle regular expression using the
R syntax.  For ease, some shorthands have been created for common pattern structures.
}
\usage{
  pattern_search <- function(Clustered_Dataframe, find_pattern = NULL, event_set = FALSE, exact = FALSE)
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

If specifying the event set structure, by setting event_set = TRUE, the following shorthands
have been created for ease of use.

\item{event*, }{will translate into "([:alnum:]*, )*" and will match an event zero or more times. To be used when specifying any number of events are allowed to occur before a specified event.}
\item{event+, }{will translate into "([:alnum:]*, )+" and will match an event one or more times. To be used when specifying at least 1 event should occur before a specified event.}

\item{, event*}{will translate into "(, [:alnum:]*)*" and will match an event zero or more times. To be used when specifying any number of events are allowed to occur after a specified event. }
\item{, event+}{will translate into "(, [:alnum:]*)+" and will match an event one or more times. To be used when specifying at least 1 event should occur after a specified event.}

\item{**}{will translate into "[[:print:]]*" and will match any printable characters which are: [:alnum:], [:punct:] and space.}

\item{eventset* }{will translate into "(\\([[:alnum:], ]*[[:alnum:]*]+\\) )*" and will match an event set structure zero or more times. Specifies that any number of event sets}
\item{eventset+ }{will translate into "(\\([[:alnum:]*, ]*[[:alnum:]*]+\\) )+" and will match an event set structure one or more times.}

\item{ eventset*}{will translate into "( \\([[:alnum:], ]*[[:alnum:]*]+\\))*" and will match an event set structure zero or more times.}
\item{ eventset+}{will translate into "( \\([[:alnum:]*, ]*[[:alnum:]*]+\\))+" and will match an event set structure one or more times.}
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


# Find all sequences which contain an event set with events K and M;
#   This allows for events to occur before, after, or between events K and M
pattern_search(pats, find_pattern = "(K, M)")
# Returns 4 sequences

# Find all sequences which contain an event set with I and M as well as
# an event set that contains L and M
pattern_search(pats, find_pattern = "(K, M), (L, M)")

# For more granular control over events within an event set, set event_set = TRUE
#   This will allow fine control over the possibility of events to occur before,
#     after, or between the events of interest
}
