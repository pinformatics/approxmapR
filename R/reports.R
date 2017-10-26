generate_reports <- function(w_sequence_dataframe,
                             html_format = TRUE,
                             truncate_patterns = FALSE){
  stopifnot("W_Sequence_Dataframe" %in% class(w_sequence_dataframe))

  formatted <-
    w_sequence_dataframe %>%
    format_sequence(compare = TRUE,
                    html_format = html_format,
                    truncate_patterns = truncate_patterns)
  # %>%
    # mutate(sequence = ifelse(pattern == "weighted_sequence", str_replace(sequence, "<\\(", " < ( "), sequence))

  rmarkdown::render("inst/rmd_w_sequence.Rmd",
                    params = list(input = formatted,
                                  title = "All Sequences"),
                    output_file = "all_sequences.html")
  patterns <-
    formatted %>%
    filter(pattern != "weighted_sequence")

  rmarkdown::render("inst/rmd_w_sequence.Rmd",
                    params = list(input = patterns,
                                  title = "Patterns"),
                    output_file = "patterns.html")

  w_sequences <-
    formatted %>%
    filter(pattern == "weighted_sequence") %>%
    select(-pattern)

  rmarkdown::render("inst/rmd_w_sequence.Rmd",
                    params = list(input = w_sequences,
                                  title = "Weighted Sequences"),
                    output_file = "weighted_sequences.html")

}
