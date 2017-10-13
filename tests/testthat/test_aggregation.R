library(approxmapR)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)

context("Aggregate Data")

inp <- read_csv("../../data/inp.csv")

# test_that("Units of time are converted to days properly",{
#   expect_equal(get_n_days("day", 5), 5)
#   expect_equal(get_n_days("week", 5), 35)
#   expect_equal(get_n_days("month", 5), 150)
#   expect_equal(get_n_days("6 months", 5), 900)
#   expect_equal(get_n_days("year", 5), 5*365)
# })


test_that("Convert to sequence returns a list", {
  expect_equal(typeof(inp %>% aggregate_sequences(format = "%Y-%m-%d") %>%
                      convert_to_sequence() %>%
                      pull(sequence)
                      ), "list")
})


inp2 <- inp %>%
  aggregate_sequences(format = "%Y-%m-%d") %>%
  convert_to_sequence() %>%
  pull(sequence)

inp3 <- inp %>%
  aggregate_sequences(format = "%Y-%m-%d", multiset = T) %>%
  convert_to_sequence() %>%
  pull(sequence)

test_that("Mutliset works", {
  expect_equal(inp2[[1]][[2]] %>% length(), 2)
  expect_equal(inp3[[1]][[2]] %>% length(), 3)
})
