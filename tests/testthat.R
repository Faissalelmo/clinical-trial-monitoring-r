# Usage: Rscript tests/testthat.R   (from the project folder)
suppressPackageStartupMessages(library(testthat))
for (f in list.files("R", full.names = TRUE)) source(f)
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
