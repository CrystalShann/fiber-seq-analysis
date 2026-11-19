#!/usr/bin/env Rscript
# Execute explicit PDF/TSV-producing chunks without retaining a knitted report.
main <- function() {
  input <- '/project/spott/cshan/fiber-seq/code/entropy/entropy_features.rmd'
  output <- tempfile('entropy_features_', fileext = '.md')
  on.exit(unlink(output), add = TRUE)
  knitr::knit(input, output = output, quiet = TRUE, envir = new.env(parent = globalenv()))
}
main()
