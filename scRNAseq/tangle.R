rmds <- dir(pattern = '\\.Rmd$')
sapply(rmds, knitr::knit, tangle = TRUE)
