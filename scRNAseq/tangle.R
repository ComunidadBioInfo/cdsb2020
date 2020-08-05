rmds <- dir(here::here("scRNAseq"), pattern = '\\.Rmd$', full.names = TRUE)
sapply(rmds, function(x) {
    knitr::knit(x, output = gsub("Rmd", "R", x), tangle = TRUE)
})
