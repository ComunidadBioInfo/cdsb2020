library("readxl")
library("here")
library("sessioninfo")

regis <- read_xlsx(here("registrados", "RegistroCDSB_TIB2020.xlsx"), skip = 3)
conf <- subset(regis, Respuesta == "Confirmado ")

emails <- with(conf, paste0(`Nombre(s)`, " ", Apellidos, " <", `Correo electrónico`, ">"))
emails <- gsub(" ", "", emails)

paste(emails, collapse = ", ")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.0.0 (2020-04-24)
#  os       macOS Catalina 10.15.5
#  system   x86_64, darwin17.0
#  ui       AQUA
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2020-07-09
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version  date       lib source
#  assertthat    0.2.1    2019-03-21 [1] CRAN (R 4.0.0)
#  backports     1.1.6    2020-04-05 [1] CRAN (R 4.0.0)
#  callr         3.4.3    2020-03-28 [1] CRAN (R 4.0.0)
#  cellranger    1.1.0    2016-07-27 [1] CRAN (R 4.0.0)
#  cli           2.0.2    2020-02-28 [1] CRAN (R 4.0.0)
#  crayon        1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
#  desc          1.2.0    2018-05-01 [1] CRAN (R 4.0.0)
#  devtools    * 2.3.0    2020-04-10 [1] CRAN (R 4.0.0)
#  digest        0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
#  dplyr         0.8.5    2020-03-07 [1] CRAN (R 4.0.0)
#  ellipsis      0.3.1    2020-05-15 [1] CRAN (R 4.0.0)
#  fansi         0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
#  fs            1.4.1    2020-04-04 [1] CRAN (R 4.0.0)
#  glue          1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
#  here        * 0.1      2017-05-28 [1] CRAN (R 4.0.0)
#  lifecycle     0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
#  magrittr      1.5      2014-11-22 [1] CRAN (R 4.0.0)
#  memoise       1.1.0    2017-04-21 [1] CRAN (R 4.0.0)
#  pillar        1.4.4    2020-05-05 [1] CRAN (R 4.0.0)
#  pkgbuild      1.0.8    2020-05-07 [1] CRAN (R 4.0.0)
#  pkgconfig     2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
#  pkgload       1.0.2    2018-10-29 [1] CRAN (R 4.0.0)
#  prettyunits   1.1.1    2020-01-24 [1] CRAN (R 4.0.0)
#  processx      3.4.2    2020-02-09 [1] CRAN (R 4.0.0)
#  ps            1.3.3    2020-05-08 [1] CRAN (R 4.0.0)
#  purrr         0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
#  R6            2.4.1    2019-11-12 [1] CRAN (R 4.0.0)
#  Rcpp          1.0.4.12 2020-06-26 [1] Github (RcppCore/Rcpp@653b4ae)
#  readxl      * 1.3.1    2019-03-13 [1] CRAN (R 4.0.0)
#  remotes       2.1.1    2020-02-15 [1] CRAN (R 4.0.0)
#  rlang         0.4.6    2020-05-02 [1] CRAN (R 4.0.0)
#  rprojroot     1.3-2    2018-01-03 [1] CRAN (R 4.0.0)
#  sessioninfo * 1.1.1    2018-11-05 [1] CRAN (R 4.0.0)
#  testthat    * 2.3.2    2020-03-02 [1] CRAN (R 4.0.0)
#  tibble        3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
#  tidyselect    1.0.0    2020-01-27 [1] CRAN (R 4.0.0)
#  usethis     * 1.6.1    2020-04-29 [1] CRAN (R 4.0.0)
#  utf8          1.1.4    2018-05-24 [1] CRAN (R 4.0.0)
#  vctrs         0.3.1    2020-06-05 [1] CRAN (R 4.0.0)
#  withr         2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.0/Resources/library
