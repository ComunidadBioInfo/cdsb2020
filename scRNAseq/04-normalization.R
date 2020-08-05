## ----setup, include=FALSE-----------------------------------------------------------------
options(htmltools.dir.version = FALSE)


## ----xaringan-themer, include=FALSE-------------------------------------------------------
library(xaringanthemer)
solarized_dark(
  code_font_family = "Fira Code",
  code_font_url    = "https://cdn.rawgit.com/tonsky/FiraCode/1.204/distr/fira_code.css"
)


## /* From https://github.com/yihui/xaringan/issues/147  */

## .scroll-output {

##   height: 80%;

##   overflow-y: scroll;

## }

## 
## /* https://stackoverflow.com/questions/50919104/horizontally-scrollable-output-on-xaringan-slides */

## pre {

##   max-width: 100%;

##   overflow-x: scroll;

## }

## 
## /* From https://github.com/yihui/xaringan/wiki/Font-Size */

## .tiny{

##   font-size: 40%

## }

## 
## /* From https://github.com/yihui/xaringan/wiki/Title-slide */

## .title-slide {

##   background-image: url(https://raw.githubusercontent.com/Bioconductor/OrchestratingSingleCellAnalysis/master/images/Workflow.png);

##   background-size: 33%;

##   background-position: 0% 100%

## }


## ----all_code, cache=TRUE-----------------------------------------------------------------
library('scRNAseq')
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)

# Control de calidad
library('scater')
is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
    quickPerCellQC(stats,
        percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[, !qc$discard]


## ----all_code2, cache=TRUE, dependson='all_code'------------------------------------------
# Estimar tamaños de librerías
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)

# Examina la distribución de los tamaños de librerías
# que acabamos de estimar
summary(lib.sf.zeisel)
hist(log10(lib.sf.zeisel), xlab = "Log10[Size factor]", col = "grey80")
ls.zeisel <- colSums(counts(sce.zeisel))
plot(
    ls.zeisel,
    lib.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)


## ----exercise_solution, cache=TRUE, dependson='all_code'----------------------------------
## Primero calcula las sumas
zeisel_sums <- colSums(counts(sce.zeisel))
identical(zeisel_sums, ls.zeisel)

## Ahora asegurate que su media sea 1 (unity mean)
zeisel_size_factors <- zeisel_sums/mean(zeisel_sums)
identical(zeisel_size_factors, lib.sf.zeisel)


## ----all_code3, cache=TRUE, dependson='all_code2'-----------------------------------------
# Normalización por circunvolución (convolution)

library('scran')
# Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)
# Calcula factores de tamaño para la decircunvolución (deconvolution)
deconv.sf.zeisel <-
    calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)

# Examina la distribución de los factores de tamaño
summary(deconv.sf.zeisel)
hist(log10(deconv.sf.zeisel), xlab = "Log10[Factor de tamaño]",
    col = "grey80")
plot(
    ls.zeisel,
    deconv.sf.zeisel,
    log = "xy",
    xlab = "Factores tamaño de librería",
    ylab = "Factores tamaño deconv"
)


## ----all_code4, cache=TRUE, dependson='all_code3'-----------------------------------------
# Factores de tamaño de librería vs 
# factor de tamaño vía circunvolución (convolution)

# Usa colores para los típos celulares que nos dieron
plot(
    lib.sf.zeisel,
    deconv.sf.zeisel,
    xlab = "Factores tamaño de librería",
    ylab = "Factores tamaño deconv",
    log = 'xy',
    pch = 16,
    col = as.integer(factor(sce.zeisel$level1class))
)
abline(a = 0, b = 1, col = "red")


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()-----------------------
options(width = 120)
sessioninfo::session_info()

