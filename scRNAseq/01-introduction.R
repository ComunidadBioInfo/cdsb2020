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


## ----install, eval = FALSE----------------------------------------------------------------
## ## Para instalar paquetes de Bioconductor
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## ## Instala los paquetes de R que necesitamos
## BiocManager::install(
##     c(
##         'SingleCellExperiment',
##         'usethis',
##         'here',
##         'scran',
##         'scater',
##         'scRNAseq',
##         'org.Mm.eg.db',
##         'AnnotationHub',
##         'ExperimentHub',
##         'BiocFileCache',
##         'DropletUtils',
##         'EnsDb.Hsapiens.v86',
##         'TENxPBMCData',
##         'BiocSingular',
##         'batchelor',
##         'uwot',
##         'Rtsne',
##         'pheatmap',
##         'fossil',
##         'ggplot2',
##         'cowplot',
##         'RColorBrewer',
##         'plotly',
##         'iSEE',
##         'pryr',
##         'spatialLIBD',
##         'sessioninfo',
##         'scPipe'
##     )
## )


## ## Si tienes las llaves de SSH configuradas

## git clone git@github.com:comunidadbioinfo/cdsb2020.git

## 
## ## o vía https

## git clone https://github.com/comunidadbioinfo/cdsb2020.git


## ----clone_repo, eval = FALSE-------------------------------------------------------------
## git2r::clone('https://github.com/comunidadbioinfo/cdsb2020',
##     'csdb2020')


## ----proj, eval = FALSE-------------------------------------------------------------------
## usethis::create_project('~/Desktop/cdsb2020_leo')


## ----create_setup, eval = FALSE-----------------------------------------------------------
## ## Crea un archivo de setup
## usethis::use_r('00-setup.R')


## ----use_git, eval = FALSE----------------------------------------------------------------
## ## Crea el repositorio de Git
## usethis::use_git()
## 
## ## Configura tu conexión a GitHub de ser necesario
## usethis::browse_github_token()
## usethis::edit_r_environ() ## y después reinicia R
## 
## ## Utiliza GitHub
## usethis::use_github() ## crea un commit, luego corre este comando
## 
## ## Empieza tus notas sobre la introducción a scRNA-seq
## usethis::use_r('01-introduction.R')


## ----'quick_intro_01', message = FALSE----------------------------------------------------
library('scRNAseq')
library('scater')
library('scran')
library('plotly')


## ----'quick_intro_02', cache = TRUE-------------------------------------------------------
sce <- scRNAseq::MacoskoRetinaData()

## ¿Qué tan grandes son los datos?
pryr::object_size(sce)

## ¿Cómo es el objeto?
sce


## ----'quick_intro_03', cache = TRUE-------------------------------------------------------
# Control de calidad.
es.mito <- grepl("^MT-", rownames(sce))
qcstats <-
    scater::perCellQCMetrics(sce, subsets = list(Mito = es.mito))
filtered <-
    scater::quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalización.
sce <- scater::logNormCounts(sce)

# Selección de genes.
dec <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(dec, prop = 0.1)

# Reducción de dimensiones.
set.seed(1234)
sce <- scater::runPCA(sce, ncomponents = 25, subset_row = hvg)
sce <- scater::runUMAP(sce, dimred = 'PCA', external_neighbors = TRUE)

# Clustering.
g <- scran::buildSNNGraph(sce, use.dimred = 'PCA')
sce$clusters <- factor(igraph::cluster_louvain(g)$membership)


## ----'quick_intro_04'---------------------------------------------------------------------
# Visualización.
scater::plotUMAP(sce, colour_by = "clusters")


## ----'quick_intro_05', eval = FALSE-------------------------------------------------------
## # Visualización interactiva.
## p <- scater::plotUMAP(sce, colour_by = "clusters")
## plotly::ggplotly(p)


## ----'quick_intro_06', eval = FALSE, echo = FALSE-----------------------------------------
## # De https://github.com/rstudio/htmltools/issues/90
## p <- scater::plotUMAP(sce, colour_by = "clusters")
## pi <- plotly::ggplotly(p)
## f <- '01-introduction_files/figure-html/quick_intro_06.html'
## htmlwidgets::saveWidget(pi, here::here(f))
## htmltools::tags$iframe(
##     src = f,
##     width = "100%",
##     height = "400",
##     scrolling = "no",
##     seamless = "seamless",
##     frameBorder = "0"
## )


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()-----------------------
options(width = 120)
sessioninfo::session_info()

