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
# Usemos datos de pbmc4k
library('BiocFileCache')
bfc <- BiocFileCache()
raw.path <-
    bfcrpath(
        bfc,
        file.path(
            "http://cf.10xgenomics.com/samples",
            "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
        )
    )
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library('DropletUtils')
library('Matrix')
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

# Anotación de los genes
library('scater')
rownames(sce.pbmc) <- uniquifyFeatureNames(rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)
library('EnsDb.Hsapiens.v86')
location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = rowData(sce.pbmc)$ID,
    column = "SEQNAME",
    keytype = "GENEID"
)


## ----all_code2, cache=TRUE, dependson='all_code'------------------------------------------
# Detección de _droplets_ con células
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

# Control de calidad
stats <-
    perCellQCMetrics(sce.pbmc, subsets = list(Mito = which(location == "MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent,
    type = "higher")
sce.pbmc <- sce.pbmc[, !high.mito]

# Normalización de los datos
library('scran')
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)
sce.pbmc <- logNormCounts(sce.pbmc)


## ----all_code3, cache=TRUE, dependson='all_code2'-----------------------------------------
# Set de datos de ejemplo: 416B ------------------------------------------------

library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Anotación de genes
library('AnnotationHub')
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(
    ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID",
    column = "SYMBOL"
)
rowData(sce.416b)$SEQNAME <- mapIds(
    ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID",
    column = "SEQNAME"
)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL,
    rowData(sce.416b)$SYMBOL)

# Control de calidad
mito <- which(rowData(sce.416b)$SEQNAME == "MT")
stats <- perCellQCMetrics(sce.416b, subsets = list(Mt = mito))
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("subsets_Mt_percent", "altexps_ERCC_percent"),
    batch = sce.416b$block
)
sce.416b <- sce.416b[, !qc$discard]

# Normalización
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)


## ----all_code4, cache=TRUE, dependson='all_code3'-----------------------------------------
# Varianza de las log-counts ---------------------------------------------------

dec.pbmc <- modelGeneVar(sce.pbmc)

# Visualicemos la relación entre la media y la varianza
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression")
curve(fit.pbmc$trend(x),
    col = "dodgerblue",
    add = TRUE,
    lwd = 2)

# Ordenemos por los genes más interesantes para checar
# los datos
dec.pbmc[order(dec.pbmc$bio, decreasing = TRUE), ]


## ----all_code5, cache=TRUE, dependson='all_code4'-----------------------------------------
# Coeficiente de variación -----------------------------------------------------
dec.cv2.pbmc <- modelGeneCV2(sce.pbmc)

# Visualicemos la relación con la media
fit.cv2.pbmc <- metadata(dec.cv2.pbmc)
plot(fit.cv2.pbmc$mean, fit.cv2.pbmc$cv2, log = "xy")
curve(fit.cv2.pbmc$trend(x),
    col = "dodgerblue",
    add = TRUE,
    lwd = 2)

# Ordenemos por los genes más interesantes para checar
# los datos
dec.cv2.pbmc[order(dec.cv2.pbmc$ratio, decreasing = TRUE), ]


## ----all_code6, cache=TRUE, dependson='all_code5'-----------------------------------------
# En la presencia de muestras técnicas añadidas (spike-ins) --------------------

dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")
# Ordenemos por los genes más interesantes para checar
# los datos
dec.spike.416b[order(dec.spike.416b$bio, decreasing = TRUE), ]

# In the absence of spike-ins --------------------------------------------------

set.seed(0010101)
dec.pois.pbmc <- modelGeneVarByPoisson(sce.pbmc)
# Ordenemos por los genes más interesantes para checar
# los datos
dec.pois.pbmc[order(dec.pois.pbmc$bio, decreasing = TRUE), ]

# Considerando factores experimentales -----------------------------------------

dec.block.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC",
    block = sce.416b$block)
dec.block.416b[order(dec.block.416b$bio, decreasing = TRUE), ]
dec.block.416b$per.block
dec.block.416b$per.block$X20160113


## ----all_code7, cache=TRUE, dependson='all_code6'-----------------------------------------
# Seleccionando los genes altamente variables (HVG) ----------------------------

# Utiliza modelGeneVar() detrás de cámaras
hvg.pbmc.var <- getTopHVGs(dec.pbmc, n = 1000)
str(hvg.pbmc.var)

# Utiliza modelGeneVarWithSpikes() detrás de cámaras
hvg.416b.var <- getTopHVGs(dec.spike.416b, n = 1000)
str(hvg.416b.var)

# O utiliza modelGeneCV2() al especificar `var.field`
hvg.pbmc.cv2 <- getTopHVGs(dec.cv2.pbmc,
    var.field = "ratio", n = 1000)
str(hvg.pbmc.cv2)


## ----all_code7b, cache=TRUE, dependson='all_code7', warning=FALSE-------------------------
if (!requireNamespace("gplots", quietly = TRUE))
    install.packages("gplots")

if (!requireNamespace("VennDiagram", quietly = TRUE))
    BiocManager::install("VennDiagram")

## Un diagrama de venn rápido pero sencillo
gplots::venn(list('var' = hvg.pbmc.var, 'cv2' = hvg.pbmc.cv2))

## Otro más bonito pero más complejo
v <- VennDiagram::venn.diagram(
    list('var' = hvg.pbmc.var, 'cv2' = hvg.pbmc.cv2),
    filename = NULL,
    fill = c('forest green', 'orange')
    
)
grid::grid.newpage()
grid::grid.draw(v)


## ----all_code8, cache=TRUE, dependson='all_code7'-----------------------------------------
# Seleccionando los HVGs usando significancia estadística ----------------------

# Utiliza modelGeneVar() detrás de cámaras
hvg.pbmc.var.2 <- getTopHVGs(dec.pbmc, fdr.threshold = 0.05)
str(hvg.pbmc.var.2)
# Utiliza modelGeneVarWithSpikes() detrás de cámaras
hvg.416b.var.2 <- getTopHVGs(dec.spike.416b,
    fdr.threshold = 0.05)
str(hvg.416b.var.2)

# O utiliza modelGeneCV2() al especificar `var.field`
hvg.pbmc.cv2.2 <- getTopHVGs(dec.cv2.pbmc,
    var.field = "ratio", fdr.threshold = 0.05)
str(hvg.pbmc.cv2.2)


## ----all_code9, cache=TRUE, dependson='all_code8'-----------------------------------------
# Seleccionando como HVGs a los genes arriba de la curva -----------------------

# Utiliza modelGeneVar() detrás de cámaras
hvg.pbmc.var.3 <- getTopHVGs(dec.pbmc, var.threshold = 0)
str(hvg.pbmc.var.3)
# Utiliza modelGeneVarWithSpikes() detrás de cámaras
hvg.416b.var.3 <- getTopHVGs(dec.spike.416b,
    var.threshold = 0)
str(hvg.416b.var.3)

# O utiliza modelGeneCV2() al especificar `var.field` y
# el valor de `var.threshold`
hvg.pbmc.cv2.3 <- getTopHVGs(dec.cv2.pbmc,
    var.field = "ratio", var.threshold = 1)
str(hvg.pbmc.cv2.2)


## ----all_code10, cache=TRUE, dependson='all_code9'----------------------------------------
# Usando todo al mismo tiempo --------------------------------------------------

dec.pbmc <- modelGeneVar(sce.pbmc)
chosen <- getTopHVGs(dec.pbmc, prop = 0.1)
str(chosen)

# Seleccionando el subconjunto de HVGs -----------------------------------------

sce.pbmc.hvg <- sce.pbmc[chosen, ]
sce.pbmc.hvg

# Especificando los HVGs en funciones posteriores ------------------------------

# Ejemplo de como especificar los HVGs en funciones posteriores
# Realizar el PCA usando solo los HVGs de "chosen"
sce.pbmc <- runPCA(sce.pbmc, subset_row = chosen)
sce.pbmc

# Brujeria ---------------------------------------------------------------------

# Agrega tu SCE principal al SCE del subconjunto de HVGs
altExp(sce.pbmc.hvg, "original") <- sce.pbmc
sce.pbmc.hvg
altExp(sce.pbmc.hvg, "original")


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()-----------------------
options(width = 120)
sessioninfo::session_info()

