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
sce.416b <- LunSpikeInData(which = "416b")

# Carga el paquete SingleCellExperiment
library('SingleCellExperiment')
# Extrae la matriz de cuentas del set de datos de 416b
counts.416b <- counts(sce.416b)
# Construye un nuevo SCE de la matriz de cuentas
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Revisa el objeto que acabamos de crear
sce

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Accesa la matriz de cuenta del compartimento (slot) "assays" 
# OJO: ¡esto puede inundar tu sesión de R!

# 1. El método general
assay(sce, "counts")[1:6, 1:3]
# 2. El método específico para accesar la matriz de cuentas "counts"
counts(sce)[1:6, 1:3]

sce <- scater::logNormCounts(sce)
# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# 1. El método general
assay(sce, "logcounts")[1:6, 1:3]
# 2. El método específico para accesar la matriz de cuentas
# transformadas "logcounts"
logcounts(sce)[1:6, 1:3]

# Asigna una nueva matriz al compartimento (slot) de "assays"
assay(sce, "counts_100") <- assay(sce, "counts") + 100
# Enumera los "assays" en el objeto
assays(sce)
assayNames(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Extrae la información de las muestras (metadata) del set de datos de 416b
colData.416b <- colData(sce.416b)
# Agrega algo de esa información a nuestro objeto de SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Revisa el objeto que acabamos de actualizar
sce
# Accesa a la información de las muestras (metadata) en nuestro SCE
colData(sce)
# Accesa una columna específica de la información de las muestras (metadata)
table(sce$block)

# Ejemplo de una función que agrega columnas nuevas al colData
sce <- scater::addPerCellQC(sce.416b)
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
colData(sce)

# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

## Agrega las cuentas normalizadas (lognorm) de nuevo
sce <- scater::logNormCounts(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Ejemplo: obtén el subconjunto de células de tipo "wild type"
# Acuérdate que las células son columnas del SCE
sce[, sce$phenotype == "wild type phenotype"]

# Accesa la información de los genes de nuestro SCE
# ¡Está vació actualmente!
rowData(sce)

# Ejemplo de una función que agrega campos nuevos en el rowData
sce <- scater::addPerFeatureQC(sce)
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Obtén la posición del cromosoma para cada gen
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
rowData(sce)$chromosome <- chromosome

# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Ejemplo: obtén el subconjunto de datos donde los genes están en el
# cromosoma 3
# NOTA: which() fue necesario para lidear con los nombres de cromosoma que son
# NA
sce[which(rowData(sce)$chromosome == "3"), ]

# Accesa la información de nuestro experimento usando metadata()
# ¡Está vació actualmente!
metadata(sce)

# La información en el metadata() es como Vegas - todo se vale
metadata(sce) <- list(favourite_genes = c("Shh", "Nck1", "Diablo"),
    analyst = c("Pete"))

# Accesa la información de nuestro experimento usando metadata() de
# nuestro objeto actualizado
metadata(sce)

# Ejemplo: agrega los componentes principales (PCs) de las logcounts
# NOTA: aprenderemos más sobre análisis de componentes principales (PCA) después
sce <- scater::runPCA(sce)
# Revisa el objeto que acabamos de actualizar
sce
# Accesa la matriz de PCA del componente (slot) reducedDims
reducedDim(sce, "PCA")[1:6, 1:3]

# Ejemplo, agrega una representación de los logcounts en t-SNE
# NOTA: aprenderemos más sobre t-SNE después
sce <- scater::runTSNE(sce)
# Revisa el objeto que acabamos de actualizar
sce
# Accesa a la matriz de t-SNE en el componente (slot) de reducedDims
head(reducedDim(sce, "TSNE"))

# Ejemplo: agrega una representación 'manual' de los logcounts en UMAP
# NOTA: aprenderemos más sobre UMAP después y de una forma más sencilla de
#       calcularla
u <- uwot::umap(t(logcounts(sce)), n_components = 2)

# Agrega la matriz de UMAP al componente (slot) reducedDims
reducedDim(sce, "UMAP") <- u

# Accesa a la matriz de UMAP desde el componente (slot) reducedDims
head(reducedDim(sce, "UMAP"))

# Enumera los resultados de reducción de dimensiones en nuestro objeto SCE
reducedDims(sce)

# Extrae la información de ERCC de nuestro SCE para el set de datos de 416b
ercc.sce.416b <- altExp(sce.416b, "ERCC")
# Inspecciona el SCE para los datos de ERCC
ercc.sce.416b

# Agrega el SCE de ERCC como un experimento alternativo a nuestro SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

# Enumera los experimentos alternativos almacenados en nuestro objeto
altExps(sce)

# El crear un subconjunto del SCE por muestra (célula) automáticamente
# obtiene el subconjunto de los experimentos alternativos
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.subset)

# Extrae los factores de tamaño (size factors)
# Estos fueron añadidos a nuestro objeto cuando corrimos
# scater::logNormCounts(sce)
head(sizeFactors(sce))

# "Automáticamente" reemplaza los factores de tamaño
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

# "Manualmente" reemplaza los factores de tamaño
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))


## ----ercc_exercise, cache = TRUE, dependson='all_code'------------------------------------
## Lee los datos de ERCC de la red
ercc_info <-
    read.delim(
        'https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt',
        as.is = TRUE,
        row.names = 2,
        check.names = FALSE
    )

## Pon los datos de ERCC en el mismo orden
m <- match(rownames(altExp(sce, "ERCC")), rownames(ercc_info))
ercc_info <- ercc_info[m, ]

## Normaliza las cuentas de ERCC
altExp(sce, "ERCC") <- scater::logNormCounts(altExp(sce, "ERCC"))


## ----ercc_solution_plots, cache = TRUE, dependson='ercc_exercise'-------------------------
for (i in seq_len(2)) {
    plot(
        log2(10 * ercc_info[, "concentration in Mix 1 (attomoles/ul)"] + 1) ~
            log2(counts(altExp(sce, "ERCC"))[, i] +
                    1),
        xlab = "cuentas log norm",
        ylab = "Mezcla 1: log2(10 * Concentración + 1)",
        main = colnames(altExp(sce, "ERCC"))[i],
        xlim = c(min(logcounts(
            altExp(sce, "ERCC")
        )), max(logcounts(
            altExp(sce, "ERCC")
        )))
    )
    abline(0, 1, lty = 2, col = 'red')
}


## ----all_code_part2, cache=TRUE-----------------------------------------------------------
# Descarga datos de ejemplo procesados con CellRanges
# Paréntesis: al usar BiocFileCache solo tenemos que descargar
#             los datos una vez.
library('BiocFileCache')
bfc <- BiocFileCache()
pbmc.url <-
    paste0(
        "http://cf.10xgenomics.com/samples/cell-vdj/",
        "3.1.0/vdj_v1_hs_pbmc3/",
        "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
    )
pbmc.data <- bfcrpath(bfc, pbmc.url)

# Extrae los archivos en un directorio temporal
untar(pbmc.data, exdir = tempdir())

# Enumera los archivos que descargamos y que extrajimos
# Estos son los archivos típicos de CellRanger
pbmc.dir <- file.path(tempdir(),
    "filtered_feature_bc_matrix")
list.files(pbmc.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library('DropletUtils')
sce.pbmc <- read10xCounts(pbmc.dir)
# Revisa el objeto que acabamos de construir
sce.pbmc

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.pbmc)

# Almacena la información de CITE-seq como un experimento alternativo
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Revisa el objeto que acabamos de actualizar
sce.pbmc

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.pbmc)

# Descarga datos de ejemplo procesados con scPipe
library('BiocFileCache')
bfc <- BiocFileCache()
sis_seq.url <-
    "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

# Extrae los archivos en un directorio temporal
unzip(sis_seq.data, exdir = tempdir())

# Enumera (algunos de) los archivos que descargamos y extrajimos
# Estos son los archivos típicos de scPipe
sis_seq.dir <- file.path(tempdir(),
    "SIS-seq_script-master",
    "data",
    "BcorKO_scRNAseq",
    "RPI10")
list.files(sis_seq.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library('scPipe')
sce.sis_seq <- create_sce_by_dir(sis_seq.dir)
# Revisa el objeto que acabamos de construir
sce.sis_seq

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.sis_seq)

# Descarga un ejemplo de un montón de archivos
library('BiocFileCache')
bfc <- BiocFileCache()
lun_counts.url <-
    paste0(
        "https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
    )
lun_counts.data <- bfcrpath(bfc, lun_counts.url)
lun_coldata.url <-
    paste0("https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.sdrf.txt")
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

# Extrae los archivos en un directorio temporal
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir = lun_counts.dir)

# Enumera los archivos que descargamos y extrajimos
list.files(lun_counts.dir)

# Lee la matriz de cuentas (para una placa)
lun.counts <- read.delim(
    file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
# Almacena la información de la longitud de los genes para después
gene.lengths <- lun.counts$Length
# Convierte los datos de cuentas de genez a una matriz
lun.counts <- as.matrix(lun.counts[, -1])

# Lee la información de las muestras (células)
lun.coldata <- read.delim(lun_coldata.data,
    check.names = FALSE,
    stringsAsFactors = FALSE)
library('S4Vectors')
lun.coldata <- as(lun.coldata, "DataFrame")

# Pon en orden la información de las muestras para que
# sea idéntico al orden en la matriz de cuentas
m <- match(colnames(lun.counts),
    lun.coldata$`Source Name`)
lun.coldata <- lun.coldata[m,]

# Construye la tabla de información de los genes
lun.rowdata <- DataFrame(Length = gene.lengths)

# Construye el objeto de SingleCellExperiment
lun.sce <- SingleCellExperiment(
    assays = list(assays = lun.counts),
    colData = lun.coldata,
    rowData = lun.rowdata
)
# Revisa el objeto que acabamos de construir
lun.sce

## ¿Qué tan grande es el objeto de R?
pryr::object_size(lun.sce)


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()-----------------------
options(width = 120)
sessioninfo::session_info()

