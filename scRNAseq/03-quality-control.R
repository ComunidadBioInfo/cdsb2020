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
## Datos
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Obtén la posición del cromosoma para cada gen
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
    ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID",
    column = "SEQNAME"
)

# Identifica los genes mitocondriales
is.mito <- which(location == "MT")

library('scater')
sce.416b <- addPerCellQC(sce.416b,
    subsets = list(Mito = is.mito))


## ----qc_metrics, cache=TRUE, dependson='all_code'-----------------------------------------
plotColData(sce.416b, x = "block", y = "detected")

plotColData(sce.416b, x = "block", y = "detected") +
    scale_y_log10()

plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)


## ----all_code_part2, cache = TRUE, dependson='all_code'-----------------------------------
# Valores de límite ejemplo
qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Obtenemos un resumen del número de células
# eliminadas por cada filtro
DataFrame(
    LibSize = sum(qc.lib),
    NExprs = sum(qc.nexprs),
    SpikeProp = sum(qc.spike),
    MitoProp = sum(qc.mito),
    Total = sum(discard)
)

qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")
qc.nexprs2 <- isOutlier(sce.416b$detected, log = TRUE,
    type = "lower")
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher")
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extraemos los límites de valores (thresholds)
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")
# Obtenemos un resumen del número de células
# eliminadas por cada filtro
DataFrame(
    LibSize = sum(qc.lib2),
    NExprs = sum(qc.nexprs2),
    SpikeProp = sum(qc.spike2),
    MitoProp = sum(qc.mito2),
    Total = sum(discard2)
)

## Más pruebas
plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)

batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)
qc.lib3 <- isOutlier(sce.416b$sum,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.nexprs3 <- isOutlier(sce.416b$detected,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher",
    batch = batch)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher",
    batch = batch)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Extraemos los límites de valores (thresholds)
attr(qc.lib3, "thresholds")
attr(qc.nexprs3, "thresholds")

# Obtenemos un resumen del número de células
# eliminadas por cada filtro
DataFrame(
    LibSize = sum(qc.lib3),
    NExprs = sum(qc.nexprs3),
    SpikeProp = sum(qc.spike3),
    MitoProp = sum(qc.mito3),
    Total = sum(discard3)
)


## ----use_case, cache=TRUE, dependson= c('all_code', 'all_code_part2')---------------------
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

plotColData(sce.grun, x = "donor", y = "altexps_ERCC_percent")

discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor)
discard.ercc2 <- isOutlier(
    sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor,
    subset = sce.grun$donor %in% c("D17", "D2", "D7")
)

plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc)
)
plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc2)
)

# Agregamos información sobre que células
# tienen valores extremos
sce.416b$discard <- discard2

# Haz esta gráfica para cada medida de
# control de calidad (QC)
plotColData(
    sce.416b,
    x = "block",
    y = "sum",
    colour_by = "discard",
    other_fields = "phenotype"
) +
    facet_wrap( ~ phenotype) +
    scale_y_log10()

# Otra gráfica de diagnóstico útil
plotColData(
    sce.416b,
    x = "sum",
    y = "subsets_Mito_percent",
    colour_by = "discard",
    other_fields = c("block", "phenotype")
) +
    facet_grid(block ~ phenotype)


## ----use_case_pbmc, cache=TRUE, dependson='all_code'--------------------------------------
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

bcrank <- barcodeRanks(counts(sce.pbmc))

# Mostremos solo los puntos únicos para acelerar
# el proceso de hacer esta gráfica
uniq <- !duplicated(bcrank$rank)
plot(
    bcrank$rank[uniq],
    bcrank$total[uniq],
    log = "xy",
    xlab = "Rank",
    ylab = "Total UMI count",
    cex.lab = 1.2
)
abline(h = metadata(bcrank)$inflection,
    col = "darkgreen",
    lty = 2)
abline(h = metadata(bcrank)$knee,
    col = "dodgerblue",
    lty = 2)
legend(
    "bottomleft",
    legend = c("Inflection", "Knee"),
    col = c("darkgreen", "dodgerblue"),
    lty = 2,
    cex = 1.2
)


set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

# Revisa ?emptyDrops para una explicación de porque hay valores NA
summary(e.out$FDR <= 0.001)

set.seed(100)
limit <- 100
all.out <-
    emptyDrops(counts(sce.pbmc), lower = limit, test.ambient = TRUE)
# Idealmente, este histograma debería verse uniforme.
# Picos grandes cerca de cero indican que los _barcodes_
# con un número total de cuentas menor a "lower" no son
# de origen ambiental.
hist(all.out$PValue[all.out$Total <= limit &
        all.out$Total > 0],
    xlab = "P-value",
    main = "",
    col = "grey80")

sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
sce.pmbc <- addPerCellQC(sce.pbmc, subsets = list(MT = is.mito))
discard.mito <-
    isOutlier(sce.pmbc$subsets_MT_percent, type = "higher")
plot(
    sce.pmbc$sum,
    sce.pmbc$subsets_MT_percent,
    log = "x",
    xlab = "Total count",
    ylab = "Mitochondrial %"
)
abline(h = attr(discard.mito, "thresholds")["higher"], col = "red")


## ----marking, cache=TRUE, dependson='use_case'--------------------------------------------
# Eliminemos las células de calidad baja
# al quedarnos con las columnas del objeto sce que NO
# queremos descartar (eso hace el !)
filtered <- sce.416b[,!discard2]
# Alternativamente, podemos marcar
# las células de baja calidad
marked <- sce.416b
marked$discard <- discard2


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()-----------------------
options(width = 120)
sessioninfo::session_info()

