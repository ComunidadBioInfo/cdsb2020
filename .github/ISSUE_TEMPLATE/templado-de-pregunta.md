---
name: Templado para preguntas
about: Describe tu pregunta
title: "[PREGUNTA] ¿Tu pregunta?"
labels: question
assignees: ''

---

## Contexto

Provee algo de contexto para tu pregunta, por ejemplo:

* el nombre de la presentación, por ejemplo `01-introduction`.
* la liga al código fuente de la presentación, por ejemplo https://github.com/comunidadbioinfo/csdb2020/blob/master/scRNAseq/00-template.Rmd#L24-L28
* la liga a algún `commit` de GitHub, ejemplo: https://github.com/ComunidadBioInfo/cdsb2020/commit/1922523982d3bfe75531d2b327e7b05ca56bd6d4
* la liga a una línea de código dentro de un `commit`, ejemplo: https://github.com/ComunidadBioInfo/cdsb2020/commit/1922523982d3bfe75531d2b327e7b05ca56bd6d4#diff-2232be352ff644d9752095cfdd3c0600R110
* la liga al código de un paquete de R, ejemplo:  https://github.com/LieberInstitute/spatialLIBD/blob/master/R/run_app.R#L51-L55

## Código para tu pregunta

Incluye el código que corriste y tus comentarios

```R
## para crear un mensaje de error
stop('hola')

## para obtener más información de que llevó al error
traceback()
```

## Pequeño ejemplo reproducible

Si copias las líneas de código que te llevaron al error, entonces puedes correr [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html) para crear un pequeño sitio web con el código que puedes facilmente copiar y pegar en este _issue_ de GitHub y el cual le ayudará a otros a que te ayuden.

```R
## para crear un mensaje de error
stop('hola')
#> Error in eval(expr, envir, enclos): hola

## para obtener más información de que llevó al error
traceback()
#> No traceback available
```


## Información de tu sesión de R

Acuérdate de incluir la información completa de tu sesión de R.

```R
options(width = 120)
sessioninfo::session_info()
```