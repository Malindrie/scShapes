
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scShapes

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/Malindrie/scShapes.svg?branch=master)](https://travis-ci.com/Malindrie/scShapes)

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/Malindrie/scShapes?branch=master&svg=true)](https://ci.appveyor.com/project/Malindrie/scShapes)
<!-- badges: end -->

We present a novel statistical framework for identifying differential
distributions in single-cell RNA-sequencing (scRNA-seq) data between
treatment conditions by modeling gene expression read counts using
generalized linear models (GLMs). We model each gene independently under
each treatment condition using the error distributions Poisson (P),
Negative Binomial (NB), Zero-inflated Poisson (ZIP) and Zero-inflated
Negative Binomial (ZINB) with log link function and model based
normalization for differences in sequencing depth. Since all four
distributions considered in our framework belong to the same family of
distributions, we first perform a Kolmogorov-Smirnov (KS) test to select
genes belonging to the family of ZINB distributions. Genes passing the
KS test will be then modeled using GLMs. Model selection is done by
calculating the Bayesian Information Criterion and likelihood ratio test
statistic.

While most methods for differential gene expression analysis aim to
detect a shift in the mean of expressed values, single cell data are
driven by over-dispersion and dropouts requiring statistical
distributions that can handle the excess zeros. By modeling gene
expression distributions, our framework can identify subtle variations
that do not involve the change in mean. It also has the flexibility to
adjust for covariates and perform multiple comparisons while explicitly
modeling the variability between samples.

## Installation

You can install the released version of scShapes with:

``` r
devtools::install_github('Malindrie/scShapes')
```

``` r
library(scShapes)
version()
#>             ____
#>           |     |
#> ____ ____ |____ |____ ____ ____ ____  ____
#> [__  |        | |   | ___| |   ||___] [__
#> ___] |___ ____| |   ||___| |___||___  ___]
#>                            |
#>                            |
#>                             Ver:0.1.0
```

## Example

This is a basic example which shows how you can use scShapes for
identifying differential distributions in single-cell RNA-seq data. For
this example data we use the human immune cells (PBMC) dataset
distributed through the
[SeuratData](https://github.com/satijalab/seurat-data) package.

``` r
library(scShapes)

#Loading and preparing data for input 
library(Seurat)
library(SeuratData)

InstallData("ifnb")
LoadData("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
```

We first filter the genes to keep only genes expressed in at least 10%
of cells:

``` r

#First extract the RNA-seq counts from the 'RNA' assay of the seurat object
ifnb.obj <- lapply(ifnb.list, function (x) as.matrix(x@assays$RNA@counts))
ifnb.filtered <- lapply(ifnb.obj, function (x) filter_counts(x, perc.zero = 0.1))
#> Removing 527 rows of genes with all zero counts
#> Removing 778 rows of genes with all zero counts
```

In order to normalize for differences in sequencing depth, the log of
the total UMI counts assigned per cell will be used as an offset in the
GLM. This function is inbuilt in the algorithm; however the user is
required to input the library sizes. We can calculate the library sizes
for the two treatment conditions as;

``` r
ifnb.lib.size <- lapply(ifnb.filtered, function (x) apply(x,2, function(y) sum(y)))
```

The ‘meta.data’ slot of the Seurat object also contains information on
the cell-types, which will be used as a covariate in the GLM model to
account for known biological variation in the data.

``` r
ifnb.variables <- lapply(ifnb.list, function (x) data.frame(
                        cell.type = factor(x@meta.data$seurat_annotations),
                        row.names = colnames(x@assays$RNA)))
```

For the purpose of this example we only run the pipeline for randomly
selected 20 common genes under both treatment conditions ‘CTRL’ and
‘STIM’.

``` r

#Randomly select 20 genes among common genes between the two treatment conditions
comm.genes <- intersect(rownames(ifnb.filtered$CTRL), rownames(ifnb.filtered$STIM))
comm.20.genes <- sample(comm.genes, 20, replace = FALSE)

#Subset the randomly selected 20 genes
ifnb.ctrl <- ifnb.filtered$CTRL[rownames(ifnb.filtered$CTRL) %in% comm.20.genes,]
ifnb.stim <- ifnb.filtered$STIM[rownames(ifnb.filtered$STIM) %in% comm.20.genes,]
```
