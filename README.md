
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

You can view the preprint of our method in [bioRxiv](/https://www.biorxiv.org/content/10.1101/2022.02.13.480299v1)

## Installation

You can install the released version of scShapes with:

``` r
devtools::install_github('Malindrie/scShapes')
```

``` r
library(scShapes)
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
library(dplyr)
library(BiocParallel)
set.seed(0xBEEF)

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
ifnb.subset <- list(CTRL = ifnb.ctrl, STIM = ifnb.stim)
```

Perform Kolmogorov-Smirnov test to select genes belonging to the family
of ZINB distributions.

``` r
ifnb.ctrl.KS <- ks_test(ifnb.subset$CTRL, cexpr=ifnb.variables$CTRL, lib.size=ifnb.lib.size$CTRL, BPPARAM=SnowParam(workers=8,type="SOCK"))
ifnb.stim.KS <- ks_test(ifnb.subset$STIM, cexpr=ifnb.variables$STIM, lib.size=ifnb.lib.size$STIM, BPPARAM=SnowParam(workers=8,type="SOCK"))

#Select genes significant from the KS test.
#By default the 'ks_sig' function performs Benjamini-Hochberg correction for multiple hypothese testing
#and selects genes significant at p-value of 0.01

ifnb.ctrl.sig.KS <- ks_sig(ifnb.ctrl.KS)
ifnb.stim.sig.KS <- ks_sig(ifnb.stim.KS)

#Subset UMI counts corresponding to the genes significant from the KS test
ifnb.sig.genes <- list(CTRL = as.data.frame(ifnb.ctrl.sig.KS$genes),
                       STIM = as.data.frame(ifnb.stim.sig.KS$genes))
ifnb.ctrl.KS <- ifnb.filtered$CTRL[rownames(ifnb.filtered$CTRL) %in% rownames(ifnb.sig.genes$CTRL),]
  ifnb.stim.KS <- ifnb.filtered$STIM[rownames(ifnb.filtered$STIM) %in% rownames(ifnb.sig.genes$STIM),]
```

Fit the 4 distributions P,NB,ZIP,ZINB for genes that belong to the ZINB
family of distributions by fitting GLM with log of the library sizes as
an offset and cell types as a covariate in the GLM.

``` r
ifnb.ctrl.fit <- fit_models(counts=ifnb.ctrl.KS, cexpr=ifnb.variables$CTRL, lib.size=ifnb.lib.size$CTRL, BPPARAM=SnowParam(workers=2,type="SOCK"))
ifnb.stim.fit <- fit_models(counts=ifnb.stim.KS, cexpr=ifnb.variables$STIM, lib.size=ifnb.lib.size$STIM, BPPARAM=SnowParam(workers=2,type="SOCK"))
```

Once the 4 distributions are fitted, we next calculate the BIC value for
each model and select the model with the least BIC value.

``` r
ifnb.ctrl.bic.val <- model_bic(ifnb.ctrl.fit)
ifnb.stim.bic.val <- model_bic(ifnb.stim.fit)

#select model with least bic value
ifnb.ctrl.lbic <- lbic_model(ifnb.ctrl.bic.val, ifnb.ctrl.KS)
ifnb.stim.lbic <- lbic_model(ifnb.stim.bic.val, ifnb.stim.KS)
```

To ensure the fit of the models selected based on the least BIC value,
additionally we perform LRT to test for model adequacy and presence of
zero-inflation.

``` r
ifnb.ctrl.gof <- gof_model(ifnb.ctrl.lbic, ifnb.variables$CTRL, ifnb.lib.size$CTRL, BPPARAM=SerialParam())
ifnb.stim.gof <- gof_model(ifnb.stim.lbic, ifnb.variables$STIM, ifnb.lib.size$STIM, BPPARAM=SerialParam())
```

Finally based on the results of the model adequacy tests, we can
identify the distribution of best fit for each gene.

``` r
ifnb.ctrl.dist.fit <- select_model(ifnb.ctrl.gof)
ifnb.stim.dist.fit <- select_model(ifnb.stim.gof)
```

Once the distribution of best fit is identified for genes of interest,
it is also possible to extract parameters of interest for the models.

``` r
ifnb.ctrl.params <- model_param (ifnb.ctrl.fit, ifnb.ctrl.dist.fit, model=NULL)
ifnb.stim.params <- model_param (ifnb.stim.fit, ifnb.stim.dist.fit, model=NULL)
```

Using above results we can now identify the differentially distributed
genes between ‘CTRL’ and ‘STIM’. First we need to subset genes that is
significant in the KS test in both conditions.

``` r
#Subset the common genes between the two groups, that pass the GOF test
ifnb.ctrl.fit <- unlist(ifnb.ctrl.dist.fit)
ifnb.stim.fit <- unlist(ifnb.stim.dist.fit)
ifnb.gof.sig <- intersect(ifnb.ctrl.fit, ifnb.stim.fit)

ifnb.dist.ctrl <- data.frame(gene = c(ifnb.ctrl.dist.fit$P_genes, ifnb.ctrl.dist.fit$NB_genes, ifnb.ctrl.dist.fit$ZIP_genes, ifnb.ctrl.dist.fit$ZINB_genes))
ifnb.dist.ctrl$dist <- c(rep("Po", length(ifnb.ctrl.dist.fit$P_genes)), rep("NB", length(ifnb.ctrl.dist.fit$NB_genes)), rep("ZIP", length(ifnb.ctrl.dist.fit$ZIP_genes)), rep("ZINB", length(ifnb.ctrl.dist.fit$ZINB_genes)))

ifnb.dist.stim <- data.frame(gene = c(ifnb.stim.dist.fit$P_genes, ifnb.stim.dist.fit$NB_genes, ifnb.stim.dist.fit$ZIP_genes, ifnb.stim.dist.fit$ZINB_genes))
ifnb.dist.stim$dist <- c(rep("Po", length(ifnb.stim.dist.fit$P_genes)), rep("NB", length(ifnb.stim.dist.fit$NB_genes)), rep("ZIP", length(ifnb.stim.dist.fit$ZIP_genes)), rep("ZINB", length(ifnb.stim.dist.fit$ZINB_genes)))

#Dataframe consisting of distributions followed by each gene passing the KS test
ifnb.gof.ctrl <- ifnb.dist.ctrl[ifnb.dist.ctrl$gene %in% ifnb.gof.sig,]
ifnb.gof.stim <- ifnb.dist.stim[ifnb.dist.stim$gene %in% ifnb.gof.sig,]

ifnb.distr <- data.frame(ctrl = ifnb.gof.ctrl$dist, row.names = ifnb.gof.ctrl$gene)
ifnb.distr$stim <- ifnb.gof.stim$dist[match(rownames(ifnb.distr), ifnb.gof.stim$gene)]
```

Using the dataframe of genes and distribution followed under each
condition now we can identify genes changing distribution between ‘CTRL’
and ‘STIM’

``` r
ifnb.DD.genes <- change_shape(ifnb.distr)
```

This will give a list of two lists with genes changing distribution
between condition and genes changing distribution from unimodal in one
condition to zero-inflated in the other condition.
