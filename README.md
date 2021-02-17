
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
library(Seurat)
library(SeuratData)
#> Registered S3 method overwritten by 'cli':
#>   method     from    
#>   print.boxx spatstat
#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used
#> -- Installed datasets ------------------------------------- SeuratData v0.2.1 --
#> v ifnb 3.1.0
#> -------------------------------------- Key -------------------------------------
#> v Dataset loaded successfully
#> > Dataset built with a newer version of Seurat than installed
#> (?) Unknown version of Seurat installed
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
