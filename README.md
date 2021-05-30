# blaseRtools

This R package includes commonly used functions for R analysis in the Blaser Lab.  Most are related to single cell analysis and are derived from monocle with some modifications.  Helper functions for general data analysis are also present.    

## Installation

You can install the released version of blaseRtools from with:

``` r
devtools::install_github("git@github.com:blaserlab/blaseRtools.git")
```

## Example

Create a biaxial UMAP plot of single cell data highlighting cell partitions:

``` r
library(blaseRtools)
bb_var_umap(cds = cds, var = "partition")
```
