# blaseRtools

This R package includes commonly used functions for R analysis in the Blaser Lab.      

## Installation

You can install the latest version of blaseRtools with:

``` r
devtools::install_github("blaserlab/blaseRtools", build_vignettes = TRUE)
```

The package uses a number of precomputed data objects which we keep in a separate repository.  It should be installed automatically as a dependency but if not it can be installed using:

``` r
devtools::install_github("blaserlab/blaseRdata")
```

All commonly used functions for the end user are prefixed "bb_".  If you load the blaseRtools package and supporting data with 

``` r
library(blaseRtools)
library(blaseRdata)
```

Then you can access functions by typing "bb_" and an autocomplete window with selections should pop up.

Functions related to the "Ape" and "Trace" classes are prefixed "Ape." and "Trace.", respectively.

## Function Modules

All functions are documented and have help pages which can be reviewed after installation.  Tutorials describing typical use-cases for each module of functions are linked when available.

* [scRNA-seq]({{ site.baseurl }}{% link pages/scRNAseq.html %}) 
    * bb_align
    * bb_annotate_npc
    * bb_cellmeta
    * bb_cluster_representation
    * bb_doubletfinder
    * bb_gene_dotplot
    * bb_gene_modules
    * bb_gene_pseudotime
    * bb_gene_umap
    * bb_gene_violinplot
    * bb_goenrichment
    * bb_gosummary
    * bb_load_tenx_targz
    * bb_monocle_regression
    * bb_pseudobulk_mf
    * bb_qc
    * bb_rejoin
    * bb_rowmeta
    * bb_seurat_anno
    * bb_triplecluster
    * bb_var_umap

* [Ape class]({{ site.baseurl }}{% link pages/Ape.html %}):  Programmatic methods for working with genbank files.
    * Ape.DNA
    * Ape.fasta
    * Ape.fimo
    * Ape.granges
    * Ape.save
    * Ape.setFeatures
    * bb_grcz11_ape
    * bb_hg38_ape
    * bb_parseape
    
* Trace class:  A container for working with range-based data from ATAC and ChIP-seq experiments.
    * bb_buff_granges
    * bb_makeTrace
    * bb_merge_narrowpeaks
    * bb_metafeature
    * bb_plot_trace_axis
    * bb_plot_trace_data
    * bb_plot_trace_feature
    * bb_plot_trace_links
    * bb_plot_trace_model
    * bb_plot_footprint
    * bb_promoter_overlap
    * bb_read_bam
    * bb_read_narrowpeak
    * Trace.data
    * Trace.features
    * Trace.gene_model
    * Trace.links
    * Trace.plot_range
    * Trace.setFeatures
    * Trace.setLinks
    * Trace.setRange
    
* Image blinding for quantitative analysis
    * bb_blind_images
    * bb_unblind_images
    
    
    
    
    

