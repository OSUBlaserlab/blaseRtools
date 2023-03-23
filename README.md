# blaseRtools

This R package includes commonly used functions for R analysis in the Blaser Lab.      

## Installation

You can install the latest version of blaseRtools with:

``` r
# using base function
install.packages('blaseRtools', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))

# or if using blaseRtemplates function
blaseRtemplates::install_one_package("blaseRtools", "new")

```



## Useage

All commonly used functions for the end user are prefixed "bb_".  If you load the blaseRtools package and supporting data with 

Then if you run ``library("blaseRtools")``` you can access functions by typing "bb_" and an autocomplete window with selections should pop up.

Functions related to the "Ape" and "Trace" classes are prefixed "Ape." and "Trace.", respectively.

## Function Modules

All functions are documented and have help pages which can be reviewed after installation.  Tutorials describing typical use-cases for each module of functions are linked when available.

* [scRNA-seq](pages/scRNAseq.html) 
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

* [Ape class](pages/Ape.html):  Programmatic methods for working with genbank files.
    * Ape.DNA
    * Ape.fasta
    * Ape.fimo
    * Ape.granges
    * Ape.save
    * Ape.setFeatures
    * bb_parseape
    * bb_make_ape_genomic
    * bb_make_ape_transcript
    
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
    
    
    
    
    

