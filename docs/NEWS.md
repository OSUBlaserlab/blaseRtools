# blaseRtools 0.0.0.9000

* Made a new gene dotplot function to allow multifactorial plotting.  Added a function to color cells by local number in bb_var_umap.

# blaseRtools 0.0.0.9001

* Made a new function, bb_align, calculates aligned umap coordinates.  These are inserted into the internal dimension slot of the CDS object.  Prealignment coordinates are added as new cell metadata columns.

# blaseRtools 0.0.0.9002

* Fixed a bug in bb_var_umap for local_n and log_local_n plotting when faceting by 1 dimension.

# blaseRtools 0.0.0.9003

* Added the option to downsample faceted plots to the bb_var_umap function. 

# blaseRtools 0.0.0.9004

* Fixed go term functions. 

# blaseRtools 0.0.0.9005

* Made it possible to explicitly order axes in bb_gene_dotplot 

# blaseRtools 0.0.0.9006-7

* Edits to violin plot to allow automatic matching of point to violin.

# blaseRtools 0.0.0.9008

* reworked print stats report

# blaseRtools 0.0.0.9009-10

* changed bb_load_multi_counts to be compatible with cellranger v6, changed qc to allow for pdx

# blaseRtools 0.0.0.9011

* added bb_renv_datapkg:  a function to load binary data packages from local sources

# blaseRtools 0.0.0.9012

* added bb_goenrichment back to namespace; edited bb_renv_datapkg to allow installation of explicit binary packages

# blaseRtools 0.0.0.9013 - 15

* updated bb_renv_datapkg to check for and install latest package version. 

# blaseRtools 0.0.0.9016-17

* new function bb_load_cloud_counts for compatibility with 10X cloud preprocessing output. 

# blaseRtools 0.0.0.9018

* added cell and feature metadata functions. 

# blaseRtools 0.0.0.9019-20

* added blind and unblind functions 

# blaseRtools 0.0.0.9021-25

* edited blind function and tbl_to_coldata and tbl_to_rowdata functions 

# blaseRtools 0.0.0.9026

* add bb_plotfootprint function from signac

# blaseRtools 0.0.0.9027-32

* added ape class and functions

# blaseRtools 0.0.0.9033

* fixed glitch in tbl to coldata and tbl to rowdata

# blaseRtools 0.0.0.9034-6

* added hg38_ape

# blaseRtools 0.0.0.9037

* fixed plyranges namespace issue

# blaseRtools 0.0.0.9038

* fixed Ape.fimo color bug

# blaseRtools 0.0.0.9039

* Added trace class and functions

# blaseRtools 0.0.0.9040

* Removed dependency on plyranges since it was causing so many conflicts and was difficult to install.
* Changed for importing functions from DESeq2 and DoubletFinder

# blaseRtools 0.0.0.9041

* Fixed load cloud counts not to err on multi-genome samples.  Now it will load all genomes and then you have to remove the ones you don't want from the cds.  Just easier to do it that way. 

# blaseRtools 0.0.0.9042

* externalized data objects for certain functions

# blaseRtools 0.0.0.9045

* added metafeature functions

# blaseRtools 0.0.0.9046

* added bb_load_tenx_targz and removed other 10X loading functions since they were very confusing.

# blaseRtools 0.0.0.9047

* added bb_cluster_representation

# blaseRtools 0.0.0.9048

* added bb_pseudobulk_mf

# blaseRtools 0.0.0.9049

* added scRNA-seq vignette

# blaseRtools 0.0.0.9050

* added bb_extract_msig

# blaseRtools 0.0.0.9051

* fixed pseudobulk mf

# blaseRtools 0.0.0.9052

* Ape vignette, fixed bugs in Ape functions and bb_renv_datapkg.

# blaseRtools 0.0.0.9053-4

* Changes to bb_gene_dotplot

# blaseRtools 0.0.0.9055

* Changed dependencies

# blaseRtools 0.0.0.9056

* Added filter_cds

# blaseRtools 0.0.0.9057

* edited dependencies

# blaseRtools 0.0.0.9058

* added bb_cds_anno and bb_cds_heatmap functions

# blaseRtools 0.0.0.9059 - 61

* edited DESCRIPTION to include ComplexHeatmap 

# blaseRtools 0.0.0.9064

* edited dependencies to include tidyverse as depends 

# blaseRtools 0.0.0.9065

* added citeseq functions

# blaseRtools 0.0.0.9066

* increased default brightness of bb_gene_umap color scale

# blaseRtools 0.0.0.9067

* fixed bug in citeseq functions where object got converted to wrong class

# blaseRtools 0.0.0.9068

* fixed spot in pseudobulk where it converted sparse to dense unnecessarily

# blaseRtools 0.0.0.9069

* added cellchat functions

# blaseRtools 0.0.0.9075-9

* edited pseudobulk mf to choose vst with large number of samples
* added back cellchat
* silenced warnings on bb_gene_umap and bb_var_umap

# blaseRtools 0.0.0.9080

* edited dependencies

# blaseRtools 0.0.0.9090 - 94

* ported many essential functions to take in Seurat and cds objects.

# blaseRtools 0.0.0.9095

* added max expression value parameter to bb_gene_umap

# blaseRtools 0.0.0.9096

* added geom split violin

# blaseRtools 0.0.0.9097

* reversed order of link plotting in trace_funcs.R
* now highest scoring link is plotted on top

# blaseRtools 0.0.0.9098-100

* added rasterize option to umap and violin plot functions 

# blaseRtools 0.0.0.9101

* added alt dims to bb_gene_umap

# blaseRtools 0.0.0.9102

* added font face option to plot trace model function

# blaseRtools 0.0.0.9103

* added split atac function
* new parameter in bb_qc to set alternate cutoffs

# blaseRtools 0.0.0.9104

* edited split citeseq withDimnames param
