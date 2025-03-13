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

# blaseRtools 0.0.0.9105

* added label option to bb_var_umap

# blaseRtools 0.0.0.9106

* added option to show alt experiment in bb_rowmeta
* debugged bb_cite_umap
* changed x, y labels in bb_gene_umap

# blaseRtools 0.0.0.9107-108

* added option to plot antibody capture in gene bubbles
* made changes in bb_aggregate to support this

# blaseRtools 0.0.0.9109

* added expression_threshold to genebubbles

# blaseRtools 0.0.0.9110

* debugged cite umap

# blaseRtools 0.0.0.9111

* debugged gene umap plotting aggregate genes

# blaseRtools 0.0.0.9112-3

* debugged bb_rowmeta and bb_aggregate experiment names parameter

# blaseRtools 0.0.0.9114

* added normalize_batch() 

# blaseRtools 0.0.0.9115

* edited bb_citeseq to optionally order cells to plot

# blaseRtools 0.0.0.9116

* removed matrix.utils dependency

# blaseRtools 0.0.0.9117-24

* removed seurat wrappers dependency
* removed renv dependency
* edited git ignore
* added vignettes to buildignore
* removed examples

# blaseRtools 0.0.0.9117-25

* added assay option to cite umap

# blaseRtools 0.0.0.9126 - 30

* major update to ape functions
* refactored
* corrected bugs
* added transcriptome APE function

# blaseRtools 0.0.0.9132

* added bb_cluster_representation2

# blaseRtools 0.0.0.9133-4

* added bb_load_tensx_h5
* edited bb_split_citeseq to recalculate size factors after removing the antibody matrix

# blaseRtools 0.0.0.9135-6

* debugged seurat_anno

# blaseRtools 0.0.0.9137-8

* pseuodime functions

# blaseRtools 0.0.0.9139

* added plot geens in pseudotime

# blaseRtools 0.0.0.9140-1

* fixed transcriptome ape function

# blaseRtools 0.0.0.9142

* added outline function to bb_var_umap

# blaseRtools 0.0.0.9143

* fixed bug where factors got lost from bb_genebubbles

# blaseRtools 0.0.0.9144

* updated blind and unblind functions to use fs::path
* added blinding vignette

# blaseRtools 0.0.0.9145-6

* major updates to trace functions

# blaseRtools 0.0.0.9147

* updated ape function to put placeholder feature in gbk files without features.

# blaseRtools 0.0.0.9148

* updated bb_doubletfinder to work with seurat v5

# blaseRtools 0.0.0.9149-51

* fixed bugs in trace functions related to trimming levels for the peaks slot.
* added option to debug the gene model track by showing the plotted transcript IDs

# blaseRtools 0.0.0.9152-3

* fixed bug where density var umap plots with two dimensions weren't being calculated correctly.


# blaseRtools 0.0.0.9154 -55

* added bb_fragment_replacement

# blaseRtools 0.0.0.9156

* edited get_seurat_clr to fix bug introduced by Seurat v5

# blaseRtools 0.0.0.9157

* fixed bug in bb_cds_anno where joining the new cell annotations back onto the cds cell metadata failed due to failure to specify the joining columns.  This was only evident with cds queries containing conflicting cell metadata columns like RNA counts etc and resulted in NA assignments.

# blaseRtools 0.0.0.9158-9

* fixed bug in bb_cite_umap where we were matching additional antibodies inappropriately

# blaseRtools 0.0.0.9160:  

* handled case when certain single cell objects do not have a gene_short_name column
* applies to gene_umap and genebubbles
* internal function get_gene_ids

# blaseRtools 0.0.0.9161 

* fixed a bug where filter_cds returned all cells instead of 0 cells

# blaseRtools 0.0.0.9162

* parameterized bb_qc to change nmads

# blaseRtools 0.0.0.9163

* added minimum segment length parameter to bb_var_umap

# blaseRtools 0.0.0.9164

* added bb_loupeR function

# blaseRtools 0.0.0.9165

* added bb_daseq

# blaseRtools 0.0.0.9166

* added hexify to bb_var_umap
* changed the interaction label used by bb_cellchat_heatmap because of duplications in the cellchat database.  This may cause cosmetic changes in old code.

# blaseRtools 0.0.0.9167

* fixed bug

# blaseRtools 0.0.0.9168

* fixed bug in bb_seurat_anno where cds was misidentified

# blaseRtools 0.0.0.9169

* added monocle anno

# blaseRtools 0.0.0.9170

* imports newest monocle3

# blaseRtools 0.0.0.9171

* fixes bug in cellchat to work with on disk matrices

# blaseRtools 0.0.0.9173

* added ImageClass

# blaseRtools 0.0.0.9174-6

* renamed as_tibble to ImageCatalog.as_tibble
* fixed scRNA-seq vignette

# blaseRtools 0.0.0.9177-81

* pseudocells
* monocle projection

# blaseRtools 0.0.0.9182

* SummarizedHeatmap

# blaseRtools 0.0.0.9183

* Updated trace functions to fill in gaps in sparse traces

# blaseRtools 0.0.0.9184

* summarized heatmap

# todo

* biaxial cite-seq plots
* consensus reference columns
