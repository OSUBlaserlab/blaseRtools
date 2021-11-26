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

# blaseRtools 0.0.0.9034-5

* added hg38_ape

