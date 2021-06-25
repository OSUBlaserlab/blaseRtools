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

# blaseRtools 0.0.0.9009

* changed bb_load_multi_counts to be compatible with cellranger v6
