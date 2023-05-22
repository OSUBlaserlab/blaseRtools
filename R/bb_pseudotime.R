#' @title Learn Graph and Calculate Pseudotime
#' @description This function determines a gene expression trajectory using `learn_graph` from monocle3 and then calculates pseudotime dimensions along this trajectory using `order_cells`.  So it is 2 functions wrapped into 1.  Usually we will not adjust the parameters for learn_graph with the possible exception of `close_loop` and `use_partition` which are also available in this function with the same defaults.  If you need to fine-tune the trajectory, use `monocle3::learn_graph` on the cds object first and then run this function to calculate pseudotime.  The graph learning will not be repeated on an object unless `force_graph` is set to `TRUE`.
#'
#' If you just want to look at the trajectory graph and not calculate pseudotime, change `calculate_pseudotime` to FALSE, or run `monocle3::learn_graph`.
#'
#' After the pseudotime values are calculated, they are handled differently than in monocle3.  In this function, they are copied from the hidden CDS slot and made an explicit cell metadata column.  Pseudotime needs a starting point or anchor.  There is no interactive option here as in monocle3.  To identify this starting point, you identify a cell metadata variable and provide it to `cluster_variable`.  This should identify a cohesive group of cells in UMAP space such as a leiden cluster, louvain cluster or partition.  Then provide a value corresponding to the cluster of interest to `cluster_value`.  The function will start pseudotime at the cell closest to the graph node in that cluster.  The pseudotime value column will be named automatically as a composite of the `cluster_variable` and `cluster_value` parameters.
#'
#' @param cds The cell data set object to calculate pseudotime upon.  Does not yet accept seurat objects.
#' @param calculate_pseudotime Logical, whether to calculate the pseudotime dimension.  If false, will only run learn_graph, Default: TRUE
#' @param cluster_variable The cell metadata column from which the pseudotime = 0 cell will be selected.
#' @param cluster_value The value of cluster_variable that identifies a cluster.  The cell closest to the root node closest to the center of this cluster will have pseudotime of 0.
#' @param use_partition Logical; If TRUE, learn_graph will construct trajectories within partitions.  If FALSE, it will connect partitions, Default: TRUE
#' @param close_loop Logical; Whether learn_graph will close looping trajectories, Default: TRUE
#' @param force_graph Logical; If TRUE, the function will recalculate the graph., Default: FALSE
#' @return A cell data set
#' @seealso
#'  \code{\link[cli]{cli_div}}, \code{\link[cli]{cli_alert}}
#'  \code{\link[monocle3]{learn_graph}}, \code{\link[monocle3]{order_cells}}, \code{\link[monocle3]{pseudotime}}
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @rdname bb_pseudotime
#' @export
#' @importFrom cli cli_div cli_alert_info
#' @importFrom monocle3 learn_graph order_cells pseudotime
#' @importFrom SummarizedExperiment colData
bb_pseudotime <-
  function(cds,
           calculate_pseudotime = TRUE,
           cluster_variable,
           cluster_value,
           use_partition = TRUE,
           close_loop = TRUE,
           force_graph = FALSE) {
    cli::cli_div(theme = list(span.emph = list(color = "orange")))
    # check and see if the principle graph needs to be calculated
    if (length(cds@principal_graph@listData) == 1) {
      cli::cli_alert_info("The graph has been learned.")
      if (force_graph) {
        cli::cli_alert_info("Re-learning graph")
        cds <-
          invisible(
            monocle3::learn_graph(
              cds = cds,
              use_partition = use_partition,
              close_loop = close_loop
            )
          )

      }
    } else {
      cli::cli_alert_info("Learning graph")
      cds <-
        invisible(
          monocle3::learn_graph(
            cds = cds,
            use_partition = use_partition,
            close_loop = close_loop
          )
        )
    }
    if (calculate_pseudotime) {
      pseudotime_colname <-
        paste0("pseudotime_", cluster_variable, "_", cluster_value)

      # identify root nodes
      cli::cli_alert_info(
        "Calculating pseudotime relative to the node nearest the center of {.emph {cluster_variable}} cluster {.emph {cluster_value}}."
      )
      cli::cli_alert_info(
        "Writing pseudotime dimensions to new cell metadata column {.emph {pseudotime_colname}}."
      )
      cds <-
        monocle3::order_cells(
          cds,
          root_pr_nodes = get_earliest_principal_node(
            cds,
            cluster_variable = cluster_variable,
            cluster_value = cluster_value
          )
        )
      SummarizedExperiment::colData(cds)[[pseudotime_colname]] <-
        monocle3::pseudotime(cds)

    }
    cds
  }

#' @importFrom SummarizedExperiment colData
#' @importFrom igraph V
#' @importFrom monocle3 principal_graph
get_earliest_principal_node <-
  function(cds, cluster_variable, cluster_value) {
    cell_ids <-
      which(SummarizedExperiment::colData(cds)[, cluster_variable] == cluster_value)

    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds),])
    root_pr_nodes <-
      igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                          (which.max(table(closest_vertex[cell_ids, ]))))]

    root_pr_nodes
  }
