# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# validate_inputs.R ----
# Input validation helpers for clusterSelector(). Moving these checks out of
# the main function keeps the factory focused on wiring UI and server logic.

#' Validate clusterSelector Inputs
#'
#' Checks that all required inputs are present and correctly typed before the
#' Shiny application is constructed.
#'
#' @param sce Full \code{SingleCellExperiment}.
#' @param sce_subsampled Subsampled \code{SingleCellExperiment}.
#' @param outputList Named list of cluster groupings.
#' @param dList List of marker pairs.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param nPlots Number of 2D SOM plots.
#' @param fsom Optional \code{FlowSOM} object.
#'
#' @return The validated \code{metadata(sce)} object invisibly.
#'
#' @keywords internal
.validateClusterSelectorInputs <- function(sce, sce_subsampled, outputList, dList,
                                           dend, dendTable, clusterPatientTable,
                                           somRasterData, somCodesName, nPlots, fsom) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment")
    }
    if (!inherits(sce_subsampled, "SingleCellExperiment")) {
        stop("'sce_subsampled' must be a SingleCellExperiment")
    }

    metaD <- S4Vectors::metadata(sce)
    if (!somCodesName %in% names(metaD)) {
        stop("'", somCodesName, "' not found in metadata(sce)")
    }
    if (is.null(metaD$SOM_stats)) {
        stop("'SOM_stats' not found in metadata(sce)")
    }
    if (is.null(metaD$map) || is.null(metaD$map$colsUsed)) {
        stop("metadata(sce)$map$colsUsed is required")
    }

    if (!"cluster_id" %in% names(SingleCellExperiment::colData(sce))) {
        stop("'cluster_id' not found in colData(sce)")
    }
    if (!"sample_id" %in% names(SingleCellExperiment::colData(sce))) {
        stop("'sample_id' not found in colData(sce)")
    }

    if (!is.list(dList) || length(dList) < nPlots) {
        stop("'dList' must be a list with at least ", nPlots, " marker pairs")
    }

    if (!inherits(dend, "dendrogram")) {
        stop("'dend' must be a dendrogram object")
    }
    if (!is.data.frame(dendTable)) {
        stop("'dendTable' must be a data frame")
    }

    if (!is.table(clusterPatientTable)) {
        stop("clusterPatientTable not a table")
    }
    if (!"cluster_id" %in% names(dimnames(clusterPatientTable))) {
        stop("cluster_id not in colnames(clusterPatientTable)")
    }

    if (!is.data.frame(somRasterData)) {
        stop("'somRasterData' must be a data frame")
    }
    if (!all(c("x", "y") %in% colnames(somRasterData))) {
        stop("'somRasterData' must contain 'x' and 'y' columns")
    }
    missing_cols <- setdiff(metaD$map$colsUsed, colnames(somRasterData))
    if (length(missing_cols) > 0) {
        stop("somRasterData is missing SOM columns: ", paste(missing_cols, collapse = ", "))
    }

    if (!is.null(fsom)) {
        if (!inherits(fsom, "FlowSOM")) {
            stop("'fsom' must be a FlowSOM object")
        }
        if (is.null(fsom$MST) || is.null(fsom$MST$l)) {
            stop("'fsom$MST$l' is required for the interactive star plot")
        }
    }

    invisible(metaD)
}


#' Initialize Output List for the Cluster Selector
#'
#' Ensures that the \code{outputList} always contains a \code{"Rest"} group
#' covering all cluster levels not already assigned to a named group.
#'
#' @param outputList Named list of cluster groupings.
#' @param clusterLevels Character or integer vector of all cluster levels.
#'
#' @return Updated named list with \code{"Rest"} populated and empty groups
#'   removed.
#'
#' @keywords internal
.initializeOutputList <- function(outputList, clusterLevels) {
    outputList[["Rest"]] <- c()
    used <- unique(unlist(outputList))
    all_levels <- clusterLevels
    outputList[["Rest"]] <- as.integer(all_levels[!all_levels %in% used])

    for (na in names(outputList)) {
        if (length(outputList[[na]]) == 0) outputList[[na]] <- NULL
    }
    outputList
}
