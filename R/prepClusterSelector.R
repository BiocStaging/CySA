# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

#' Default Marker Pairs for 2D Plots
#'
#' Builds a default list of marker pairs from a \code{SingleCellExperiment}
#' object.
#'
#' @param sce A \code{SingleCellExperiment}.
#'
#' @return A list of character pairs.
#'
#' @keywords internal
.defaultDList <- function(sce) {
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed
    rn <- if (!is.null(colsUsed) && length(colsUsed) >= 12) colsUsed else rownames(sce)
    if (length(rn) < 12) {
        stop("sce must have at least 12 row names to build default dList")
    }
    list(
        d1 = c(rn[1], rn[2]),
        d2 = c(rn[3], rn[4]),
        d3 = c(rn[5], rn[6]),
        d4 = c(rn[7], rn[8]),
        d5 = c(rn[9], rn[10]),
        d6 = c(rn[11], rn[12])
    )
}


#' Prepare Data for the Cluster Selector Shiny App
#'
#' Subsamples a \code{SingleCellExperiment} object and builds the inputs
#' required by \code{\link{clusterSelector}}.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing the full dataset.
#' @param somFile Path to the SOM object. Used to cache/restore the subsampled
#'   version.
#' @param dList Optional list of marker pairs for 2D plots. If \code{NULL},
#'   defaults to the first 12 row names of \code{sce}.
#' @param total_cells_to_sample Number of cells to subsample in total.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with \code{sce}, \code{sce_subsampled}, and \code{dList}.
#'
#' @examples
#' sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
#' prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
#' names(prepped)
#'
#' @export
prepClusterSelectorData <- function(sce,
                                    somFile = NULL,
                                    dList = NULL,
                                    total_cells_to_sample = 100000,
                                    somCodesName = "SOM_codes",
                                    assay = "exprs",
                                    seed = 123) {
    if (length(total_cells_to_sample) != 1 || !is.numeric(total_cells_to_sample) || total_cells_to_sample < 1) {
        stop("'total_cells_to_sample' must be a positive scalar")
    }

    if (is.null(dList)) {
        dList <- .defaultDList(sce)
    }

    # Cache file handling
    cache_file <- NULL
    if (!is.null(somFile)) {
        cache_file <- paste0(tools::file_path_sans_ext(somFile), ".subsampled.RData")
    }

    if (!is.null(cache_file) && file.exists(cache_file)) {
        env <- new.env()
        load(cache_file, envir = env)
        sce_subsampled <- env$sce_subsampled
    } else {
        sce_subsampled <- withr::with_seed(seed, {
            cd <- as.data.frame(SingleCellExperiment::colData(sce))
            proportions_df <- dplyr::group_by(cd, .data$sample_id, .data$cluster_id)
            proportions_df <- dplyr::summarise(proportions_df, group_size = dplyr::n(), .groups = "drop")
            proportions_df <- dplyr::ungroup(proportions_df)
            proportions_df <- dplyr::mutate(
                proportions_df,
                total_size = sum(.data$group_size),
                proportion = .data$group_size / .data$total_size,
                n_to_sample = ceiling(.data$proportion * total_cells_to_sample + 10)
            )

            sampling_indices <- purrr::map2(
                proportions_df$n_to_sample,
                proportions_df$group_size,
                ~ {
                    if (.x > .y) {
                        return(seq_len(.y))
                    }
                    sample(.y, .x)
                }
            )
            sampling_indices <- purrr::set_names(
                sampling_indices,
                paste(proportions_df$sample_id, proportions_df$cluster_id, sep = "_")
            )
            sampling_indices <- purrr::imap(sampling_indices, ~ {
                sid <- stringr::str_replace(.y, "(^.*)_(.*)", "\\1")
                cid <- stringr::str_replace(.y, "(^.*)_(.*)", "\\2")
                sce[
                    ,
                    SingleCellExperiment::colData(sce)$sample_id == sid &
                        SingleCellExperiment::colData(sce)$cluster_id == cid
                ][, .x]
            })

            sce_subsampled <- do.call(SummarizedExperiment::cbind, sampling_indices)

            if (!is.null(cache_file)) {
                save(sce_subsampled, file = cache_file)
            }

            sce_subsampled
        })
    }

    if (length(SummarizedExperiment::assays(sce_subsampled)) == 0) {
        SummarizedExperiment::assays(sce_subsampled)[[1]] <-
            SummarizedExperiment::assays(sce)[[1]][, SummarizedExperiment::colData(sce_subsampled)$id]
        names(SummarizedExperiment::assays(sce_subsampled)) <- assay
    }

    list(sce = sce, sce_subsampled = sce_subsampled, dList = dList)
}


#' Build SOM Codes From Cell-Level Expression Data
#'
#' Computes per-cluster mean expression for all markers in \code{sce}.
#' Markers already present in the existing SOM codes metadata slot are left
#' unchanged; missing markers are calculated from the cell-level assay.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id} in
#'   \code{colData} and SOM codes stored in \code{metadata(sce)[[somCodesName]]}.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param markerList Optional character vector of markers to include. If
#'   \code{NULL}, all markers except \code{"label"} and \code{"TIME"} are used.
#' @param assay Name of the assay to use.
#'
#' @return The updated \code{SingleCellExperiment} with refreshed SOM codes.
#'
#' @keywords internal
.buildSOMCodes <- function(sce,
                           somCodesName = "SOM_codes",
                           markerList = NULL,
                           assay = "exprs") {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]

    if (is.null(markerList)) {
        markerList <- setdiff(rownames(sce), c("label", "TIME"))
    }

    # Only recompute markers that are not already present in SOM codes.
    markers_to_compute <- setdiff(markerList, colnames(somCodes))
    if (length(markers_to_compute) == 0) {
        return(sce)
    }

    exprs <- SummarizedExperiment::assay(sce, assay)[markers_to_compute, , drop = FALSE]
    cluster_id <- as.integer(SingleCellExperiment::colData(sce)$cluster_id)
    n_clusters <- nrow(somCodes)

    # Sum expression per cluster, then divide by cluster counts.
    cluster_sums <- rowsum(t(exprs), cluster_id, na.rm = TRUE)
    cluster_counts <- tabulate(cluster_id, nbins = n_clusters)
    cluster_counts[cluster_counts == 0] <- 1  # avoid division by zero; will be 0 anyway
    cluster_means <- cluster_sums / cluster_counts

    # rowsum may drop empty clusters; ensure full n_clusters rows.
    if (nrow(cluster_means) < n_clusters) {
        full <- matrix(0, nrow = n_clusters, ncol = length(markers_to_compute))
        colnames(full) <- markers_to_compute
        present <- as.integer(rownames(cluster_means))
        full[present, ] <- as.matrix(cluster_means)
        cluster_means <- full
    }

    new_codes <- cbind(somCodes, cluster_means)
    # Reorder to requested marker order, keeping existing columns first if needed.
    S4Vectors::metadata(sce)[[somCodesName]] <- new_codes[, markerList, drop = FALSE]
    sce
}


#' Build SOM Summary Statistics
#'
#' Computes per-cluster summary statistics (median, mean, third quartile, max,
#' cell count) based on distances between SOM codes and their assigned cells.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id} in
#'   \code{colData} and SOM codes in \code{metadata(sce)[[somCodesName]]}.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#'
#' @return A \code{data.frame} of per-cluster statistics.
#'
#' @keywords internal
.buildSOMStats <- function(sce,
                           somCodesName = "SOM_codes",
                           assay = "exprs") {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]
    colData <- SingleCellExperiment::colData(sce)
    exprs <- SummarizedExperiment::assay(sce, assay)
    stats <- data.frame(row.names = seq_len(nrow(somCodes)))

    for (lv in seq_len(nrow(somCodes))) {
        cells <- exprs[colnames(somCodes), colData$cluster_id == lv, drop = FALSE]
        if (ncol(cells) == 0) {
            stats[lv, c("median", "mean", "rdQu", "max", "n", "id")] <-
                c(NA, NA, NA, NA, 0, lv)
            next
        }
        d <- sqrt(colSums((t(cells) - somCodes[lv, , drop = TRUE])^2))
        s1 <- summary(d)
        stats[lv, "median"] <- s1[["Median"]]
        stats[lv, "mean"] <- s1[["Mean"]]
        stats[lv, "rdQu"] <- s1[["3rd Qu."]]
        stats[lv, "max"] <- s1[["Max."]]
        stats[lv, "n"] <- ncol(cells)
        stats[lv, "id"] <- lv
    }
    stats
}


#' Build SOM Cluster Stats
#'
#' Computes per-cluster summary statistics for each marker relative to the SOM
#' code, returning a long-format data frame suitable for downstream plotting.
#'
#' @param sce A \code{SingleCellExperiment} with \code{cluster_id} in
#'   \code{colData} and SOM codes in \code{metadata(sce)[[somCodesName]]}.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param assay Name of the assay to use.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A long \code{data.frame} with columns \code{cluster}, \code{stat},
#'   \code{variable}, and \code{value}.
#'
#' @keywords internal
.buildSOMClusterStats <- function(sce,
                                  somCodesName = "SOM_codes",
                                  assay = "exprs",
                                  verbose = TRUE) {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]
    markers <- colnames(somCodes)
    exprs <- SummarizedExperiment::assay(sce, assay)[markers, , drop = FALSE]
    cluster_id <- as.integer(SingleCellExperiment::colData(sce)$cluster_id)
    n_clusters <- nrow(somCodes)

    stats_names <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    n_stats <- length(stats_names)
    n_markers <- length(markers)

    # Precompute cluster membership as a list of cell indices once.
    cell_idx <- seq_along(cluster_id)
    idx_by_cluster <- split(cell_idx, cluster_id)

    result <- lapply(seq_len(n_clusters), function(lv) {
        if (verbose && lv %% 100 == 0) {
            cat(file = stderr(), "cluster ", lv, format(Sys.time(), "%a %b %d %X %Y"), "\n")
        }
        idx <- idx_by_cluster[[as.character(lv)]]
        if (is.null(idx) || length(idx) == 0) {
            df <- data.frame(
                cluster = rep(lv, n_markers * n_stats),
                stat = rep(stats_names, each = n_markers),
                variable = rep(markers, n_stats),
                value = 0
            )
            return(df)
        }
        cells <- exprs[, idx, drop = FALSE] - somCodes[lv, , drop = TRUE]
        # Compute summaries per marker (row). matrixStats is much faster than apply+summary.
        qs <- matrixStats::rowQuantiles(
            cells,
            probs = c(0, 0.25, 0.5, 0.75, 1),
            type = 7,
            drop = FALSE
        )
        means <- matrixStats::rowMeans2(cells)
        df <- data.frame(
            cluster = rep(lv, n_markers * n_stats),
            stat = rep(stats_names, each = n_markers),
            variable = rep(markers, n_stats),
            value = c(qs[, 1], qs[, 2], qs[, 3], means, qs[, 4], qs[, 5])
        )
        df$value[is.na(df$value)] <- 0
        df
    })
    do.call(rbind, result)
}


#' Build SOM Raster Data
#'
#' Builds a data frame representing the SOM grid, including the SOM codes and
#' per-node cell counts.
#'
#' @param sce A \code{SingleCellExperiment} with SOM codes in
#'   \code{metadata(sce)[[somCodesName]]} and optional \code{SOM_stats}.
#'
#' @return A \code{data.frame}.
#'
#' @keywords internal
.buildSOMRasterData <- function(sce, somCodesName = "SOM_codes") {
    somCodes <- S4Vectors::metadata(sce)[[somCodesName]]
    nDim <- sqrt(nrow(somCodes))
    if (nDim != as.integer(nDim)) {
        stop("SOM grid is not square: nrow(somCodes) is not a perfect square")
    }
    nDim <- as.integer(nDim)
    somRasterData <- data.frame(
        x = rep(seq_len(nDim), nDim),
        y = rep(seq_len(nDim), each = nDim)
    )
    somRasterData <- cbind(somRasterData, somCodes)

    somStats <- S4Vectors::metadata(sce)$SOM_stats
    if (!is.null(somStats) && is.data.frame(somStats) && "n" %in% colnames(somStats)) {
        n_values <- somStats$n
    } else {
        n_values <- rep(NA_integer_, nrow(somCodes))
    }
    somRasterData$N <- n_values
    colnames(somRasterData) <- make.names(colnames(somRasterData))
    somRasterData
}


#' Build Dendrogram Table
#'
#' Creates a table describing the hierarchical clustering dendrogram of the SOM
#' codes, optionally building a collapsible tree network widget.
#'
#' @param sce A \code{SingleCellExperiment} with SOM codes in
#'   \code{metadata(sce)[[somCodesName]]}.
#' @param useColTree Logical indicating whether to build the collapsible tree.
#'
#' @return A list with \code{dend}, \code{dendTable}, and \code{colTree}.
#'
#' @keywords internal
.buildDendTable <- function(sce, somCodesName = "SOM_codes", useColTree = FALSE) {
    hc <- stats::hclust(stats::dist(S4Vectors::metadata(sce)[[somCodesName]]))
    dend <- stats::as.dendrogram(hc)
    dendTable <- data.frame(
        parent = NA,
        child = "rt",
        nleaf = .nDendrogramLeaves(dend),
        indexString = ""
    )

    colTree <- NULL
    if (useColTree) {
        recDendTable <- function(dend, parent = "rt", myindexString = "") {
            if (.isDendrogramLeaf(dend)) {
                return()
            }
            for (di in seq_along(dend)) {
                indexString <- paste0(myindexString, "[[", di, "]]")
                nleaf <- .nDendrogramLeaves(dend[[di]])
                dendTable <<- rbind(
                    dendTable,
                    data.frame(
                        parent = parent,
                        child = paste0(parent, ".", di),
                        nleaf = nleaf,
                        indexString = indexString
                    )
                )
                recDendTable(dend[[di]], paste0(parent, ".", di), indexString)
            }
        }
        recDendTable(dend)
        colTree <- collapsibleTree::collapsibleTreeNetwork(dendTable, attribute = "nleaf", inputId = "node")
    }

    list(dend = dend, dendTable = dendTable, colTree = colTree)
}

#' Test if an Object is a Dendrogram Leaf
#'
#' @param x A dendrogram object.
#'
#' @return Logical.
#'
#' @keywords internal
.isDendrogramLeaf <- function(x) {
    !is.null(attr(x, "leaf")) && isTRUE(attr(x, "leaf"))
}

#' Count Leaves in a Dendrogram
#'
#' @param x A dendrogram object.
#'
#' @return Integer.
#'
#' @keywords internal
.nDendrogramLeaves <- function(x) {
    if (.isDendrogramLeaf(x)) {
        return(1L)
    }
    sum(vapply(x, .nDendrogramLeaves, integer(1)))
}


#' Prepare Cluster Selector Outputs
#'
#' Wrapper that takes a post-SOM \code{SingleCellExperiment}, builds all
#' supporting data structures required by \code{\link{clusterSelector}}, and
#' saves a \code{sce.shiny.RData} bundle compatible with the original
#' \code{prepCytoApp} output.
#'
#' @param currentRoot Character string specifying the project root directory.
#'   Output files are written to \code{currentRoot/_data/}.
#' @param sce A \code{SingleCellExperiment} with SOM clustering. If
#'   \code{NULL}, the function reads \code{currentRoot/_data/sce.RDS}.
#' @param useColTree Logical indicating whether to create a collapsible tree.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param clusterColName Column name containing cluster IDs.
#' @param total_cells_to_sample Number of cells to subsample for the Shiny app.
#' @param dList Optional list of marker pairs for 2D plots.
#' @param assay Name of the assay to use.
#' @param verbose Logical indicating whether to print progress messages.
#' @param force Logical indicating whether to force recomputation.
#' @param outputList Optional output list to include in the saved bundle.
#' @param seed Random seed for subsampling.
#'
#' @return The path \code{currentRoot}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsample the SCE using \code{prepClusterSelectorData()}, or reuse a
#'     previously subsampled SCE file if \code{subsampledSCEFile} is provided.
#'   \item Update \code{cluster_id} to a factor with all SOM clusters as levels.
#'   \item Build the cluster-by-patient table.
#'   \item Recompute missing SOM code means.
#'   \item Build SOM summary stats and cluster-level stats.
#'   \item Build SOM raster data.
#'   \item Build hierarchical clustering dendrogram and table.
#'   \item Save \code{sce.shiny.RData}.
#' }
#'
#' @param subsampledSCEFile Optional path to a previously subsampled SCE object
#'   (for example \code{sce.subsampled.RDS}). If provided and the file exists,
#'   it is reused instead of subsampling the full SCE. This avoids the memory
#'   cost of subsampling large datasets.
#'
#' @export
prepClusterSelector <- function(currentRoot,
                                sce = NULL,
                                useColTree = FALSE,
                                somCodesName = "SOM_codes",
                                clusterColName = "cluster_id",
                                total_cells_to_sample = 100000,
                                dList = NULL,
                                assay = "exprs",
                                verbose = TRUE,
                                force = FALSE,
                                outputList = NULL,
                                seed = 123,
                                subsampledSCEFile = NULL) {
    somFile <- file.path(currentRoot, "_data", "sce.RDS")
    outFile <- file.path(currentRoot, "_data", "sce.shiny.RData")

    if (!file.exists(somFile)) {
        stop("SCE file not found: ", somFile)
    }

    if (file.exists(outFile) && !force) {
        if (verbose) message("Output exists and force = FALSE; returning currentRoot")
        return(currentRoot)
    }

    if (is.null(sce)) {
        if (verbose) message("Reading SCE")
        sce <- readRDS(somFile)
    }

    if (!is.null(subsampledSCEFile) && file.exists(subsampledSCEFile)) {
        if (verbose) message("Using provided subsampled SCE file: ", subsampledSCEFile)
        if (is.null(dList)) {
            dList <- .defaultDList(sce)
        }
        sce_subsampled <- readRDS(subsampledSCEFile)
        if (!is.null(sce_subsampled) && length(SummarizedExperiment::assays(sce_subsampled)) == 0) {
            if (verbose) message("Restoring assay data into subsampled SCE")
            SummarizedExperiment::assays(sce_subsampled)[[1]] <-
                SummarizedExperiment::assays(sce)[[1]][, SummarizedExperiment::colData(sce_subsampled)$id]
            names(SummarizedExperiment::assays(sce_subsampled)) <- assay
        }
    } else {
        if (verbose) message("Subsampling SCE")
        prep <- prepClusterSelectorData(
            sce = sce,
            somFile = somFile,
            dList = dList,
            total_cells_to_sample = total_cells_to_sample,
            somCodesName = somCodesName,
            assay = assay,
            seed = seed
        )
        sce_subsampled <- prep$sce_subsampled
        dList <- prep$dList
    }

    if (verbose) message("Building cluster patient table")
    n_clusters <- nrow(S4Vectors::metadata(sce)[[somCodesName]])
    sce$cluster_id <- factor(as.character(SingleCellExperiment::colData(sce)[[clusterColName]]),
                             levels = seq_len(n_clusters))
    clusterPatientTable <- table(as.data.frame(SingleCellExperiment::colData(sce)[, c("sample_id", clusterColName)]))

    if (verbose) message("Updating SOM codes")
    sce <- .buildSOMCodes(sce, somCodesName = somCodesName, assay = assay)

    if (is.null(S4Vectors::metadata(sce)$SOM_stats)) {
        if (verbose) message("Building SOM stats")
        S4Vectors::metadata(sce)$SOM_stats <- .buildSOMStats(sce, somCodesName = somCodesName, assay = assay)
    } else {
        if (verbose) message("Reusing existing SOM_stats")
    }

    if (verbose) message("Building SOM cluster stats")
    S4Vectors::metadata(sce)[["som_clster_stats"]] <- .buildSOMClusterStats(
        sce,
        somCodesName = somCodesName,
        assay = assay,
        verbose = verbose
    )

    if (verbose) message("Building SOM raster data")
    somRasterData <- .buildSOMRasterData(sce, somCodesName = somCodesName)
    somRasterObj <- NULL

    if (verbose) message("Building dendrogram table")
    dend_res <- .buildDendTable(sce, somCodesName = somCodesName, useColTree = useColTree)
    dend_obj <- dend_res$dend
    dend_table <- dend_res$dendTable

    if (verbose) message("Saving outputs")
    save(
        sce, sce_subsampled,
        outputList, dList,
        dend = dend_obj,
        dendTable = dend_table,
        clusterPatientTable,
        somRasterData,
        somRasterObj,
        file = outFile,
        compress = FALSE
    )

    if (verbose) message("Saved: ", outFile)
    currentRoot
}
