# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# utils_clusterSelector.R ----
# Internal app builder. Validates inputs, pre-computes cached objects, and
# returns the UI/server pair for clusterSelector().

#' Build clusterSelector UI and Server
# '
# ' Creates the UI/server pair for the clusterSelector Shiny app.
#'
#' @inheritParams clusterSelector
#'
#' @return A named list with elements \code{ui} and \code{server}.
#'
#' @keywords internal
#' @name INTERNAL_buildClusterSelectorApp
.buildClusterSelectorApp <- function(sce,
                                     sce_subsampled,
                                     outputList = list(),
                                     colTree = NULL,
                                     dList,
                                     dend,
                                     dendTable,
                                     clusterPatientTable,
                                     somCodesName = "SOM_codes",
                                     nPlots = 6,
                                     somRasterData,
                                     somRasterObj,
                                     fsom = NULL,
                                     env = environment(),
                                     verbose = FALSE) {
    # nocov start
    # for (idx in seq_along(dList)) {
    #     assign(paste0("d", idx, ".1"), dList[[idx]][1])
    #     assign(paste0("d", idx, ".2"), dList[[idx]][2])
    # }

    metaD <- .validateClusterSelectorInputs(
        sce, sce_subsampled, outputList, dList, dend, dendTable,
        clusterPatientTable, somRasterData, somCodesName, nPlots, fsom
    )

    rnSCE <- rownames(sce)

    jscode <- "shinyjs.closeWindow = function() { window.close(); }"

    if ("outputList" %in% ls(envir = env)) {
        outputList <- get("outputList", envir = env)
    }
    cluster_levels <- if (is.factor(sce$cluster_id)) levels(sce$cluster_id) else sort(unique(as.integer(sce$cluster_id)))
    outputList <- .initializeOutputList(outputList, cluster_levels)
    colsUsed <- metaD$map$colsUsed
    keepCols <- intersect(c("x", "y", "id", colsUsed), colnames(somRasterData))
    somRasterData <- somRasterData[, keepCols]
    if (!is.null(somRasterObj)) {
        if (inherits(somRasterObj, c("RasterBrick", "RasterStack"))) {
            somRasterObj <- raster::subset(somRasterObj, colsUsed)
        } else if (is.array(somRasterObj)) {
            # Subset array layers by name
            layer_idx <- which(dimnames(somRasterObj)[[3]] %in% colsUsed)
            somRasterObj <- somRasterObj[, , layer_idx, drop = FALSE]
        }
    }

    # Pre-build the SOM raster heatmap as a ggplot object once.
    # The overlay (selected nodes) is added per-render as ggplot layers.
    # Using ggplot avoids the lattice/base-graphics incompatibility that
    # broke the previous recordPlot/replayPlot approach.
    baseRasterGgplot <- .buildBaseRasterGgplot(somRasterData, colsUsed)

    # Pre-compute per-channel axis limits once from the full assay matrix.
    # This avoids scanning assays(sce)[[1]] every time dimSelection changes.
    assay_mat <- SummarizedExperiment::assays(sce)[[1]]
    channelLimits <- lapply(rownames(sce), function(ch) {
        vals <- assay_mat[ch, ]
        c(min = min(vals), max = max(vals))
    })
    names(channelLimits) <- rownames(sce)

    assign(x = "outputList", value = outputList, envir = env)


    # nPlots = 15
    # nPlots = 6


    ui <- shinydashboard::dashboardPage(
        shinydashboard::dashboardHeader(title = "Cluster Selector"),
        sidebar = .buildClusterSelectorSidebar(sce, outputList, nPlots),
        body = shinydashboard::dashboardBody(
            .buildFirstBodyRow(colTree, nPlots),
            .buildSOM2DPlotsBox(nPlots, colsUsed, outputList,
                                markers = rnSCE, dList = dList),
            .buildSixStaticSOMBox(markers = rnSCE, dList = dList),
            .buildFlowSOMPieBox(),
            .buildFlowSOMStarsBox(fsom),
            .buildStatsBox(metaD, somCodesName, outputList),
            .buildSOMRasterBox(),
            .buildViolinBox(metaD[[somCodesName]], colsUsed),
            .buildUpSetBox(outputList)
        )
    )

    server <- .buildClusterSelectorServer(
        sce = sce, sce_subsampled = sce_subsampled,
        outputList = outputList, colTree = colTree, dList = dList,
        dend = dend, dendTable = dendTable,
        clusterPatientTable = clusterPatientTable,
        somCodesName = somCodesName, nPlots = nPlots,
        somRasterData = somRasterData, somRasterObj = somRasterObj,
        fsom = fsom, env = env, verbose = verbose
    )

    list(ui = ui, server = server)
    # nocov end
}

