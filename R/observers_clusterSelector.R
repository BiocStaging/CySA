# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# observers_clusterSelector.R ----
# Shiny observers for the clusterSelector app. Separating observers from
# outputs makes the server wiring easier to read and aligns with Bioconductor's
# Shiny app packaging guidelines.

#' Register clusterSelector Observers
#'
#' Attaches all `shiny::observe*` calls used by the clusterSelector app.
#' Output renderers live in `.registerClusterSelectorOutputs()`.
#'
#' @param input,output,session Shiny objects.
#' @param rv Reactive values object holding `outputList`.
#' @param rsUsed,rsUsed_d Reactive values for the current node selection.
#' @param activePlot Reactive value tracking the active SOM plot.
#' @param dListRV Reactive value holding the list of marker pairs.
#' @param dimSelection Reactive value holding per-plot axis limits.
#' @param triggerRedraw Reactive value used to force plot redraws.
#' @param selectedUpdate2,selectedUpdate Reactive counters for invalidation.
#' @param updatedoutputList Function that updates UI choices from outputList.
#' @param sample2PlotDb Debounced reactive of samples to plot.
#' @param dfPlot Reactive data frame for the scatter plot.
#' @param groupsInput Debounced reactive of sample-group selections.
#' @param inputClusterNumber Debounced reactive of typed cluster ids.
#' @param zoomFunc Function handling plotly zoom relayout events.
#' @param sce Full `SingleCellExperiment`.
#' @param sce_subsampled Subsampled `SingleCellExperiment`.
#' @param metaD `metadata(sce)` list.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param colsUsed Character vector of SOM columns used.
#' @param nPlots Number of 2D SOM plots.
#' @param verbose Logical indicating whether to print debug messages.
#'
#' @keywords internal
#' @name INTERNAL_registerClusterSelectorObservers
.registerClusterSelectorObservers <- function(
        input, output, session,
        rv, rsUsed, rsUsed_d, activePlot, dListRV, dimSelection,
        triggerRedraw, selectedUpdate2, selectedUpdate, updatedoutputList,
        sample2PlotDb, dfPlot, groupsInput, inputClusterNumber,
        zoomFunc, sce, sce_subsampled, metaD, dend, dendTable,
        clusterPatientTable, somCodesName, colsUsed, nPlots, verbose,
        fsom = NULL, env
) {
    # Keep external env in sync for callers that still read env$outputList.
    shiny::observe({
        assign(x = "outputList", value = rv$outputList, envir = env)
    })

    # Force redraw when the implicit "selected" group is used for colouring.
    shiny::observe({
        rs <- rsUsed_d()
        if ("selected" %in% input$colorbyGroups) {
            shiny::isolate(selectedUpdate2(selectedUpdate2() + 1))
        }
    })

    # Update active plot when the user picks a preset pair.
    shiny::observeEvent(input$dimPairSelect, {
        val <- input$dimPairSelect
        if (!is.null(val) && val %in% seq_len(nPlots)) {
            activePlot(as.integer(val))
        }
    })

    # Observe typed cluster ids and update the current selection.
    shiny::observe({
        ic <- inputClusterNumber()
        ic <- suppressWarnings(as.integer(ic))
        ic <- ic[!is.na(ic)]
        valid <- as.integer(intersect(as.character(ic),
                                      colnames(clusterPatientTable)))
        if (length(valid) > 0L) {
            shiny::isolate(rsUsed(valid))
        }
    })

    # Sync clusterNumbers text with the current selection.
    shiny::observe({
        rs <- as.integer(unlist(rsUsed()))
        rs <- sort(rs[!is.na(rs)])
        outputList <- rv$outputList
        updateOL <- FALSE
        if (!"selected" %in% names(outputList)) updateOL <- TRUE
        outputList$selected <- rs
        rv$outputList <- outputList
        if (updateOL) updatedoutputList()
        shiny::updateTextInput(
            inputId = "clusterNumbers",
            value = paste(rs, collapse = ", ")
        )
    })

    # Apply a named group from the current selection.
    shiny::observeEvent(input$applyName, {
        outputList <- rv$outputList
        cName <- input$clusterName
        if (is.null(cName) || !nzchar(cName)) return(NULL)
        cList <- inputClusterNumber()
        outputList[[cName]] <- cList
        outputList <- .rebuildOutputList(outputList, levels(sce$cluster_id))
        rv$outputList <- outputList
        updatedoutputList()
        currentSelection <- input$colorbyGroups
        shiny::updateSelectInput(
            inputId = "colorbyGroups",
            choices = names(outputList),
            selected = c(currentSelection, cName)
        )
    })

    # Remove a single named group.
    shiny::observeEvent(input$rmGrp, {
        cName <- input$clusterNameRM
        if (is.null(cName) || !nzchar(cName) ||
            !cName %in% names(rv$outputList)) {
            return(NULL)
        }
        outputList <- rv$outputList
        outputList[[cName]] <- NULL
        outputList <- .rebuildOutputList(outputList, levels(sce$cluster_id))
        rv$outputList <- outputList
        updatedoutputList()
    })

    # Remove multiple named groups.
    shiny::observeEvent(input$rmGroups, {
        rmGroups <- input$groupRM
        if (is.null(rmGroups) || length(rmGroups) == 0L) return(NULL)
        outputList <- rv$outputList
        outputList[rmGroups] <- NULL
        rs <- shiny::isolate(rsUsed())
        rmCluster <- unique(unlist(rv$outputList[rmGroups]))
        if (!is.null(rs)) {
            rsUsed(setdiff(rs, rmCluster))
        }
        outputList <- .rebuildOutputList(outputList, levels(sce$cluster_id))
        rv$outputList <- outputList
        updatedoutputList()
    })

    # Activate a saved group in the clusterNumbers input.
    shiny::observeEvent(input$clusterNameSelect, {
        outputList <- rv$outputList
        listNames <- names(outputList) %in% input$clusterNameSelect
        combinedSoms <- unique(unlist(outputList[listNames]))
        if (length(combinedSoms) > 0L) {
            shiny::updateTextInput(
                session = session,
                inputId = "clusterNumbers",
                value = paste(combinedSoms, collapse = ", ")
            )
        }
    })

    # Update group1/group2 choices when the grouping variable changes.
    shiny::observeEvent(input$groupsVar, {
        groupsVar <- input$groupsVar
        if (is.null(groupsVar) ||
            !groupsVar %in% colnames(metaD$experiment_info)) {
            return(NULL)
        }
        levs <- levels(metaD$experiment_info[, groupsVar])
        shiny::updateSelectInput(session = session, inputId = "group1",
                                 choices = levs)
        shiny::updateSelectInput(session = session, inputId = "group2",
                                 choices = levs)
    })

    # Keep group1 and group2 mutually exclusive.
    shiny::observe({
        groupsVar <- input$groupsVar
        if (is.null(groupsVar) ||
            !groupsVar %in% colnames(metaD$experiment_info)) {
            return(NULL)
        }
        levs <- levels(metaD$experiment_info[, groupsVar])
        grp1 <- input$group1
        grp2 <- shiny::isolate(input$group2)
        levs <- levs[!levs %in% grp1]
        shiny::updateSelectInput(
            session = session, inputId = "group2",
            choices = levs, selected = grp2
        )
    })

    shiny::observe({
        groupsVar <- input$groupsVar
        if (is.null(groupsVar) ||
            !groupsVar %in% colnames(metaD$experiment_info)) {
            return(NULL)
        }
        levs <- levels(metaD$experiment_info[, groupsVar])
        grp2 <- input$group2
        grp1 <- shiny::isolate(input$group1)
        levs <- levs[!levs %in% grp2]
        shiny::updateSelectInput(
            session = session, inputId = "group1",
            choices = levs, selected = grp1
        )
    })

    # Update subsampled data frame when sample filter changes.
    shiny::observe({
        sampleIds <- sample2PlotDb()
        if (is.null(sampleIds)) return(NULL)
        selected <- sce_subsampled[, sce_subsampled$sample_id %in% sampleIds]
        dfPlot(data.frame(t(SummarizedExperiment::assays(selected)[[1]])))
    })

    # Plotly selection observers for SOM plots.
    lapply(seq_len(nPlots), function(i) {
        shiny::observeEvent(
            suppressWarnings(
                .safeEventData(verbose = verbose, "plotly_selected",
                               source = paste0("somData", i))
            ),
            {
                if (verbose) message("som", i, " touched")
                activePlot(i)
                rs <- shiny::isolate(rsUsed_d())
                shiny::req(rs)
                d <- .safeEventData(
                    verbose = verbose,
                    "plotly_selected",
                    source = paste0("somData", i)
                )
                if (is.null(d)) return(NULL)
                d <- .inputSelect(d, rs, shiny::isolate(input$selectMode))
                shiny::isolate(rsUsed(d))
            }
        )
    })

    # Zoom observers for SOM plots.
    lapply(seq_len(nPlots), function(plotIdx) {
        shiny::observe({
            zoom <- suppressWarnings(
                .safeEventData(verbose = verbose, "plotly_relayout",
                               source = paste0("somData", plotIdx))
            )
            zoomFunc(zoom, plotIdx)
        })
    })

    # somDataMain selection observer.
    shiny::observeEvent(
        suppressWarnings(
            .safeEventData(verbose = verbose, "plotly_selected",
                           source = "somDataMain")
        ),
        {
            d <- .safeEventData(verbose = verbose, "plotly_selected",
                                source = "somDataMain")
            rs <- shiny::isolate(rsUsed_d())

            if (is.null(d)) return(NULL)
            rs <- shiny::isolate(rsUsed_d())
            shiny::req(rs)
            d <- .inputSelect(d, rs, shiny::isolate(input$selectMode))
            shiny::isolate(rsUsed(d))
        }
    )

    # somDataMain zoom observer.
    shiny::observe({
        zoom <- suppressWarnings(
            .safeEventData(verbose = verbose, "plotly_relayout",
                           source = "somDataMain")
        )
        zoomFunc(zoom, 1L)
    })

    # Dimension-reduction plot selection observers.
    for (src in c("tsne", "umap", "pca")) {
        local({
            src_local <- src
            shiny::observeEvent(
                suppressWarnings(
                    .safeEventData(verbose = verbose, "plotly_selected",
                                   source = src_local)
                ),
                {
                    if (verbose) message(src_local, " touched")
                    rs <- shiny::isolate(rsUsed_d())
                    if (is.null(rs)) return(NULL)
                    d <- .safeEventData(
                        verbose = verbose,
                        "plotly_selected",
                        source  = src_local   # ← use local copy
                    )
                    if (is.null(d)) return(NULL)
                    d <- .inputSelect(d, rs, shiny::isolate(input$selectMode))
                    shiny::isolate(rsUsed(d))
                }
            )
        })
    }

    # Dendrogram selection observers.
    shiny::observeEvent(
        suppressWarnings(
            .safeEventData(verbose = verbose, "plotly_selected",
                           source = "dendPlotly")
        ),
        {
            d <- .safeEventData(
                verbose = verbose,
                "plotly_selected",
                source = "dendPlotly"
            )
            if (is.null(d) || nrow(d) == 0L || is.null(d$customdata)) {
                return(NULL)
            }
            clicked <- as.integer(d$customdata)
            rs <- shiny::isolate(rsUsed_d())
            mode <- shiny::isolate(input$selectMode)
            new_rs <- .applySelectMode(rs, clicked, mode)
            shiny::isolate(rsUsed(new_rs))
        }
    )

    shiny::observeEvent(
        suppressWarnings(
            .safeEventData(verbose = verbose, "plotly_click",
                           source = "dendPlotly")
        ),
        {
            d <- .safeEventData(
                verbose = verbose,
                "plotly_click",
                source = "dendPlotly"
            )
            # message("DIAG d = ", if (is.null(d)) "NULL" else paste(capture.output(str(d)), collapse = " | "))
            # message("DIAG nrow(d) = ", paste(capture.output(print(nrow(d))), collapse = " | "))
            if (is.null(d) || nrow(d) == 0L || is.null(d$customdata)) {
                # message("DIAG: guard triggered, returning early")
                return(NULL)
            }
            clicked <- as.integer(d$customdata)
            rs <- shiny::isolate(rsUsed_d())
            mode <- shiny::isolate(input$selectMode)
            new_rs <- .applySelectMode(rs, clicked, mode)
            shiny::isolate(rsUsed(new_rs))
        }
    )

    # SOM raster grid selection observer.
    shiny::observeEvent(
        suppressWarnings(
            .safeEventData(verbose = verbose, "plotly_selected",
                           source = "somGrid")
        ),
        {
            rs <- shiny::isolate(rsUsed())
            if (is.null(rs)) return(NULL)
            d <- .safeEventData(
                verbose = verbose,
                "plotly_selected",
                source = "somGrid"
            )
            if (is.null(d)) return(NULL)
            d <- .inputSelect(d, rs, shiny::isolate(input$selectMode))
            shiny::isolate(rsUsed(d))
        }
    )

    # Scatter-plot rectangular selection observer.
    shiny::observeEvent(
        suppressWarnings(
            .safeEventData(verbose = verbose, "plotly_selected",
                           source = "scatterPlot")
        ),
        {
            rs <- shiny::isolate(rsUsed_d())
            dimSel <- dimSelection()
            sampleIds <- shiny::isolate(input$samples2plot)
            d <- .safeEventData(
                verbose = verbose,
                "plotly_selected",
                source = "scatterPlot"
            )

            plotIdx <- activePlot()
            if (is.null(d)) return(NULL)
            shiny::req(d$curveNumber)
            d <- d[d$curveNumber == 1L, , drop = FALSE]
            if (nrow(d) == 0L) return(NULL)
            minx <- min(d$x)
            maxx <- max(d$x)
            miny <- min(d$y)
            maxy <- max(d$y)
            df_local <- dfPlot()
            dims <- dimSel[[plotIdx]]$dims
            ids <- which(
                df_local[, make.names(dims[1L]), drop = TRUE] > minx &
                    df_local[, make.names(dims[1L]), drop = TRUE] < maxx &
                    df_local[, make.names(dims[2L]), drop = TRUE] > miny &
                    df_local[, make.names(dims[2L]), drop = TRUE] < maxy
            )
            ids <- SingleCellExperiment::colData(
                sce_subsampled[, sce_subsampled$sample_id %in% sampleIds]
            )[ids, "cluster_id"]
            clicked <- unique(as.integer(ids))
            mode <- shiny::isolate(input$selectMode)
            new_rs <- .applySelectMode(rs, clicked, mode)
            shiny::isolate(rsUsed(new_rs))
        }
    )

    # FlowSOM star plot selection observers (only when fsom is available).
    if (!is.null(fsom)) {
        shiny::observeEvent(
            .safeEventData(verbose = verbose, "plotly_selected",
                           source = "flowSOMStars"),
            {
                d <- .safeEventData(
                    verbose = verbose,
                    "plotly_selected",
                    source = "flowSOMStars"
                )
                if (is.null(d) || nrow(d) == 0L || is.null(d$customdata)) {
                    return(NULL)
                }
                clicked <- as.integer(d$customdata)
                rs <- shiny::isolate(rsUsed_d())
                mode <- shiny::isolate(input$selectMode)
                new_rs <- .applySelectMode(rs, clicked, mode)
                shiny::isolate(rsUsed(new_rs))
            }
        )

        shiny::observeEvent(
            .safeEventData(verbose = verbose, "plotly_click",
                           source = "flowSOMStars"),
            {
                d <- .safeEventData(
                    verbose = verbose,
                    "plotly_click",
                    source = "flowSOMStars"
                )
                if (is.null(d) || is.null(d$customdata)) return(NULL)
                clicked <- as.integer(d$customdata)
                rs <- shiny::isolate(rsUsed_d())
                mode <- shiny::isolate(input$selectMode)
                new_rs <- .applySelectMode(rs, clicked, mode)
                shiny::isolate(rsUsed(new_rs))
            }
        )
    }

    # Close the app and return the current outputList.
    shiny::observeEvent(input$close, {
        shinyjs::js$closeWindow()
        shiny::stopApp(rv$outputList)
    })
}
