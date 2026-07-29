# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server_clusterSelector.R ----
# Server factory for the clusterSelector Shiny app.
# This file builds the `server` function that is passed to shiny::shinyApp().

#' Build clusterSelector Server
#'
#' Creates the Shiny server function for the clusterSelector app.
#'
#' @param sce Full \code{SingleCellExperiment}.
#' @param sce_subsampled Subsampled \code{SingleCellExperiment}.
#' @param outputList Named list of cluster groupings.
#' @param colTree Optional collapsible tree object.
#' @param dList List of marker pairs for 2D plots.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param nPlots Number of 2D SOM plots to display.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somRasterObj Raster object for SOM visualization.
#' @param fsom Optional \code{FlowSOM} object.
#' @param env Environment used to store mutable state (legacy argument).
#' @param verbose Logical indicating whether to show detailed runtime messages.
#'
#' @return A \code{shiny} server function.
#'
#' @keywords internal
.buildClusterSelectorServer <- function(sce,
                                        sce_subsampled,
                                        outputList,
                                        colTree,
                                        dList,
                                        dend,
                                        dendTable,
                                        clusterPatientTable,
                                        somCodesName,
                                        nPlots,
                                        somRasterData,
                                        somRasterObj,
                                        fsom,
                                        env,
                                        verbose) {
    function(input, output, session) {
    # nocov start
        rnSCE <- rownames(sce)
        metaD <- S4Vectors::metadata(sce)
        colsUsed <- metaD$map$colsUsed

        ## Reactive values -----------------------------------------------------
        rsUsed <- shiny::reactiveVal(c(1))
        triggerRedraw <- shiny::reactiveVal(1)
        selectedPoints <- shiny::reactiveVal(NULL)
        activePlot <- shiny::reactiveVal(1)
        dListRV <- shiny::reactiveVal(dList)
        selectedUpdate2 <- shiny::reactiveVal(value = 0)
        selectedUpdate <- shiny::reactiveVal(value = 0)
        rv <- shiny::reactiveValues(outputList = outputList)

        ## Static / cached data ------------------------------------------------
        df <- as.data.frame(t(SummarizedExperiment::assays(sce_subsampled)[[1]]))
        colnames(df) <- make.names(rownames(sce_subsampled))
        dfPlot <- shiny::reactiveVal(df)

        rN <- colnames(metaD[[somCodesName]])
        nNsub <- colnames(S4Vectors::metadata(sce_subsampled)[[somCodesName]])
        sce_subsampledRN <- rownames(sce_subsampled)
        sceRN <- rownames(sce)
        sceCN <- colnames(sce)

        # Pre-compute per-channel axis limits once from the full assay matrix.
        assay_mat <- SummarizedExperiment::assays(sce)[[1]]
        channelLimits <- lapply(rownames(sce), function(ch) {
            vals <- assay_mat[ch, ]
            c(min = min(vals), max = max(vals))
        })
        names(channelLimits) <- rownames(sce)

        # Pre-build the SOM raster heatmap as a ggplot object once.
        baseRasterGgplot <- .buildBaseRasterGgplot(somRasterData, colsUsed)

        # in-session base-plot cache: avoids re-running plotSOMScatter when only
        # the selection changes. Key = ch1 + ch2 + colorVar + sizeVar.
        .bpCache <- new.env(parent = emptyenv())
        getBasePlot <- function(ch1, ch2, colorVar, sizeVar) {
            key <- paste(ch1, ch2, colorVar, sizeVar, sep = "\x01")
            if (exists(key, envir = .bpCache, inherits = FALSE)) {
                return(.bpCache[[key]])
            }
            pp1 <- plotSOMScatter(
                x         = sce,
                chs       = c(ch1, ch2),
                pointSize = sizeVar,
                color_by  = colorVar,
                zeros     = TRUE,
                xRN       = sceRN,
                xCN       = sceCN
            )
            .bpCache[[key]] <- pp1
            pp1
        }

        updatedoutputList <- function() {
            .updateOutputListInputs(session, input, rv$outputList, metaD)
        }

        zoomFunc <- function(zoom, plotIdx) {
            dimSelectionInternal <- dimSelection()
            shiny::req(dimSelectionInternal)
            if (all(!is.null(zoom), "xaxis.range[0]" %in% names(zoom), na.rm = TRUE)) {
                rezoom <- FALSE
                if (all(c(
                    dimSelectionInternal[[plotIdx]]$xzoom[1] > zoom$`xaxis.range[0]`,
                    !is.null(dimSelectionInternal[[plotIdx]]$xzoom[1])
                ), na.rm = TRUE)) {
                    rezoom <- TRUE
                }
                if (rezoom) {
                    dimSelectionInternal[[plotIdx]]$xzoom <- c(NULL, NULL)
                    dimSelectionInternal[[plotIdx]]$yzoom <- c(NULL, NULL)
                } else {
                    dimSelectionInternal[[plotIdx]]$xzoom <- c(zoom$`xaxis.range[0]`, zoom$`xaxis.range[1]`)
                    dimSelectionInternal[[plotIdx]]$yzoom <- c(zoom$`yaxis.range[0]`, zoom$`yaxis.range[1]`)
                }
                dimSelection(dimSelectionInternal)
                shiny::isolate(triggerRedraw(triggerRedraw() + 1))
            }
        }

        # activeDims ------------------------------------------------------------
        activeDims <- shiny::reactive({
            x <- input$currentDimX
            y <- input$currentDimY
            shiny::req(x, y)
            c(x, y)
        }) %>% shiny::debounce(250)

        dimSelection <- shiny::reactiveVal(list())
        shiny::observe({
            dims <- activeDims()
            # message("DIAG[dimSelection] activeDims() = ", if (is.null(dims)) "NULL" else paste(dims, collapse = ", "))
            shiny::req(dims[1L], dims[2L])
            lim1 <- channelLimits[[dims[1L]]]
            lim2 <- channelLimits[[dims[2L]]]
            shiny::req(lim1, lim2)
            existing <- shiny::isolate(dimSelection())
            # message("DIAG[dimSelection] existing length = ", length(existing),
            #         " | existing[[1]]$dims = ", if (length(existing) >= 1L) paste(existing[[1L]]$dims, collapse = ", ") else "none",
            #         " | identical to new dims? ", if (length(existing) >= 1L) identical(existing[[1L]]$dims, dims) else "n/a")
            if (length(existing) >= 1L && identical(existing[[1L]]$dims, dims)){
                # message("DIAG[dimSelection] EARLY RETURN -- dims unchanged, dimSelection() not updated")
                return()
            }
            dimSelection(list(list(
                dims  = dims,
                xlim  = c(lim1["min"], lim1["max"]),
                ylim  = c(lim2["min"], lim2["max"]),
                xzoom = c(NULL, NULL),
                yzoom = c(NULL, NULL)
            )))
            # message("DIAG[dimSelection] wrote new dimSelection() for dims = ", paste(dims, collapse = ", "))

        })

        # Dimension-reduction reactives ----------------------------------------
        dimRedSelection <- shiny::reactive({
            input$dimRedSelection
        }) %>% shiny::debounce(1000)

        tsne <- shiny::reactive({
            dimRedCols <- dimRedSelection()
            perplexity <- input$perplexity
            tsneFunc(dimRedSelection = dimRedCols, perplexity = perplexity, sce, somCodesName)
        }) %>%
            shiny::bindCache(dimRedSelection(), input$perplexity) %>%
            shiny::debounce(1000)

        umap <- shiny::reactive({
            dimRedCols <- input$dimRedSelection
            pumap <- umap::umap.defaults
            pumap$n_neighbors <- input$n_neighbors
            umap::umap(metaD[[somCodesName]][, dimRedCols], config = pumap)
        }) %>%
            shiny::bindCache(input$dimRedSelection, input$n_neighbors) %>%
            shiny::debounce(1000)

        pca <- shiny::reactive({
            dimRedCols <- input$dimRedSelection
            stats::prcomp(t(metaD[[somCodesName]][, dimRedCols]), scale = FALSE, rank. = 2)
        }) %>%
            shiny::bindCache(input$dimRedSelection) %>%
            shiny::debounce(1000)

        # Plot-producing reactives ---------------------------------------------
        somBasePlot <- shiny::reactive({
            dims     <- activeDims()
            colorVar <- input$somColorVar
            sizeVar  <- input$somSizeVar
            if (is.null(colorVar)) colorVar <- "n"
            if (is.null(sizeVar))  sizeVar  <- "max"
            shiny::req(dims)
            plotSOMScatter(
                x         = sce,
                chs       = dims,
                pointSize = sizeVar,
                color_by  = colorVar,
                zeros     = TRUE,
                xRN       = sceRN,
                xCN       = sceCN
            )
        }) %>% shiny::bindCache(activeDims(), input$somColorVar, input$somSizeVar)

        somProjectionDf <- shiny::reactive({
            pp1    <- somBasePlot()
            dimSel <- dimSelection()
            shiny::req(pp1, length(dimSel) >= 1L)
            buildProjectionDf(pp1, 1L, dimSel, sce, somCodesName = somCodesName)
        }) %>% shiny::bindCache(activeDims())

        tsnePlot <- shiny::reactive({
            selectedUpdate2()
            rs <- rsUsed_d()
            shiny::req(rs)
            triggerRedraw()
            tsneObj <- tsne()
            df <- as.data.frame(tsneObj$Y)
            names(df) <- c("tsne1", "tsne2")
            drawProjection(df, rs, colorbyGroups = input$colorbyGroups,
                           sce = sce, outputList = rv$outputList)
        })

        umapPlot <- shiny::reactive({
            umapObj <- umap()
            selectedUpdate2()
            rs <- rsUsed_d()
            shiny::req(rs)
            triggerRedraw()
            df <- as.data.frame(umapObj$layout)
            names(df) <- c("umap1", "umap2")
            drawProjection(df, rs, colorbyGroups = input$colorbyGroups,
                           sce = sce, outputList = rv$outputList)
        })

        pcaPlot <- shiny::reactive({
            pcaObj <- pca()
            selectedUpdate2()
            rs <- rsUsed_d()
            shiny::req(rs)
            triggerRedraw()
            shiny::req(pcaObj)
            df <- as.data.frame(pcaObj$rotation)
            colnames(df) <- c("pc1", "pc2")
            drawProjection(df, rs, colorbyGroups = input$colorbyGroups,
                           sce = sce, outputList = rv$outputList)
        })

        scatterPlot <- shiny::reactive({
            rs <- rsUsed()
            shiny::req(rs)
            dimSel <- dimSelection()
            sampleIds <- input$samples2plot
            shiny::req(sampleIds)
            plotIdx <- activePlot()
            shiny::req(dimSel, length(dimSel) >= 1L)
            cidIdx <- colData(sce_subsampled)$cluster_id %in% rs &
                colData(sce_subsampled)$sample_id %in% sampleIds
            if (sum(cidIdx) < 1) return(NULL)
            chs <- dimSel[[1L]]$dims
            chx <- make.names(chs[1])
            chy <- make.names(chs[2])
            xVals <- df[cidIdx, chx, drop = TRUE]
            yVals <- df[cidIdx, chy, drop = TRUE]
            pctl <- input$scatterPercentile
            if (is.null(pctl)) pctl <- 0.99
            tailP <- (1 - pctl) / 2
            xlimP <- stats::quantile(xVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
            ylimP <- stats::quantile(yVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
            pp <- plotCytoScatter(
                rownms = sce_subsampledRN,
                x      = sce_subsampled[, cidIdx],
                chs    = chs,
                bins   = 200
            ) +
                ggplot2::xlim(xlimP) +
                ggplot2::ylim(ylimP)
            pp
        })

        somRasterPlot <- shiny::reactive({
            rs <- as.integer(unlist(rsUsed()))
            rs <- rs[!is.na(rs) & rs > 0L]
            shiny::req(length(rs) > 0L)
            somRasterData[somRasterData$id %in% rs, c("x", "y"), drop = FALSE]
        })

        dendPlot <- shiny::reactive({
            rs <- rsUsed()
            shiny::req(rs)
            selectedPoints()
            labCol <- rep("blue", length(labels(dend)))
            labCol[which(labels(dend) %in% rs)] <- "red"
            dend %>%
                dendextend::set("leaves_col", labCol) %>%
                dendextend::set("leaves_pch", 15) %>%
                dendextend::set("leaves_cex", 1) %>%
                dendextend::set("labels_cex", 0.5)
        })

        dendPlotlyData <- shiny::reactive({
            g      <- dendextend::as.ggdend(dend)
            labels <- g$labels
            if (!is.null(dendTable) && all(c("id", "label") %in% colnames(dendTable))) {
                labels <- dplyr::left_join(
                    labels,
                    dendTable[, c("id", "label")],
                    by = "label"
                )
            } else {
                labels$id <- suppressWarnings(as.integer(as.character(labels$label)))
            }
            list(segments = g$segments, labels = labels)
        })

        countBarPlot <- shiny::reactive({
            cst <- input$compareStatsTo
            outputList <- rv$outputList
            rs <- rsUsed()
            shiny::req(rs)
            grpInp <- groupsInput()
            shiny::req(clusterPatientTable)
            countBarPlotFunc(rs, clusterPatientTable, cst, sce, outputList, grpInp)
        })

        PercentBarPlot <- shiny::reactive({
            outputList <- rv$outputList
            rs <- rsUsed()
            shiny::req(rs)
            shiny::req(clusterPatientTable)
            relativeToCol <- input$relativeTo
            grpInp <- groupsInput()
            PercentBarPlotFunc(sce, relativeToCol, clusterPatientTable, rs,
                               outputList, group, grpInp)
        })

        groupsInput <- shiny::reactive({
            groupsVar <- input$groupsVar
            if (is.null(groupsVar) || length(groupsVar) == 0 ||
                !groupsVar %in% colnames(metaD$experiment_info)) {
                return(NULL)
            }
            grp1 <- metaD$experiment_info[
                metaD$experiment_info[, groupsVar] %in% input$group1, "sample_id"
            ]
            grp2 <- metaD$experiment_info[
                metaD$experiment_info[, groupsVar] %in% input$group2, "sample_id"
            ]
            list(group1 = grp1, group2 = grp2)
        }) %>% shiny::debounce(1000)

        inputClusterNumber <- shiny::reactive({
            txt <- input$clusterNumbers
            if (is.null(txt) || !nzchar(txt)) return(integer(0))
            vals <- stringr::str_split(txt, ",")[[1]] %>%
                trimws() %>%
                as.integer()
            vals[!is.na(vals)]
        }) %>% shiny::debounce(1000)


        violinPlotSelection <- shiny::reactive({
            input$violinSelection
        }) %>% shiny::debounce(1000)

        vlnPlot <- shiny::reactive({
            input$applyName
            input$rmGrp
            violinSel <- violinPlotSelection()
            upsetSel <- input$upsetSelection
            outputList <- rv$outputList
            shiny::req(outputList)
            markers <- colnames(sce@metadata[[somCodesName]])
            if (length(upsetSel) < 3) upsetSel <- names(outputList)
            markers <- intersect(markers, violinSel)
            somCodes <- sce@metadata[[somCodesName]]
            nNodes <- nrow(somCodes)
            parts <- lapply(upsetSel, function(na) {
                rows <- outputList[[na]]
                rows <- rows[rows != 0]
                if (length(rows) < 2) return(NULL)
                wide <- as.data.frame(somCodes[rows, markers, drop = FALSE])
                wide$somNode <- factor(rows, levels = seq_len(nNodes))
                long <- tidyr::pivot_longer(wide, cols = tidyselect::all_of(markers),
                                            names_to = "marker",
                                            values_to = "expr")
                long$grpName <- na
                long
            })
            data <- data.table::rbindlist(parts)
            if (nrow(data) == 0) {
                return(
                    ggplot2::ggplot() +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "No groups with >= 2 SOM nodes")
                )
            }
            data$marker <- factor(data$marker, levels = markers)
            data$grpName <- factor(data$grpName, levels = upsetSel)
            nb.cols <- length(unique(data$marker))
            mycolors <- grDevices::colorRampPalette(
                RColorBrewer::brewer.pal(8, "Set2")
            )(nb.cols)
            ggplot2::ggplot(data, ggplot2::aes(factor(marker), expr, fill = marker)) +
                ggplot2::geom_violin(scale = "width", adjust = 1, trim = TRUE) +
                ggplot2::scale_y_continuous(
                    expand = c(0, 0), position = "right",
                    labels = function(x) {
                        c(rep(x = "", times = length(x) - 2),
                          x[length(x) - 1], "")
                    }
                ) +
                ggplot2::facet_grid(rows = ggplot2::vars(grpName),
                                     scales = "free", switch = "y") +
                cowplot::theme_cowplot(font_size = 12) +
                ggplot2::theme(
                    legend.position = "none",
                    panel.spacing = grid::unit(0, "lines"),
                    plot.title = ggplot2::element_text(hjust = 0.5),
                    panel.background = ggplot2::element_rect(fill = NA, color = "black"),
                    strip.background = ggplot2::element_blank(),
                    strip.text = ggplot2::element_text(face = "bold"),
                    strip.text.y.left = ggplot2::element_text(angle = 0),
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
                ) +
                ggplot2::scale_fill_manual(values = mycolors) +
                ggplot2::ggtitle("Marker on x-axis") +
                ggplot2::xlab("Marker") +
                ggplot2::ylab("Expression Level")
        })

        VlnPlot2 <- shiny::reactive({
            input$applyName
            input$rmGrp
            upsetSel <- input$upsetSelection
            violinSel <- violinPlotSelection()
            outputList <- rv$outputList
            shiny::req(outputList)
            plotViolin2Func(sce, somCodesName, violinSel, upsetSel, outputList)
        })

        upSetPlot <- shiny::reactive({
            selectedUpdate()
            input$applyName
            upsetSel <- input$upsetSelection
            input$rmGrp
            outputList <- rv$outputList
            upsetPlotFunc(upsetSel, outputList, sce)
        })

        flowSOMPiePlot <- shiny::reactive({
            rs <- rsUsed()
            shiny::req(rs)
            ncol_val <- input$flowSOMPieCols
            maxPies <- input$flowSOMPieMax
            .buildFlowSOMPiePlot(metaD[[somCodesName]], rs, colsUsed, ncol = ncol_val, maxPies = maxPies)
        })

        rsUsed_d <- rsUsed %>% shiny::debounce(1500)
        sample2PlotDb <- shiny::reactive(input$samples2plot) %>% shiny::debounce(500)

        # Verbose logging observers (only when verbose = TRUE) -------------------
        if (isTRUE(verbose)) {
            # Log when dimension reduction computations happen
            shiny::observe({
                tsne()
                cat(sprintf("[CySA] t-SNE computed at %s\n", format(Sys.time(), "%H:%M:%S")), file = stderr())
            }) %>% shiny::bindCache(dimRedSelection(), input$perplexity)

            shiny::observe({
                umap()
                cat(sprintf("[CySA] UMAP computed at %s\n", format(Sys.time(), "%H:%M:%S")), file = stderr())
            }) %>% shiny::bindCache(input$dimRedSelection, input$n_neighbors)

            shiny::observe({
                pca()
                cat(sprintf("[CySA] PCA computed at %s\n", format(Sys.time(), "%H:%M:%S")), file = stderr())
            }) %>% shiny::bindCache(input$dimRedSelection)

            # Log when SOM base plot is rebuilt
            shiny::observe({
                somBasePlot()
                dims <- activeDims()
                cat(sprintf("[CySA] SOM base plot rebuilt (%s vs %s) at %s\n",
                            dims[1], dims[2], format(Sys.time(), "%H:%M:%S")), file = stderr())
            }) %>% shiny::bindCache(activeDims(), input$somColorVar, input$somSizeVar)
        }

        # Register observers and outputs -----------------------------------------
        # Capture env explicitly so observers/outputs can write back to it even
        # when the enclosing factory environment is not on the search path (e.g.
        # during shiny::testServer()).
        env_local <- env

        .registerClusterSelectorObservers(
            input, output, session,
            rv, rsUsed, rsUsed_d, activePlot, dListRV, dimSelection,
            triggerRedraw, selectedUpdate2, selectedUpdate, updatedoutputList,
            sample2PlotDb, dfPlot, groupsInput, inputClusterNumber,
            zoomFunc, sce, sce_subsampled, metaD, dend, dendTable,
            clusterPatientTable, somCodesName, colsUsed, nPlots, verbose,
            fsom = fsom, env = env_local
        )

        .registerClusterSelectorOutputs(
            input, output, session,
            rv, rsUsed, rsUsed_d, activePlot, dListRV, dimSelection,
            triggerRedraw, selectedUpdate2, selectedUpdate, updatedoutputList,
            colTree, somRasterData, somRasterObj, baseRasterGgplot,
            somRasterPlot, somBasePlot, somProjectionDf,
            tsnePlot, umapPlot, pcaPlot, scatterPlot,
            dendPlot, dendPlotlyData, countBarPlot, PercentBarPlot,
            vlnPlot, VlnPlot2, upSetPlot, flowSOMPiePlot, groupsInput,
            channelLimits, getBasePlot, fsom,
            sce, sce_subsampled, sce_subsampledRN, sceRN, sceCN,
            metaD, dend, dendTable, clusterPatientTable, somCodesName,
            colsUsed, rnSCE, nPlots, df, dfPlot, env = env_local, verbose
        )
    }
    # nocov end
}
