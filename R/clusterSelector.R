# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

utils::globalVariables(c(
    ".", ".data", "..ncount..", "counts", "groups", "id", "marker", "expr",
    "grpName", "Percent", "n", "value", "variable", "cluster", "colGrp",
    "x", "y", "xend", "yend", "customdata", "label", "tsne1", "tsne2", "umap1", "umap2",
    "pc1", "pc2", "rowElement", "cellCounts", "cellPercentages",
    "dend", "dendTable", "clusterPatientTable", "somRasterData",
    "somRasterObj", "sce_subsampled", "metaD", "sceRN", "sceCN",
    "chx", "nNsub", "sce_subsampledRN", "df", "dfPlot",
    "outputList", "dimSelection", "activePlot", "rsUsed",
    "inputClusterNumber", "violinPlotSelection", "groupsInput",
    "countBarPlot", "PercentBarPlot", "dendPlot", "tsnePlot",
    "umapPlot", "pcaPlot", "scatterPlot", "somRasterPlot", "flowSOMPiePlot",
    "vlnPlot", "VlnPlot2", "upSetPlot", "selectedUpdate", "selectedUpdate2",
    "sample2PlotDb", "choicesRV",
    "sample_id", "group", "val", "somNode", "N", "thrdQu",
    "scatterPercentile", "somColorVar", "somSizeVar"
))

#' Cluster Selector Shiny Application
#'
#' Creates an interactive Shiny dashboard for selecting and visualizing
#' clusters from a SingleCellExperiment object.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing the full dataset.
#' @param sce_subsampled A subsampled \code{SingleCellExperiment} for
#'   performance-sensitive plots.
#' @param outputList A named list of cluster groupings.
#' @param colTree Optional collapsible tree object.
#' @param dList List of marker pairs for 2D plots.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param nPlots Number of 2D SOM plots to display.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somRasterObj Raster object for SOM visualization. If \code{NULL},
#'   the SOM grid dimensions are inferred from \code{somRasterData}.
#' @param fsom Optional \code{FlowSOM} object. If provided, an interactive
#'   FlowSOM star plot is shown.
#' @param env Environment used to store mutable state (legacy argument).
#' @param verbose Logical indicating whether to show detailed runtime messages.
#'
#' @return A \code{\link[shiny]{shinyApp}} object. Launch the app with
#'   \code{shiny::runApp(app)}.
#'
#' @examples
#' sce <- CySA_example_sce()
#' prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
#' som_codes <- S4Vectors::metadata(sce)$SOM_codes
#' dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
#' dendTable <- data.frame(
#'   id = seq_len(nrow(som_codes)),
#'   label = rownames(som_codes),
#'   stringsAsFactors = FALSE
#' )
#' clusterPatientTable <- table(
#'   sample_id = sce$sample_id,
#'   cluster_id = sce$cluster_id
#' )
#' markers <- S4Vectors::metadata(sce)$map$colsUsed
#' somRasterData <- data.frame(
#'   x = rep(seq_len(5), length.out = nrow(som_codes)),
#'   y = rep(seq_len(2), each = nrow(som_codes) / 2),
#'   id = seq_len(nrow(som_codes))
#' )
#' for (m in markers) {
#'   somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
#' }
#' arr <- array(
#'   data = seq_len(10 * 10 * length(markers)),
#'   dim = c(10, 10, length(markers))
#' )
#' somRasterObj <- raster::brick(arr)
#' names(somRasterObj) <- markers
#'
#' app <- clusterSelector(
#'   sce = prepped$sce,
#'   sce_subsampled = prepped$sce_subsampled,
#'   dList = prepped$dList,
#'   dend = dend,
#'   dendTable = dendTable,
#'   clusterPatientTable = clusterPatientTable,
#'   somRasterData = somRasterData,
#'   somRasterObj = somRasterObj
#' )
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#'
#' @export
clusterSelector <- function(sce, # main input has to contain:
                            sce_subsampled, # subsampled sce object
                            outputList = list(), # list of named nodes
                            colTree = NULL, # Tree object to plot
                            dList,
                            dend,
                            dendTable,
                            clusterPatientTable,
                            somCodesName = "SOM_codes", # SOM_codes.1
                            nPlots = 6,
                            somRasterData,
                            somRasterObj,
                            fsom = NULL,
                            env = environment(),
                            verbose = FALSE) {
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
    outputList <- .initializeOutputList(outputList, levels(sce$cluster_id))
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

    ## UI ----
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

    #  server  ----
    server <- function(input, output, session) {
        rsUsed <- reactiveVal(c(1))
        triggerRedraw <- reactiveVal(1)
        selectedPoints <- reactiveVal(NULL)
        activePlot <- reactiveVal(1)
        df <- as.data.frame(t(assays(sce_subsampled)[[1]]))
        colnames(df) <- make.names(rownames(sce_subsampled))
        dfPlot <- shiny::reactiveVal(df)
        # dfall = data.frame(t(assays(sce)[[1]]))
        # rN = rownames(sce)
        metaD <- S4Vectors::metadata(sce)
        rN <- colnames(metaD[[somCodesName]])
        nNsub <- colnames(S4Vectors::metadata(sce_subsampled)[[somCodesName]])
        rownames(sce_subsampled)
        sce_subsampledRN <- rownames(sce_subsampled)
        sceRN <- rownames(sce)
        sceCN <- colnames(sce)
        chx <- channels(sce)
        if (!is.null("colTree")) {
            output$plot <- collapsibleTree::renderCollapsibleTree({
                colTree
            })
        }
        # activeDims ----
        # Debounced to absorb the brief window when a preset-pair change sends two
        # updateSelectInput messages (one for X, one for Y) in the same flush cycle.
        activeDims <- shiny::reactive({
            x <- input$currentDimX
            y <- input$currentDimY
            shiny::req(x, y)
            c(x, y)
        }) %>% shiny::debounce(250)
        ## ── dListRV: reactive store for the pair list ─────────────────────────────
        dListRV <- shiny::reactiveVal(dList)

        ## ── in-session base-plot cache ────────────────────────────────────────────
        ## Avoids re-running plotSOMScatter (the expensive part) when only rs changes.
        ## Key: ch1 + ch2 + colorVar + sizeVar.  Lives for the Shiny session.
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

        ## ── helper: rebuild remove-picker choices from a pair list ────────────────
        .pairChoices <- function(lst) {
            if (length(lst) == 0L) return(character(0L))
            labels <- vapply(lst,
                             function(p) paste(p[1L], "\u00b7", p[2L]), character(1L))
            stats::setNames(seq_along(lst), labels)
        }

        ## ── dynamic preset picker in the SOM 2D plots box ─────────────────────────
        output$dimPairSelectUI <- shiny::renderUI({
            current <- dListRV()
            if (length(current) == 0L) return(NULL)
            pair_labels <- vapply(current,
                                  function(p) paste(p[1L], "\u00b7", p[2L]), character(1L))
            shiny::selectInput(
                inputId   = "dimPairSelect",
                label     = "Preset pair",
                choices   = stats::setNames(seq_along(current), pair_labels),
                selected  = "1",
                multiple  = FALSE,
                selectize = FALSE,   # native <select>: no overflow in narrow column
                width     = "100%"
            )
        })

        ## ── preset pair → update axis selects ────────────────────────────────────
        shiny::observeEvent(input$dimPairSelect, {
            idx     <- suppressWarnings(as.integer(input$dimPairSelect))
            current <- dListRV()
            if (is.na(idx) || idx < 1L || idx > length(current)) return(NULL)
            pair <- current[[idx]]
            shiny::updateSelectInput(session, "currentDimX", selected = pair[1L])
            shiny::updateSelectInput(session, "currentDimY", selected = pair[2L])
        }, ignoreInit = TRUE, ignoreNULL = TRUE)

        ## ── Add pair ──────────────────────────────────────────────────────────────
        shiny::observeEvent(input$staticSomAdd, {
            shiny::req(input$staticSomX, input$staticSomY)
            new_pair <- c(input$staticSomX, input$staticSomY)
            current  <- dListRV()

            ## Reject exact duplicates (same X and Y in same order)
            is_dup <- any(vapply(current,
                                 function(p) identical(p, new_pair), logical(1L)))
            if (is_dup) {
                shiny::showNotification(
                    paste("Pair", new_pair[1L], "\u00b7", new_pair[2L], "already exists."),
                    type = "warning", duration = 3L
                )
                return(NULL)
            }

            new_list <- c(current, list(new_pair))
            dListRV(new_list)

            shiny::updateSelectInput(session, "staticSomRemove",
                                     choices  = .pairChoices(new_list),
                                     selected = as.character(length(new_list))   # select the new entry
            )
        })

        ## ── Remove pair ───────────────────────────────────────────────────────────
        shiny::observeEvent(input$staticSomRemoveBtn, {
            current <- dListRV()
            idx     <- suppressWarnings(as.integer(input$staticSomRemove))
            if (is.na(idx) || idx < 1L || idx > length(current)) return(NULL)

            new_list <- current[-idx]
            dListRV(new_list)

            new_sel <- as.character(max(1L, idx - 1L))
            shiny::updateSelectInput(session, "staticSomRemove",
                                     choices  = .pairChoices(new_list),
                                     selected = new_sel
            )
        })

        ## ── Dynamic plot grid UI ──────────────────────────────────────────────────
        output$staticSomGrid <- shiny::renderUI({
            current <- dListRV()
            n_cols  <- as.integer(input$staticSomCols)
            if (is.na(n_cols) || n_cols < 1L) n_cols <- 3L
            height  <- paste0(
                if (is.null(input$staticSomHeight)) 240L else input$staticSomHeight,
                "px"
            )

            if (length(current) == 0L) {
                return(htmltools::tags$p(
                    "No pairs configured. Use the controls above to add pairs.",
                    style = "color:#aaa; font-style:italic; padding:12px;"
                ))
            }

            col_width <- 12L %/% n_cols

            panels <- lapply(seq_along(current), function(i) {
                pair <- current[[i]]
                shiny::column(
                    width = col_width,
                    htmltools::tags$p(
                        paste(pair[1L], "\u00b7", pair[2L]),
                        style = paste(
                            "font-size:11px; font-weight:600;",
                            "text-align:center; margin:4px 0 2px; color:#555;"
                        )
                    ),
                    shiny::plotOutput(paste0("staticSomDyn", i), height = height)
                )
            })

            ## Chunk into rows of n_cols
            row_starts <- seq(1L, length(panels), by = n_cols)
            rows <- lapply(row_starts, function(s) {
                do.call(shiny::fluidRow,
                        panels[s:min(s + n_cols - 1L, length(panels))])
            })
            do.call(htmltools::tagList, rows)
        })

        ## ── Dynamic renderPlot registration ──────────────────────────────────────
        ## Re-runs only when dListRV() changes (pair added / removed).
        ## Re-registering an existing output ID is safe in Shiny (replaces renderer).
        ## Changing staticSomCols/Height only affects the UI grid, not this observer.
        shiny::observe({
            current <- dListRV()

            lapply(seq_along(current), function(idx) {
                local({
                    i      <- idx
                    pair_i <- current[[i]]
                    out_id <- paste0("staticSomDyn", i)

                    output[[out_id]] <- shiny::renderPlot({
                        rs       <- rsUsed_d()
                        shiny::req(rs)
                        colorVar <- if (is.null(input$somColorVar)) "n"   else input$somColorVar
                        sizeVar  <- if (is.null(input$somSizeVar))  "max" else input$somSizeVar
                        pctl     <- if (is.null(input$scatterPercentile)) 0.99 else input$scatterPercentile

                        ## Expensive base scatter — served from cache when args unchanged
                        pp1 <- getBasePlot(pair_i[1L], pair_i[2L], colorVar, sizeVar)
                        shiny::req(pp1)

                        lim1 <- channelLimits[[pair_i[1L]]]
                        lim2 <- channelLimits[[pair_i[2L]]]
                        if (is.null(lim1)) lim1 <- c(min = 0, max = 1)
                        if (is.null(lim2)) lim2 <- c(min = 0, max = 1)

                        dimSel_i <- list(list(
                            dims  = pair_i,
                            xlim  = c(lim1["min"], lim1["max"]),
                            ylim  = c(lim2["min"], lim2["max"]),
                            xzoom = c(NULL, NULL),
                            yzoom = c(NULL, NULL)
                        ))

                        tailP     <- (1 - pctl) / 2
                        somCodes  <- metaD[[somCodesName]]
                        xlimP <- if (pair_i[1L] %in% colnames(somCodes))
                            stats::quantile(somCodes[, pair_i[1L]],
                                            probs = c(tailP, 1 - tailP), na.rm = TRUE)
                        else NULL
                        ylimP <- if (pair_i[2L] %in% colnames(somCodes))
                            stats::quantile(somCodes[, pair_i[2L]],
                                            probs = c(tailP, 1 - tailP), na.rm = TRUE)
                        else NULL

                        ## Cheap: adds red overlay to cached base plot
                        ggsomPlot(pp1, 1L, rs, dimSel_i,
                                  sce   = sce,
                                  metaD = metaD,
                                  xlim  = xlimP,
                                  ylim  = ylimP)
                    })

                    ## Suspend while box is collapsed — no renders until user opens it
                    shiny::outputOptions(output, out_id, suspendWhenHidden = TRUE)
                })
            })
        })
        # obs dimPairSelect ----
        # When the user clicks a preset radio button, sync the X/Y selects.
        # ignoreInit = TRUE keeps startup free of a redundant update.
        shiny::observeEvent(input$dimPairSelect, {
            idx <- suppressWarnings(as.integer(input$dimPairSelect))
            if (is.na(idx) || idx < 1L || idx > length(dList)) return(NULL)
            pair <- dList[[idx]]
            shiny::updateSelectInput(session, "currentDimX", selected = pair[1L])
            shiny::updateSelectInput(session, "currentDimY", selected = pair[2L])
        }, ignoreInit = TRUE, ignoreNULL = TRUE)

        # somBasePlot ----
        # Single cached SOM scatter base plot. Only recomputes when dims, color, or
        # size variable change — NOT when rs or colorbyGroups change.
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

        # somProjectionDf ----
        # Pre-join of SOM stats into the 2D coordinate frame.  Cached per dim pair so
        # that switching colorbyGroups or rs does not redo the ggplot_build + join.
        somProjectionDf <- shiny::reactive({
            pp1    <- somBasePlot()
            dimSel <- dimSelection()
            shiny::req(pp1, length(dimSel) >= 1L)
            buildProjectionDf(pp1, 1L, dimSel, sce,
                              somCodesName = somCodesName)   # <-- add this
        }) %>% shiny::bindCache(activeDims())

        # output$somDataMain ----
        output$somDataMain <- plotly::renderPlotly({
            colorbyGroups <- input$colorbyGroups
            selectedUpdate2()
            showGroups <- input$showGroups
            dimSel     <- dimSelection()
            shiny::req(length(dimSel) >= 1L)
            rs <- rsUsed_d()
            shiny::req(rs)
            triggerRedraw()

            dims <- dimSel[[1L]]$dims
            pp1  <- somBasePlot()
            pDf  <- if (isTRUE(showGroups)) somProjectionDf() else NULL

            pctl  <- input$scatterPercentile
            if (is.null(pctl)) pctl <- 0.99
            tailP     <- (1 - pctl) / 2
            somCodes  <- metaD[[somCodesName]]
            xlimP <- if (dims[1L] %in% colnames(somCodes))
                stats::quantile(somCodes[, dims[1L]], probs = c(tailP, 1 - tailP), na.rm = TRUE)
            else NULL
            ylimP <- if (dims[2L] %in% colnames(somCodes))
                stats::quantile(somCodes[, dims[2L]], probs = c(tailP, 1 - tailP), na.rm = TRUE)
            else NULL

            somPlot(
                pp1, 1L, rs, colorbyGroups, showGroups,
                dimSelection = dimSel,
                sce          = sce,
                metaD        = metaD,
                outputList   = rv$outputList,
                projectionDf = pDf,
                xlim         = xlimP,
                ylim         = ylimP,
                source       = "somDataMain"   # explicit source — no longer derived from plotIdx
            )
        })

        # obs somDataMain selection ----
        shiny::observeEvent(
            safe_event_data(verbose = verbose, "plotly_selected", source = "somDataMain"),
            {
                d <- safe_event_data(verbose = verbose, "plotly_selected", source = "somDataMain")
                if (is.null(d)) return(NULL)
                rs <- shiny::isolate(rsUsed_d())
                shiny::req(rs)
                d <- .inputSelect(d, rs, shiny::isolate(input$selectMode))
                shiny::isolate(rsUsed(d))
            }
        )

        # obs somDataMain zoom ----
        shiny::observe({
            zoom <- safe_event_data(
                verbose = verbose, "plotly_relayout", source = "somDataMain"
            )
            zoomFunc(zoom, 1L)
        })

        # force redraw if selected group is used
        selectedUpdate2 <- reactiveVal(value = 0)

        # reactiveValues replacement for external env$outputList
        rv <- reactiveValues(outputList = outputList)

        # Keep external env in sync for callers that still read env$outputList
        shiny::observe({
            assign(x = "outputList", value = rv$outputList, envir = env)
        })

        shiny::observe({
            rs <- rsUsed_d()
            if ("selected" %in% input$colorbyGroups) {
                isolate(selectedUpdate2(selectedUpdate2() + 1))
            }
        })

        ## flowSOMPie dynamic height ----
        ## Computes total plot height from (number of rows) × (pixels per pie).
        ## This keeps each pie a consistent size as the selection grows or shrinks.
        output$flowSOMPieUI <- shiny::renderUI({
            rs       <- rsUsed()
            n_nodes  <- max(1L, length(rs))
            n_cols   <- max(1L, as.integer(input$flowSOMPieCols))
            pie_px   <- max(50L, as.integer(
                if (is.null(input$flowSOMPieSize)) 150L else input$flowSOMPieSize
            ))

            n_rows       <- ceiling(n_nodes / n_cols)
            ## Add 60 px for the shared legend at the bottom of the ggplot
            total_height <- n_rows * pie_px + 60L

            shiny::plotOutput(
                "flowSOMPie",
                height = paste0(total_height, "px")
            ) %>% shinyjqui::jqui_resizable()
        })

        choicesRV <- reactiveValues(trigger = 1)

        observe({
            choicesRV$trigger
        })
        ## Six static SOM views ----
        ## Each pair gets its own cached base-plot reactive + a renderPlot output.
        ## outputOptions(suspendWhenHidden) means the plots are never computed while
        ## the box is collapsed, eliminating six ggsomPlot calls at startup.

            n_static <- min(6L, length(dList))

            ## Non-reactive per-pair dim-selection slots (axis limits fixed at startup).
            staticDimSelections <- lapply(seq_len(n_static), function(idx) {
                pair <- dList[[idx]]
                lim1 <- channelLimits[[pair[1L]]]
                lim2 <- channelLimits[[pair[2L]]]
                if (is.null(lim1)) lim1 <- c(min = 0, max = 1)
                if (is.null(lim2)) lim2 <- c(min = 0, max = 1)
                list(
                    dims  = pair,
                    xlim  = c(lim1["min"], lim1["max"]),
                    ylim  = c(lim2["min"], lim2["max"]),
                    xzoom = c(NULL, NULL),
                    yzoom = c(NULL, NULL)
                )
            })

            ## One cached base-plot reactive per pair.
            ## Cache key: pair index + color/size variables.
            ## Invalidates only when the user changes somColorVar or somSizeVar —
            ## NOT when rs or colorbyGroups change.
            staticBasePlots <- lapply(seq_len(n_static), function(idx) {
                shiny::reactive({
                    colorVar <- input$somColorVar
                    sizeVar  <- input$somSizeVar
                    if (is.null(colorVar)) colorVar <- "n"
                    if (is.null(sizeVar))  sizeVar  <- "max"
                    pair <- dList[[idx]]
                    plotSOMScatter(
                        x         = sce,
                        chs       = pair,
                        pointSize = sizeVar,
                        color_by  = colorVar,
                        zeros     = TRUE,
                        xRN       = sceRN,
                        xCN       = sceCN
                    )
                }) %>% shiny::bindCache(idx, input$somColorVar, input$somSizeVar)
            })

            ## renderPlot (static, not plotly) outputs — one per pair.
            lapply(seq_len(n_static), function(idx) {
                output[[paste0("staticSom", idx)]] <- shiny::renderPlot({
                    rs <- rsUsed_d()
                    shiny::req(rs)
                    pp1 <- staticBasePlots[[idx]]()
                    shiny::req(pp1)

                    dimSel   <- list(staticDimSelections[[idx]])
                    somCodes <- metaD[[somCodesName]]
                    pair     <- dList[[idx]]
                    pctl     <- input$scatterPercentile
                    if (is.null(pctl)) pctl <- 0.99
                    tailP <- (1 - pctl) / 2

                    xlimP <- if (pair[1L] %in% colnames(somCodes))
                        stats::quantile(somCodes[, pair[1L]],
                                        probs = c(tailP, 1 - tailP), na.rm = TRUE)
                    else NULL
                    ylimP <- if (pair[2L] %in% colnames(somCodes))
                        stats::quantile(somCodes[, pair[2L]],
                                        probs = c(tailP, 1 - tailP), na.rm = TRUE)
                    else NULL

                    ggsomPlot(pp1, 1L, rs, dimSel,
                              sce   = sce,
                              metaD = metaD,
                              xlim  = xlimP,
                              ylim  = ylimP)
                })

                ## Suspend while the box is collapsed — no wasted renders.
                shiny::outputOptions(
                    output, paste0("staticSom", idx),
                    suspendWhenHidden = TRUE
                )
            })

        # observe clusterNameSelect ----
        # print clusters based on selection
        shiny::observeEvent(input$clusterNameSelect, {
            outputList <- rv$outputList
            listNames <- names(outputList) %in% input$clusterNameSelect
            combinedSoms <- outputList[listNames] %>%
                unlist() %>%
                unique()
            updateTextInput(session, "clusterNumbers", value = paste(combinedSoms, collapse = ", "))
        })

        # observeEvent groupsVar ----
        shiny::observeEvent(input$groupsVar, {
            # browser()
            groupsVar <- input$groupsVar
            if (!input$groupsVar %in% colnames(metaD$experiment_info)) {
                return(NULL)
            }
            levs <- levels(metaD$experiment_info[, input$groupsVar])
            updateSelectInput(session = session, inputId = "group1", choices = levs)
            updateSelectInput(session = session, inputId = "group2", choices = levs)
        })

        # observe group ----
        shiny::observe({
            if (!input$groupsVar %in% colnames(metaD$experiment_info)) {
                return(NULL)
            }
            levs <- levels(metaD$experiment_info[, input$groupsVar])
            grp1 <- input$group1
            grp2 <- isolate(input$group2)
            levs <- levs[!levs %in% grp1]
            updateSelectInput(session = session, inputId = "group2", choices = levs, selected = grp2)
            # save(file = "dev2.RData", list = c("groupsVar", "levs", "grp1"))
            #
        })
        shiny::observe({
            if (!input$groupsVar %in% colnames(metaD$experiment_info)) {
                return(NULL)
            }
            levs <- levels(metaD$experiment_info[, input$groupsVar])
            grp2 <- input$group2
            grp1 <- isolate(input$group1)
            levs <- levs[!levs %in% grp2]
            updateSelectInput(session = session, inputId = "group1", choices = levs, selected = grp1)
            # save(file = "dev2.RData", list = c("groupsVar", "levs", "grp1"))
            #
        })
        # delayed group input for t-test ----
        # returns sample_ids from potentially other metadata factorials
        groupsInput <- reactive({
            if (!input$groupsVar %in% colnames(metaD$experiment_info)) {
                return(NULL)
            }
            grp1 <- metaD$experiment_info[
                metaD$experiment_info[, input$groupsVar] %in% input$group1, "sample_id"
            ]
            grp2 <- metaD$experiment_info[
                metaD$experiment_info[, input$groupsVar] %in% input$group2, "sample_id"
            ]
            list(group1 = grp1, group2 = grp2)
        }) %>% debounce(1000)

        # t.test ----
        shiny::observe({
            empty <- FALSE
            groupsVar <- input$groupsVar
            if (!input$groupsVar %in% colnames(metaD$experiment_info)) empty <- TRUE
            grpInp <- groupsInput()
            if (is_empty(grpInp)) empty <- TRUE
            if (any(lapply(grpInp, is_empty) %>% unlist())) empty <- TRUE
            rs <- rsUsed()
            req(rs)
            relativeToCol <- input$relativeTo
            numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
            expInfo <- metaD$experiment_info[, numCols, drop = FALSE]
            rownames(expInfo) <- metaD$experiment_info$sample_id
            outputList <- rv$outputList

            if (length(rs) < 1) empty <- TRUE
            if (empty) {
                output$ttestResult <- shiny::renderPrint("no data")
                return(NULL)
            }
            if (relativeToCol == "none") {
                rSums <- rep(1, nrow(clusterPatientTable))
            } else {
                rSums <- .computeRelativeCounts(
                    clusterPatientTable, rs, relativeToCol, expInfo, outputList
                ) / 100
                rSums <- rSums[rownames(clusterPatientTable)]
            }
            names(rSums) <- rownames(clusterPatientTable)
            cD <- SingleCellExperiment::colData(sce)
            cD <- cD[cD$cluster_id %in% rs, ]
            cellCounts <- cD %>%
                as_tibble() %>%
                group_by(sample_id) %>%
                count()
            x <- cellCounts %>%
                filter(sample_id %in% grpInp$group1) %>%
                ungroup()
            x$rsums <- rSums[x$sample_id]
            x$val <- x$n / x$rsums
            x <- x %>% pull(val)
            y <- cellCounts %>%
                filter(sample_id %in% grpInp$group2) %>%
                ungroup()
            y$rsums <- rSums[y$sample_id]
            y$val <- y$n / y$rsums
            y <- y %>% pull(val)
            shiny::req(x)
            shiny::req(y)
            if (sd(x) == 0 & sd(y) == 0) {
                output$ttestResult <- shiny::renderPrint("not enough data")
                return(NULL)
            }
            tt <- stats::t.test(x, y)
            output$ttestResult <- shiny::renderPrint(tt)
        })

        # updatedoutputList function ----
        updatedoutputList <- function() {
            .updateOutputListInputs(session, input, rv$outputList, metaD)
        }

        # observe applyclusterNumbers ---
        shiny::observeEvent(input$applyclusterNumbers, {
            outputList <- rv$outputList
            # isolate(input$clusterNumbers)
            updatedoutputList()
        })

        # obs rmGrp ----
        shiny::observeEvent(input$rmGrp, {
            outputList <- rv$outputList
            cName <- input$clusterNameRM
            cList <- inputClusterNumber()
            if (cName %in% names(outputList)) {
                outputList[[cName]] <- NULL
                outputList[["Rest"]] <- c()
                used <- unique(unlist(outputList))
                all_levels <- levels(sce$cluster_id)
                outputList[["Rest"]] <- as.integer(all_levels[!all_levels %in% used])
                for (na in names(outputList)) {
                    if (length(outputList[[na]]) == 0) outputList[[na]] <- NULL
                }
                # outputList = append(outputList, list(cList ))
                # names(outputList)[length(outputList)] = cName
                rv$outputList <- outputList
                updatedoutputList()
            }
        })


        # obs applyName ----
        shiny::observeEvent(input$applyName, {
            outputList <- rv$outputList
            cName <- input$clusterName
            cList <- inputClusterNumber()
            outputList <- rv$outputList
            outputList[[cName]] <- cList
            outputList[["Rest"]] <- c()
            used <- unique(unlist(outputList))
            all_levels <- levels(sce$cluster_id)
            outputList[["Rest"]] <- as.integer(all_levels[!all_levels %in% used])
            for (na in names(outputList)) {
                if (length(outputList[[na]]) == 0) outputList[[na]] <- NULL
            }
            # browser()
            # outputList = append(outputList, list(cList ))
            # names(outputList)[length(outputList)] = cName
            rv$outputList <- outputList
            updatedoutputList()
            currentSelection <- input$colorbyGroups
            updateSelectInput(
                inputId = "colorbyGroups",
                choices = names(outputList),
                selected = c(currentSelection, cName)
            )
        })

        shiny::observeEvent(input$rmGroups, {
            outputList <- rv$outputList
            rmGroups <- input$groupRM
            rmCluster <- outputList[rmGroups] %>%
                unlist() %>%
                unique()
            rs <- isolate(rsUsed())
            if (is.null(rs)) return(NULL)
            rsUsed(setdiff(rs, rmCluster))
            # keep outputList in sync with the current selection
            rv$outputList <- .rebuildOutputList(rv$outputList, levels(sce$cluster_id), removed = NULL)
            updatedoutputList()
        })

        # output dend ----
        output$dend <- renderPlot({
            dendPlot() %>%
                plot(main = "dendrogram")
        })

        dendPlot <- reactive({
            rs <- rsUsed()
            req(rs)
            selectedPoints()
            labCol <- rep("blue", length(labels(dend)))
            labCol[which(labels(dend) %in% (rs))] <- "red"
            dend %>%
                dendextend::set("leaves_col", labCol) %>%
                dendextend::set("leaves_pch", 15) %>%
                dendextend::set("leaves_cex", 1) %>%
                dendextend::set("labels_cex", 0.5)
        })

        # interactive dendrogram ----
        dendPlotlyData <- shiny::reactive({
            g      <- dendextend::as.ggdend(dend)
            labels <- g$labels

            ## Map character leaf labels → integer cluster ids via dendTable.
            ## Falls back to direct as.integer() when dendTable is not available or
            ## when the labels are already pure integer strings.
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

        output$dendPlotly <- plotly::renderPlotly({
            dd <- dendPlotlyData()
            rs <- rsUsed()
            shiny::req(dd, rs)

            lab_df  <- dd$labels
            seg_df  <- dd$segments
            rs_int  <- as.integer(rs)

            ## Colour leaf nodes: red = selected, black = unselected
            node_col <- ifelse(
                !is.na(lab_df$id) & lab_df$id %in% rs_int,
                "red", "black"
            )

            ## Extend y slightly below 0 so leaf points are never clipped.
            y_min <- min(c(seg_df$y, seg_df$yend, lab_df$y), na.rm = TRUE)
            y_max <- max(c(seg_df$y, seg_df$yend, lab_df$y), na.rm = TRUE)
            y_pad <- (y_max - y_min) * 0.08

            plotly::plot_ly(source = "dendPlotly") %>%
                ## Trace 0: branch segments (not selectable, no customdata)
                plotly::add_segments(
                    data       = seg_df,
                    x          = ~x,    y    = ~y,
                    xend       = ~xend, yend = ~yend,
                    line       = list(color = "black", width = 1L),
                    hoverinfo  = "none",
                    showlegend = FALSE
                ) %>%
                ## Trace 1: leaf points — customdata holds the integer cluster id
                plotly::add_markers(
                    data       = lab_df,
                    x          = ~x,
                    y          = ~y,
                    customdata = ~id,          # integer id, not the character label
                    marker     = list(
                        color = node_col,
                        size  = 10L,
                        line  = list(color = "white", width = 1L)
                    ),
                    hoverinfo  = "none",
                    showlegend = FALSE
                ) %>%
                ## Trace 2: leaf text labels
                plotly::add_text(
                    data         = lab_df,
                    x            = ~x,
                    y            = ~y,
                    text         = ~label,
                    textposition = "bottom center",
                    textfont     = list(size = 9L, color = node_col),
                    hoverinfo    = "none",
                    showlegend   = FALSE
                ) %>%
                plotly::layout(
                    dragmode     = "select",
                    showlegend   = FALSE,
                    xaxis = list(
                        visible  = FALSE,
                        showgrid = FALSE,
                        zeroline = FALSE
                    ),
                    yaxis = list(
                        visible  = FALSE,
                        showgrid = FALSE,
                        zeroline = FALSE,
                        range    = c(y_min - y_pad, y_max + y_pad)
                    ),
                    plot_bgcolor  = "rgba(0,0,0,0)",
                    paper_bgcolor = "rgba(0,0,0,0)",
                    margin        = list(l = 5, r = 5, t = 5, b = 5)
                ) %>%
                plotly::event_register("plotly_selected") %>%
                plotly::event_register("plotly_click")
        })

        shiny::observeEvent(
            safe_event_data(verbose = verbose, "plotly_selected", source = "dendPlotly"),
            {
                d <- safe_event_data(verbose = verbose, "plotly_selected", source = "dendPlotly")
                if (is.null(d) || nrow(d) == 0L) return(NULL)

                ## Keep only the leaf-point trace (curveNumber == 1)
                if ("curveNumber" %in% names(d)) {
                    d <- d[d$curveNumber == 1L, , drop = FALSE]
                }
                if (nrow(d) == 0L || is.null(d$customdata)) return(NULL)

                clicked <- as.integer(unlist(d$customdata))
                clicked <- clicked[!is.na(clicked) & clicked > 0L]
                if (length(clicked) == 0L) return(NULL)

                rs   <- shiny::isolate(rsUsed_d())
                mode <- shiny::isolate(input$selectMode)
                new_rs <- switch(
                    EXPR = mode,
                    "remove others" = clicked,
                    "add"           = unique(c(rs, clicked)),
                    "remove"        = setdiff(rs, clicked),
                    clicked
                )
                shiny::isolate(rsUsed(new_rs))
            }
        )

        shiny::observeEvent(
            safe_event_data(verbose = verbose, "plotly_click", source = "dendPlotly"),
            {
                d <- safe_event_data(verbose = verbose, "plotly_click", source = "dendPlotly")
                if (is.null(d) || nrow(d) == 0L) return(NULL)

                ## Keep only the leaf-point trace (curveNumber == 1)
                if ("curveNumber" %in% names(d)) {
                    d <- d[d$curveNumber == 1L, , drop = FALSE]
                }
                if (nrow(d) == 0L || is.null(d$customdata)) return(NULL)

                clicked <- as.integer(unlist(d$customdata))
                clicked <- clicked[!is.na(clicked) & clicked > 0L]
                if (length(clicked) == 0L) return(NULL)

                rs   <- shiny::isolate(rsUsed_d())
                mode <- shiny::isolate(input$selectMode)
                new_rs <- switch(
                    EXPR = mode,
                    "remove others" = clicked,
                    "add"           = unique(c(rs, clicked)),
                    "remove"        = setdiff(rs, clicked),
                    clicked
                )
                shiny::isolate(rsUsed(new_rs))
            }
        )


        # inputClusterNumber <- reactive ----
        inputClusterNumber <- reactive({
            str_split(input$clusterNumbers, ",")[[1]] %>% as.integer()
        }) %>% debounce(1000)

        # violinPlotSelection <- reactive ----
        violinPlotSelection <- reactive({
            input$violinSelection
        }) %>% debounce(1000)

        # inputClusterNumber() ----
        # Guard: only write rsUsed when the parsed value is genuinely different.
        # Without the setequal check, reactiveVal fires even for identical values,
        # which would bounce updateTextInput → inputClusterNumber → rsUsed forever.
        shiny::observe({
            ic <- inputClusterNumber()
            ic <- suppressWarnings(as.integer(ic))
            ic <- ic[!is.na(ic)]
            ic <- as.integer(intersect(as.character(ic), colnames(clusterPatientTable)))
            if (!setequal(ic, shiny::isolate(rsUsed()))) {
                shiny::isolate(rsUsed(ic))
            }
        })

        # rsUsed_d
        rsUsed_d <- rsUsed %>% debounce(1500)

        # updateTextInput clusterNumbers ----
        shiny::observe({
            rs <- as.integer(unlist(rsUsed()))
            rs <- sort(rs[!is.na(rs)])
            outputList <- rv$outputList
            updateOL <- FALSE
            if (!"selected" %in% names(outputList)) updateOL <- TRUE
            outputList$selected <- rs
            rv$outputList <- outputList
            if (updateOL) updatedoutputList()
            new_text <- paste(rs, collapse = ", ")
            # Guard: skip updateTextInput when the displayed text already matches
            if (!identical(shiny::isolate(input$clusterNumbers), new_text)) {
                updateTextInput(inputId = "clusterNumbers", value = new_text)
            }
        })

        shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_selected", source = "somGrid"), {
            message("somGrid touched")
            # browser()
            rs <- isolate(rsUsed())
            if (is.null(rs)) return(NULL)
            d <- safe_event_data(verbose = verbose, "plotly_selected", source = "somGrid")
            if (is.null(d)) {
                return(NULL)
            }
            d <- .inputSelect(d, rs, isolate(input$selectMode))
            isolate(rsUsed(d))
        })


        shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_selected", source = "tsne"), {
            message("tsne touched")
            rs <- isolate(rsUsed())
            if (is.null(rs)) return(NULL)
            d <- safe_event_data(verbose = verbose, "plotly_selected", source = "tsne")
            if (is.null(d)) {
                return(NULL)
            }
            # browser()
            # save(file = "/pasteur/appa/scratch/bernd/event.RData", list = c("d", "rs"))
            #
            d <- .inputSelect(d, rs, isolate(input$selectMode))
            isolate(rsUsed(d))
        })
        shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_selected", source = "umap"), {
            message("umap touched")
            rs <- isolate(rsUsed())
            if (is.null(rs)) return(NULL)
            d <- safe_event_data(verbose = verbose, "plotly_selected", source = "umap")
            if (is.null(d)) {
                return(NULL)
            }
            d <- .inputSelect(d, rs, isolate(input$selectMode))
            isolate(rsUsed(d))
        })
        shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_selected", source = "pca"), {
            message("pca touched")
            rs <- isolate(rsUsed())
            if (is.null(rs)) return(NULL)
            d <- safe_event_data(verbose = verbose, "plotly_selected", source = "pca")
            if (is.null(d)) {
                return(NULL)
            }
            d <- .inputSelect(d, rs, isolate(input$selectMode))
            isolate(rsUsed(d))
        })

        # output$selected  ----
        output$selected <- renderPrint({

            rs <- rsUsed()
            req(rs)
            rs
        })

        # cellCounts ----
        output$cellCounts <- DT::renderDT(options = list(
            lengthChange = FALSE,
            scrollX = TRUE
        ), {
            cst <- input$compareStatsTo
            sN <- input$singleNode %>% as.integer()
            rs <- rsUsed()
            req(rs)

            outputList <- rv$outputList

            # browser()
            req(clusterPatientTable)
            rSums <- rowSums(clusterPatientTable[, rs, drop = FALSE])
            names(rSums) <- rownames(clusterPatientTable)
            rSums <- data.table::as.data.table(t(rSums))
            if (!cst == "none") {
                if (cst %in% colnames(metaD$experiment_info)) {
                    # expInfo = metaD$experiment_info
                    numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
                    expInfo <- metaD$experiment_info[, numCols, drop = FALSE]
                    eI <- expInfo[, !colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"), drop = FALSE]
                    eI <- apply(eI, 2, as.numeric) %>% as.data.frame()
                    rownames(eI) <- metaD$experiment_info$sample_id
                    if (is.na(eI) %>% any()) {
                        stop("NAs produced when converting experiment_info to numeric")
                    }
                    # ctStats = eI[,cst]
                    ctStats <- data.table::as.data.table(t(eI[colnames(rSums), cst]))
                    colnames(ctStats) <- colnames(rSums)
                    rSums <- data.table::rbindlist(list(rSums, ctStats))
                    # rbind(rSums, eI[colnames(rSums),cst])
                    rownames(rSums) <- c("selection", cst)
                } else if (cst %in% names(outputList)) {
                    ctStats <- data.table::as.data.table(t(rowSums(clusterPatientTable[, outputList[[cst]], drop = FALSE])))
                    colnames(ctStats) <- rownames(clusterPatientTable)
                    rSums <- data.table::rbindlist(list(rSums, ctStats))
                    # rSums = data.table::rbindlist(list(rSums, data.table::as.data.table(t())))
                    rownames(rSums) <- c("selection", cst)
                } else {
                    # browser()
                }
            }
            rn <- rownames(rSums)
            ctStats <- data.table::as.data.table(t(rowSums(clusterPatientTable[, sN, drop = FALSE])))
            colnames(ctStats) <- rownames(clusterPatientTable)
            rSums <- data.table::rbindlist(list(rSums, ctStats))
            # rbind(rSums, rowSums(clusterPatientTable[,sN,drop=FALSE]))
            rownames(rSums) <- c(rn, paste("SOM node", sN))
            # save(file = "/pasteur/appa/scratch/bernd/Rtest3.Rdata",
            #      list = c("clusterPatientTable", "rs", "cst", "rSums", "sN"))
            # cp =load(file = "/pasteur/appa/scratch/bernd/Rtest3.Rdata")
            rSums
        })

        # countBar ----
        output$CountBar <- renderPlot({
            countBarPlot()
        })

        countBarPlot <- reactive({
            cst <- input$compareStatsTo
            # browser()
            outputList <- rv$outputList
            rs <- rsUsed()
            req(rs)
            groupsInput <- groupsInput()
            req(clusterPatientTable)
            countBarPlotFunc(rs, clusterPatientTable, cst, sce, outputList, groupsInput)
        })

        # cellPercentages ----
        output$cellPercentages <- renderPrint({

            outputList <- rv$outputList
            rs <- rsUsed()
            req(rs)

            req(clusterPatientTable)
            relativeToCol <- input$relativeTo
            expInfo <- metadata(sce)$experiment_info
            rownames(expInfo) <- expInfo$sample_id
            expInfo <- expInfo[, !colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"), drop = FALSE]
            eI <- apply(expInfo, 2, as.numeric)

            rSums <- compute_relative_counts(
                clusterPatientTable, rs, relativeToCol, expInfo, outputList
            )
            names(rSums) <- rownames(clusterPatientTable)
            noquote(formatC(signif(rSums, digits = 2), digits = 2, format = "fg", flag = "#"))
        })

        # PercentBar ----
        output$PercentBar <- renderPlot({
            PercentBarPlot()
        })

        PercentBarPlot <- reactive({
            outputList <- rv$outputList
            rs <- rsUsed()
            req(rs)
            req(clusterPatientTable)
            relativeToCol <- input$relativeTo
            groupsInput <- groupsInput()
            PercentBarPlotFunc(sce, relativeToCol, clusterPatientTable, rs, outputList, group, groupsInput)
        })


        ### zooming ----
        zoomFunc <- function(zoom, plotIdx) {
            dimSelectionInternal <- dimSelection()
            req(dimSelectionInternal)
            if (all(!is.null(zoom), "xaxis.range[0]" %in% names(zoom), na.rm = TRUE)) {
                # browser()
                rezoom <- FALSE
                if (all(c(
                    dimSelectionInternal[[plotIdx]]$xzoom[1] > zoom$`xaxis.range[0]`,
                    !is.null(dimSelectionInternal[[plotIdx]]$xzoom[1])
                ), na.rm = TRUE)) {
                    rezoom <- TRUE
                }
                if (rezoom) {

                    # browser()
                    dimSelectionInternal[[plotIdx]]$xzoom <- c(NULL, NULL)
                    dimSelectionInternal[[plotIdx]]$yzoom <- c(NULL, NULL)
                } else {
                    # browser()
                    dimSelectionInternal[[plotIdx]]$xzoom <- c(zoom$`xaxis.range[0]`, zoom$`xaxis.range[1]`)
                    dimSelectionInternal[[plotIdx]]$yzoom <- c(zoom$`yaxis.range[0]`, zoom$`yaxis.range[1]`)
                }
                dimSelection(dimSelectionInternal)
                isolate(triggerRedraw(triggerRedraw() + 1))
            }
        }

        # observe samples2plot. ---
        sample2PlotDb <- reactive(input$samples2plot) %>% debounce(500)

        shiny::observe({
            sampleIds <- sample2PlotDb()
            selected <- sce_subsampled[, sce_subsampled$sample_id %in% sampleIds]
            dfPlot(data.frame(t(assays(selected)[[1]])))
        })

        shiny::observeEvent(dimSelection(), {
            dimSelection <- dimSelection()
            for (idx in seq_along(dimSelection)) {
            }
        })

        shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_selected", source = "scatterPlot"), {

            sc <- sce
            rs <- isolate(rsUsed_d())
            dimSelection <- dimSelection()
            sampleIds <- isolate(input$samples2plot)
            # browser()
            # req(rs)
            d <- safe_event_data(verbose = verbose, "plotly_selected", source = "scatterPlot")
            plotIdx <- activePlot()
            if (is.null(d)) {
                return(NULL)
            }
            req(d$curveNumber)
            d <- d[d$curveNumber == 1, ]
            # d = d$pointNumber+1
            # "view", , "add"
            minx <- min(d$x)
            maxx <- max(d$x)
            miny <- min(d$y)
            maxy <- max(d$y)
            ids <- which(dfPlot()[, dimSelection[[plotIdx]]$dims[1] %>% make.names()] > minx &
                             dfPlot()[, dimSelection[[plotIdx]]$dims[1] %>% make.names()] < maxx &
                             dfPlot()[, dimSelection[[plotIdx]]$dims[2] %>% make.names()] > miny &
                             dfPlot()[, dimSelection[[plotIdx]]$dims[2] %>% make.names()] < maxy)
            ids <- colData(sce_subsampled[, sce_subsampled$sample_id %in% sampleIds])[ids, "cluster_id"]
            # activePlot(1)
            # rs = rsUsed()
            # req(rs)
            # d <- safe_event_data(verbose = verbose, "plotly_selected", source = "somData1")
            # if(is.null(d)){return(NULL)}
            ids <- switch(EXPR = isolate(input$selectMode),
                          "remove others" = intersect(ids, rs),
                          "add" = unique(c(ids, rs)),
                          "remove" = rs[!rs %in% ids],
                          ids
            )
            isolate(rsUsed(ids))
        })



        dimRedSelection <- reactive({
            retVal <- input$dimRedSelection
            retVal
        }) %>% debounce(1000)

        # tsne ----
        tsne <- reactive({
            # browser()
            dimRedCols <- dimRedSelection()
            perplexity <- input$perplexity
            tsne <- tsneFunc(dimRedSelection = dimRedCols, perplexity = perplexity, sce, somCodesName)
            return(tsne)
        }) %>%
            bindCache(dimRedSelection(), input$perplexity) %>%
            debounce(1000)

        umap <- reactive({
            dimRedCols <- input$dimRedSelection
            pumap <- umap::umap.defaults
            pumap$n_neighbors <- input$n_neighbors
            um <- umap::umap(metaD[[somCodesName]][, dimRedCols], config = pumap)
            return(um)
        }) %>%
            bindCache(input$dimRedSelection, input$n_neighbors) %>%
            debounce(1000)

        pca <- reactive({
            dimRedCols <- input$dimRedSelection
            pca <- prcomp(t(metaD[[somCodesName]][, dimRedCols]), scale = FALSE, rank. = 2)
            return(pca)
        }) %>%
            bindCache(input$dimRedSelection) %>%
            debounce(1000)

        output$tsne <- renderPlotly({
            p3 <- tsnePlot()
            showLegend <- input$showlegend
            # browser()
            retVal <- quiet_ggplotly(p3, source = paste0("tsne"), tooltip = "text")

            if (showLegend) {
                retVal <- retVal %>% plotly::layout(legend = list(
                    x = 0, y = -3,
                    xanchor = "left",
                    yanchor = "bottom",
                    orientation = "h"
                ))
            } else {
                retVal <- retVal %>% layout(
                    showlegend = FALSE
                )
            }

            retVal <- retVal %>%
                layout(dragmode = "select") %>%
                event_register("plotly_selected") %>%
                event_register("plotly_relayout")
            retVal
        })

        tsnePlot <- reactive({
            selectedUpdate2()
            rs <- rsUsed_d()
            # browser()
            req(rs)
            triggerRedraw()
            tsne <- tsne()
            # dimSelection =  dimSelection()

            df <- as.data.frame(tsne$Y)
            names(df) <- c("tsne1", "tsne2")
            p3 <- drawProjection(df, rs, colorbyGroups = input$colorbyGroups, sce = sce, outputList = outputList)
            return(p3)
        })

        output$umap <- renderPlotly({
            p3 <- umapPlot()
            showLegend <- input$showlegend
            retVal <- quiet_ggplotly(p3, source = paste0("umap"), tooltip = "text")
            if (showLegend) {
                retVal <- retVal %>% plotly::layout(legend = list(
                    x = 0, y = -3,
                    xanchor = "left",
                    yanchor = "bottom",
                    orientation = "h"
                ))
            } else {
                retVal <- retVal %>% layout(
                    showlegend = FALSE
                )
            }

            retVal <- retVal %>%
                layout(dragmode = "select") %>%
                event_register("plotly_selected") %>%
                event_register("plotly_relayout")
            retVal
        })


        umapPlot <- reactive({
            umap <- umap()
            # plot(um$layout)
            selectedUpdate2()
            rs <- rsUsed_d()
            req(rs)
            triggerRedraw()
            df <- as.data.frame(umap$layout)
            names(df) <- c("umap1", "umap2")

            p3 <- drawProjection(df, rs, colorbyGroups = input$colorbyGroups, sce = sce, outputList = outputList)
            return(p3)
        })

        output$pca <- renderPlotly({
            pca <- pca()
            rs <- rsUsed_d()
            req(rs)
            triggerRedraw()
            req(pca)
            df <- as.data.frame(pca$rotation)
            colnames(df) <- c("pc1", "pc2")
            showLegend <- input$showlegend

            p3 <- pcaPlot()
            retVal <- quiet_ggplotly(p3, source = paste0("pca"), tooltip = "text")

            if (showLegend) {
                retVal <- retVal %>% plotly::layout(legend = list(
                    x = 0, y = -3,
                    xanchor = "left",
                    yanchor = "bottom",
                    orientation = "h"
                ))
            } else {
                retVal <- retVal %>% layout(
                    showlegend = FALSE
                )
            }
            retVal <- retVal %>%
                layout(dragmode = "select") %>%
                event_register("plotly_selected") %>%
                event_register("plotly_relayout")
            return(retVal)
            # retVal %>%
            #   add_trace(x=df$pc1, y=df$pc2, text=~paste("cluster: ", df$cluster),
            #             hoverinfo = 'text', mode = "none")
        })

        pcaPlot <- reactive({
            pca <- pca()
            selectedUpdate2()
            rs <- rsUsed_d()
            req(rs)
            triggerRedraw()
            req(pca)
            df <- as.data.frame(pca$rotation)
            colnames(df) <- c("pc1", "pc2")
            p3 <- drawProjection(df, rs, colorbyGroups = input$colorbyGroups, sce = sce, outputList = outputList)

            return(p3)
        })


        output$somClusters <- renderText({

            rs <- rsUsed()
            req(rs)
            rs <- as.integer(intersect(rs, colnames(clusterPatientTable)))
            # the row it should be plotted
            as.integer(rs) / 40 + 1
            rowElement <- list()
            lapply(40:1, function(x) {
                rowElement[[x]] <- paste(rs[as.integer(rs / 40) + 1 == x], collapse = ", ")
            }) %>%
                unlist() %>%
                paste(collapse = "\n")
            # paste(sort(rs,decreasing = TRUE), collapse = ", ")
        })


        ### output$scatter ----
        output$scatter <- renderPlotly({
            rs <- rsUsed()
            req(rs)
            dimSelection <- dimSelection()
            sampleIds <- input$samples2plot
            # browser()
            plotIdx <- activePlot()
            cidIdx <- colData(sce_subsampled)$cluster_id %in% rs & colData(sce_subsampled)$sample_id %in% sampleIds
            if (length(cidIdx) < 1) {
                return(NULL)
            }
            if (!any(cidIdx)) return(plotly::plotly_empty())

            pp <- scatterPlot()
            if (is.null(pp)) return(plotly::plotly_empty())

            # here we add a grid for selecting points, spanning the current percentile zoom
            chs <- dimSelection[[plotIdx]]$dims
            chx <- make.names(chs[1])
            chy <- make.names(chs[2])
            xVals <- df[cidIdx, chx, drop = TRUE]
            yVals <- df[cidIdx, chy, drop = TRUE]
            pctl <- input$scatterPercentile
            if (is.null(pctl)) pctl <- 0.99
            tailP <- (1 - pctl) / 2
            xlimP <- quantile(xVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
            ylimP <- quantile(yVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
            xpoints <- seq(from = xlimP[1], to = xlimP[2], length.out = 100) %>%
                rep(100) %>%
                sort()
            ypoints <- seq(from = ylimP[1], to = ylimP[2], length.out = 100) %>% rep(100)
            # browser()
            pp <- quiet_ggplotly(pp, source = "scatterPlot") %>%
                add_trace(
                    x = xpoints,
                    y = ypoints,
                    mode = "markers",
                    type = "scatter",
                    fill = "none",
                    fillcolor = "#e763fa",
                    opacity = 0.01 # ,     # size=1.9,
                ) %>%
                layout(dragmode = "select") %>%
                event_register("plotly_selected") %>%
                event_register("plotly_relayout")
            return(pp)
        })

        scatterPlot <- reactive({

            rs <- rsUsed()
            req(rs)
            dimSelection <- dimSelection()
            sampleIds <- input$samples2plot
            req(sampleIds)
            # browser()
            plotIdx <- activePlot()
            req(dimSelection, length(dimSelection) >= plotIdx)
            cidIdx <- colData(sce_subsampled)$cluster_id %in% rs & colData(sce_subsampled)$sample_id %in% sampleIds
            if (sum(cidIdx) < 1) {
                return(NULL)
            }

            chs <- dimSelection[[plotIdx]]$dims
            chx <- make.names(chs[1])
            chy <- make.names(chs[2])
            xVals <- df[cidIdx, chx, drop = TRUE]
            yVals <- df[cidIdx, chy, drop = TRUE]

            pctl <- input$scatterPercentile
            if (is.null(pctl)) pctl <- 0.99
            tailP <- (1 - pctl) / 2
            xlimP <- quantile(xVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
            ylimP <- quantile(yVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)

            pp <- plotCytoScatter(
                rownms = sce_subsampledRN,
                x = sce_subsampled[, cidIdx],
                chs = chs,
                bins = 200
            ) +
                xlim(xlimP) +
                ylim(ylimP)

            return(pp)
        })


        # somRaster ----
        output$somRaster <- renderPlot({
            xy <- somRasterPlot()
            req(xy)
            req(baseRasterGgplot)
            baseRasterGgplot +
                geom_point(data = xy, aes(x = x, y = y), color = "red", size = 1, inherit.aes = FALSE) +
                geom_point(data = xy, aes(x = x, y = y), color = "red", shape = 3, size = 1, inherit.aes = FALSE)
        })

        # somRasterSelect ----
        output$somRasterSelect <- renderPlotly({
            rs <- rsUsed()
            if (is.null(rs)) return(NULL)
            p3 <- .buildSOMRasterSelectPlot(somRasterData, rs)
            req(inherits(p3, "ggplot"))
            quiet_ggplotly(p3, source = "somGrid", tooltip = "") %>%
                plotly::layout(showlegend = FALSE, dragmode = "select") %>%
                plotly::event_register("plotly_selected") %>%
                plotly::event_register("plotly_relayout")
        })

        somRasterPlot <- reactive({
            rs <- as.integer(unlist(rsUsed()))
            rs <- rs[!is.na(rs) & rs > 0L]
            req(length(rs) > 0L)
            # id-based subsetting: safe regardless of row ordering
            somRasterData[somRasterData$id %in% rs, c("x", "y"), drop = FALSE]
        })

        # flowSOMPie ----
        output$flowSOMPie <- renderPlot({
            p <- flowSOMPiePlot()
            req(p)
            p
        })

        flowSOMPiePlot <- reactive({
            rs <- rsUsed()
            req(rs)
            ncol_val <- input$flowSOMPieCols          # NULL-safe: returns NULL before widget exists
            .buildFlowSOMPiePlot(metaD[[somCodesName]], rs, colsUsed, ncol = ncol_val)
        })

        # flowSOMStars ----
        if (!is.null(fsom)) {
            output$flowSOMStars <- renderPlotly({
                req(fsom)
                ret <- FlowSOM::PlotStars(
                    fsom,
                    backgroundValues = fsom$metaclustering,
                    list_insteadof_ggarrange = TRUE
                )
                p <- ret$tree
                coords <- fsom$MST$l
                coord_df <- data.frame(
                    x = coords[, 1],
                    y = coords[, 2],
                    id = seq_len(nrow(coords))
                )
                p_click <- p + geom_point(
                    data = coord_df,
                    aes(x = x, y = y, customdata = id),
                    alpha = 0,
                    size = 8,
                    inherit.aes = FALSE
                )
                quiet_ggplotly(p_click, source = "flowSOMStars", tooltip = "") %>%
                    layout(dragmode = "select") %>%
                    event_register("plotly_selected") %>%
                    event_register("plotly_click")
            })

            shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_selected", source = "flowSOMStars"), {
                d <- safe_event_data(verbose = verbose, "plotly_selected", source = "flowSOMStars")
                if (is.null(d) || nrow(d) == 0 || is.null(d$customdata)) {
                    return(NULL)
                }
                clicked <- as.integer(d$customdata)
                rs <- isolate(rsUsed_d())
                mode <- isolate(input$selectMode)
                new_rs <- switch(EXPR = mode,
                                 "remove others" = clicked,
                                 "add" = unique(c(rs, clicked)),
                                 "remove" = setdiff(rs, clicked),
                                 clicked
                )
                isolate(rsUsed(new_rs))
            })

            shiny::observeEvent(safe_event_data(verbose = verbose, "plotly_click", source = "flowSOMStars"), {
                d <- safe_event_data(verbose = verbose, "plotly_click", source = "flowSOMStars")
                if (is.null(d) || is.null(d$customdata)) {
                    return(NULL)
                }
                clicked <- as.integer(d$customdata)
                rs <- isolate(rsUsed_d())
                mode <- isolate(input$selectMode)
                new_rs <- switch(EXPR = mode,
                                 "remove others" = clicked,
                                 "add" = unique(c(rs, clicked)),
                                 "remove" = setdiff(rs, clicked),
                                 clicked
                )
                isolate(rsUsed(new_rs))
            })
        }

        ## VlnPlot ----

        output$VlnPlot <- renderPlot({
            vlnPlot()
        })
        shiny::observe({
        })
        vlnPlot <- reactive({
            # req(input$violinBox)  # only compute while violin box is expanded
            # observe
            # this changes the groups
            input$applyName
            input$rmGrp
            violinSelection <- violinPlotSelection()
            upsetSelection <- input$upsetSelection
            outputList <- rv$outputList

            req(outputList)

            # save(file = "vln.Rdata", list=ls())
            ##                  Cell Idents  Feat     Expr
            # Idents = named groups of clusters
            # Feat = marker
            # Expr = SOM value
            # cell = SOM node
            markers <- colnames(sce@metadata[[somCodesName]])
            if (length(upsetSelection) < 3) {
                upsetSelection <- names(outputList)
            }
            markers <- intersect(markers, violinSelection)
            somCodes <- sce@metadata[[somCodesName]]
            nNodes <- nrow(somCodes)
            parts <- lapply(upsetSelection, function(na) {
                rows <- outputList[[na]]
                rows <- rows[rows != 0]
                if (length(rows) < 2) {
                    return(NULL)
                }
                wide <- as.data.frame(somCodes[rows, markers, drop = FALSE])
                wide$somNode <- factor(rows, levels = seq_len(nNodes))
                long <- tidyr::pivot_longer(wide,
                                            cols = markers,
                                            names_to = "marker", values_to = "expr"
                )
                long$grpName <- na
                long
            })
            data <- data.table::rbindlist(parts)
            if (nrow(data) == 0) {
                return(
                    ggplot() + theme_void() +
                        labs(title = "No groups with \u2265 2 SOM nodes")
                )
            }

            data$marker <- factor(data$marker, levels = markers)
            data$grpName <- factor(data$grpName, levels = upsetSelection)

            nb.cols <- length(unique(data$marker))
            mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

            p <- ggplot(data, aes(factor(marker), expr, fill = marker)) +
                geom_violin(scale = "width", adjust = 1, trim = TRUE) +
                scale_y_continuous(expand = c(0, 0), position = "right", labels = function(x) {
                    c(rep(x = "", times = length(x) - 2), x[length(x) - 1], "")
                }) +
                facet_grid(rows = vars(grpName), scales = "free", switch = "y") +
                theme_cowplot(font_size = 12) +
                theme(
                    legend.position = "none", panel.spacing = unit(0, "lines"),
                    plot.title = element_text(hjust = 0.5),
                    panel.background = element_rect(fill = NA, color = "black"),
                    strip.background = element_blank(),
                    strip.text = element_text(face = "bold"),
                    strip.text.y.left = element_text(angle = 0),
                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
                ) +
                scale_fill_manual(values = mycolors) +
                ggtitle("Marker on x-axis") +
                xlab("Marker") +
                ylab("Expression Level")

            p
        })
        ## VlnPlot2 ----
        output$VlnPlot2 <- renderPlot({
            VlnPlot2()
        })


        VlnPlot2 <- reactive({
            # req(input$violinBox)  # only compute while violin box is expanded
            # observe
            # this changes the groups
            input$applyName
            input$rmGrp
            upsetSelection <- input$upsetSelection
            violinSelection <- violinPlotSelection()
            outputList <- rv$outputList
            req(outputList)
            # save(file = "vln.Rdata", list=ls())
            ##                  Cell Idents  Feat     Expr
            # Idents = named groups of clusters
            # Feat = marker
            # Expr = SOM value
            # cell = SOM node
            p <- plotViolin2Func(sce, somCodesName, violinSelection, upsetSelection, outputList)
            return(p)
        })

        ## close app ----
        shiny::observeEvent(input$close, {
            shinyjs::js$closeWindow()
            assign(x = "outputList", value = rv$outputList, envir = env)
            shiny::stopApp(rv$outputList)
        })

        ## upset plot ----
        output$UpSet <- renderPlot({
            upSetPlot()
        })

        selectedUpdate <- reactiveVal(value = 0)

        shiny::observe({
            rs <- rsUsed_d()
            if ("selected" %in% input$upsetSelection) {
                isolate(selectedUpdate(selectedUpdate() + 1))
            }
        })

        upSetPlot <- reactive({
            # req(input$upsetBox)  # only compute while UpSet box is expanded
            # inputList = reactiveValuesToList(input)
            # save(file = "UpSet.RData", list = c(ls()))
            # cp = load("UpSet.RData")
            selectedUpdate()
            input$applyName
            upsetSelection <- input$upsetSelection
            input$rmGrp
            outputList <- rv$outputList
            upsetPlotFunc(upsetSelection, outputList, sce)
        })


        output$downloadPlots <- downloadHandler(
            filename = function() {
                req(input$clusterNameSelect)
                paste(input$clusterNameSelect, ".pdf", sep = "")
            },
            content = function(file) {
                dimSelection <- dimSelection()
                rs <- c(1, 2, 3, 4, 8)
                rs <- rsUsed()
                req(rs)
                pdf(file = file)

                plotIdx <- 1

                message("printing dendPlot", class(dendPlot()))
                dendPlot() %>%
                    plot(main = "dendrogram")
                message("printing countBarPlot")
                print(countBarPlot())

                message("printing PercentBarPlot")
                print(PercentBarPlot())


                pctl <- input$scatterPercentile
                if (is.null(pctl)) pctl <- 0.99
                tailP <- (1 - pctl) / 2
                somCodes <- metaD[[somCodesName]]
                dlColorVar <- input$somColorVar
                dlSizeVar <- input$somSizeVar
                if (is.null(dlColorVar)) dlColorVar <- "n"
                if (is.null(dlSizeVar)) dlSizeVar <- "max"
               {
                    message("printing somScatter (current dims)")
                    dimSel   <- dimSelection()
                    dlDimSel <- if (length(dimSel) >= 1L) dimSel else list(list(dims = activeDims()))
                    dims     <- dlDimSel[[1L]]$dims

                    pctl  <- input$scatterPercentile
                    if (is.null(pctl)) pctl <- 0.99
                    tailP    <- (1 - pctl) / 2
                    somCodes <- metaD[[somCodesName]]

                    dlColorVar <- input$somColorVar;  if (is.null(dlColorVar)) dlColorVar <- "n"
                    dlSizeVar  <- input$somSizeVar;   if (is.null(dlSizeVar))  dlSizeVar  <- "max"

                    pp1   <- plotSOMScatter(
                        x         = sce,
                        chs       = dims,
                        pointSize = dlSizeVar,
                        color_by  = dlColorVar,
                        xRN       = sceRN, xCN = sceCN
                    )
                    xlimP <- if (dims[1L] %in% colnames(somCodes))
                        stats::quantile(somCodes[, dims[1L]], probs = c(tailP, 1 - tailP), na.rm = TRUE)
                    else NULL
                    ylimP <- if (dims[2L] %in% colnames(somCodes))
                        stats::quantile(somCodes[, dims[2L]], probs = c(tailP, 1 - tailP), na.rm = TRUE)
                    else NULL

                    print(ggsomPlot(pp1, 1L, rs, dlDimSel,
                                    sce = sce, metaD = metaD,
                                    xlim = xlimP, ylim = ylimP))
                }
                message("printing tsnePlot")
                print(tsnePlot())
                message("printing umapPlot")
                print(umapPlot())
                message("printing pcaPlot")
                print(pcaPlot())

                message("printing scatterPlot")
                print(scatterPlot())
                message("printing somRasterPlot")
                xy <- somRasterPlot()
                if (!is.null(xy) && !is.null(baseRasterGgplot)) {
                    print(baseRasterGgplot +
                              geom_point(data = xy, aes(x = x, y = y), color = "red", size = 1, inherit.aes = FALSE) +
                              geom_point(data = xy, aes(x = x, y = y), color = "red", shape = 3, size = 1, inherit.aes = FALSE))
                }

                message("printing vlnPlot")
                dlOutputList <- rv$outputList
                dlUpsetSel <- input$upsetSelection
                dlViolinSel <- input$violinSelection
                if (length(dlUpsetSel) < 3) dlUpsetSel <- names(dlOutputList)
                if (is.null(dlViolinSel)) dlViolinSel <- colsUsed
                pVln <- plotViolinFunc(sce, somCodesName, dlUpsetSel, dlOutputList, dlViolinSel)
                if (!is.null(pVln)) print(pVln)

                message("printing VlnPlot2")
                pVln2 <- plotViolin2Func(sce, somCodesName, dlViolinSel, dlUpsetSel, dlOutputList)
                if (!is.null(pVln2)) print(pVln2)

                message("printing upSetPlot")
                pUpset <- upsetPlotFunc(dlUpsetSel, dlOutputList, sce)
                if (!is.null(pUpset)) print(pUpset)
                dev.off()
                message("done printing")
            }
        )


        dimSelection <- shiny::reactiveVal(list())
        shiny::observe({
            dims <- activeDims()
            shiny::req(dims[1L], dims[2L])
            lim1 <- channelLimits[[dims[1L]]]
            lim2 <- channelLimits[[dims[2L]]]
            shiny::req(lim1, lim2)
            # Preserve existing zoom when user only changes color/size, not dims
            existing <- shiny::isolate(dimSelection())
            if (length(existing) >= 1L && identical(existing[[1L]]$dims, dims)) return()
            dimSelection(list(list(
                dims  = dims,
                xlim  = c(lim1["min"], lim1["max"]),
                ylim  = c(lim2["min"], lim2["max"]),
                xzoom = c(NULL, NULL),
                yzoom = c(NULL, NULL)
            )))
        })
    }


    shiny::shinyApp(ui = ui, server = server)
}
