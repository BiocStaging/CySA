# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# outputs_clusterSelector.R ----
# Shiny output renderers for the clusterSelector app. Outputs are registered
# separately from observers to follow Bioconductor's Shiny app packaging
# guidelines.

#' Register clusterSelector Outputs
#'
#' Attaches all `output$*` renderers used by the clusterSelector app. Observer
#' logic lives in `.registerClusterSelectorObservers()`.
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
#' @param colTree Optional collapsible tree object.
#' @param somRasterData Data frame for SOM raster visualization.
#' @param somRasterObj Raster object for SOM visualization.
#' @param baseRasterGgplot Pre-built ggplot raster background.
#' @param somRasterPlot Reactive selected-node coordinates.
#' @param somBasePlot Reactive base SOM scatter plot.
#' @param somProjectionDf Reactive projection data frame for showGroups mode.
#' @param tsnePlot,umapPlot,pcaPlot Reactive dimension-reduction plots.
#' @param scatterPlot Reactive 2-D scatter plot.
#' @param dendPlot Reactive dendrogram object.
#' @param dendPlotlyData Reactive dendrogram plotting data.
#' @param countBarPlot,PercentBarPlot Reactive stats bar plots.
#' @param vlnPlot,VlnPlot2 Reactive violin plots.
#' @param upSetPlot Reactive UpSet plot.
#' @param flowSOMPiePlot Reactive FlowSOM marker pie plot.
#' @param groupsInput Debounced reactive of sample-group selections.
#' @param channelLimits Named list of per-channel axis limits.
#' @param getBasePlot Function that returns a cached base SOM ggplot.
#' @param fsom Optional `FlowSOM` object.
#' @param sce Full `SingleCellExperiment`.
#' @param sce_subsampled Subsampled `SingleCellExperiment`.
#' @param sce_subsampledRN,sceRN,sceCN Row/column name helpers.
#' @param metaD `metadata(sce)` list.
#' @param dend Dendrogram object.
#' @param dendTable Data frame for dendrogram navigation.
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param colsUsed Character vector of SOM columns used.
#' @param rnSCE `rownames(sce)`.
#' @param nPlots Number of 2D SOM plots.
#' @param df Static subsampled data frame.
#' @param dfPlot Reactive subsampled data frame.
#' @param env Environment used to store mutable state.
#' @param verbose Logical indicating whether to print debug messages.
#'
#' @keywords internal
#' @name INTERNAL_registerClusterSelectorOutputs
.registerClusterSelectorOutputs <- function(
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
  colsUsed, rnSCE, nPlots, df, dfPlot, env, verbose
) {
  # Collapsible sample tree.
  if (!is.null(colTree)) {
    output$plot <- collapsibleTree::renderCollapsibleTree({
      colTree
    })
  }

  # Dimension-pair preset picker rendered from the reactive dListRV.
  output$dimPairSelectUI <- shiny::renderUI({
    dList <- dListRV()
    if (length(dList) == 0L) {
      return(shiny::selectInput("dimPairSelect", "Preset pair", choices = character(0L)))
    }
    labels <- vapply(dList, function(p) paste(p[1L], "-", p[2L]), character(1L))
    choices <- stats::setNames(seq_along(dList), labels)
    shiny::selectInput(
      "dimPairSelect",
      "Preset pair",
      choices = choices,
      selected = shiny::isolate(activePlot())
    )
  })

  # Main 2-D scatter plot.
  output$scatter <- plotly::renderPlotly({
    # message("DIAG[scatter] renderPlotly block ENTERED at ", Sys.time())

    rs <- rsUsed()
    # message("DIAG[scatter] rs = ", if (is.null(rs)) "NULL" else if (length(rs) == 0) "empty" else paste(rs, collapse = ","))
    shiny::req(rs)

    dimSel <- dimSelection()
    sampleIds <- input$samples2plot
    # message("DIAG[scatter] dimSel length = ", length(dimSel),
    #         " | dimSel[[1]]$dims = ", if (length(dimSel) >= 1L) paste(dimSel[[1L]]$dims, collapse = ", ") else "none",
    #         " | sampleIds = ", if (is.null(sampleIds)) "NULL" else paste(sampleIds, collapse = ","))
    shiny::req(dimSel, length(dimSel) >= 1L)

    cidIdx <- SingleCellExperiment::colData(sce_subsampled)$cluster_id %in% rs &
      SingleCellExperiment::colData(sce_subsampled)$sample_id %in% sampleIds
    # message("DIAG[scatter] any(cidIdx) = ", any(cidIdx), " | sum(cidIdx) = ", sum(cidIdx))
    if (!any(cidIdx)) {
      return(plotly::plotly_empty())
    }

    pp <- scatterPlot()
    # message("DIAG[scatter] is.null(pp) = ", is.null(pp))
    if (is.null(pp)) {
      return(plotly::plotly_empty())
    }

    chs <- dimSel[[1L]]$dims
    # message("DIAG[scatter] chs (dims used) = ", if (is.null(chs)) "NULL" else paste(chs, collapse = ", "))
    chx <- make.names(chs[1L])
    chy <- make.names(chs[2L])
    xVals <- df[cidIdx, chx, drop = TRUE]
    yVals <- df[cidIdx, chy, drop = TRUE]
    # message("DIAG[scatter] length(xVals) = ", length(xVals), " | length(yVals) = ", length(yVals))

    pctl <- input$scatterPercentile
    if (is.null(pctl)) pctl <- 0.99
    tailP <- (1 - pctl) / 2
    xlimP <- stats::quantile(xVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
    ylimP <- stats::quantile(yVals, probs = c(tailP, 1 - tailP), na.rm = TRUE)
    # message("DIAG[scatter] xlimP = ", paste(xlimP, collapse = ", "), " | ylimP = ", paste(ylimP, collapse = ", "))

    xpoints <- sort(rep(seq(from = xlimP[1L], to = xlimP[2L], length.out = 100L), 100L))
    ypoints <- rep(seq(from = ylimP[1L], to = ylimP[2L], length.out = 100L), 100L)

    quiet_ggplotly(pp, source = "scatterPlot") %>%
      plotly::add_trace(
        x = xpoints, y = ypoints,
        mode = "markers", type = "scatter",
        fill = "none", fillcolor = "#e763fa", opacity = 0.01
      ) %>%
      plotly::layout(dragmode = "select") %>%
      plotly::event_register("plotly_selected") %>%
      plotly::event_register("plotly_relayout")
  })

  # SOM 2D plots.
  lapply(seq_len(nPlots), function(i) {
    local({
      plotIdxLocal <- i

      output[[paste0("somData", plotIdxLocal)]] <- plotly::renderPlotly({
        colorbyGroups <- input$colorbyGroups

        selectedUpdate2()
        showGroups <- input$showGroups

        dimSel <- dimSelection()
        rs <- rsUsed_d()
        shiny::req(rs)
        triggerRedraw()

        plotIdx <- if (plotIdxLocal == nPlots) activePlot() else plotIdxLocal
        shiny::req(length(dimSel) >= plotIdx)
        dims <- dimSel[[plotIdx]]$dims
        colorVar <- input$somColorVar
        sizeVar <- input$somSizeVar
        if (is.null(colorVar)) colorVar <- "n"
        if (is.null(sizeVar)) sizeVar <- "max"

        pp1 <- getBasePlot(dims[1L], dims[2L], colorVar, sizeVar)
        if (showGroups) {
          projectionDf <- somProjectionDf()
          req(projectionDf)
          p3 <- drawProjection(
            projectionDf, rs,
            colorbyGroups = colorbyGroups,
            sce = sce, outputList = rv$outputList
          )
        } else {
          p3 <- ggsomPlot(
            pp1, plotIdx, rs,
            dimSelection = dimSel,
            somCodesName = somCodesName, sce = sce, metaD = metaD
          )
        }
        if (is.null(p3)) {
          return(NULL)
        }

        source_id <- paste0("somData", plotIdxLocal)
        quiet_ggplotly(p3, source = source_id, tooltip = "") %>%
          plotly::layout(showlegend = FALSE, dragmode = "select") %>%
          plotly::event_register("plotly_selected") %>%
          plotly::event_register("plotly_relayout")
      })
    })
  })

  # SOM 2D Main Plot ----
  # Single active plot showing the currently selected dimension pair.
  output$somDataMain <- plotly::renderPlotly({
    colorbyGroups <- input$colorbyGroups
    selectedUpdate2()
    showGroups <- input$showGroups
    dimSel <- dimSelection()
    shiny::req(length(dimSel) >= 1L)
    rs <- rsUsed_d()
    shiny::req(rs)
    triggerRedraw()

    dims <- dimSel[[1L]]$dims
    pp1 <- somBasePlot()
    pDf <- if (isTRUE(showGroups)) somProjectionDf() else NULL

    pctl <- input$scatterPercentile
    if (is.null(pctl)) pctl <- 0.99
    tailP <- (1 - pctl) / 2
    somCodes <- metaD[[somCodesName]]
    xlimP <- if (dims[1L] %in% colnames(somCodes)) {
      stats::quantile(somCodes[, dims[1L]], probs = c(tailP, 1 - tailP), na.rm = TRUE)
    } else {
      NULL
    }
    ylimP <- if (dims[2L] %in% colnames(somCodes)) {
      stats::quantile(somCodes[, dims[2L]], probs = c(tailP, 1 - tailP), na.rm = TRUE)
    } else {
      NULL
    }

    somPlot(
      pp1, 1L, rs, colorbyGroups, showGroups,
      dimSelection = dimSel,
      sce = sce,
      metaD = metaD,
      outputList = rv$outputList,
      projectionDf = pDf,
      xlim = xlimP,
      ylim = ylimP,
      source = "somDataMain"
    )
  })

  # Static SOM Grid UI ----
  output$staticSomGrid <- shiny::renderUI({
    current <- dListRV()
    n_cols <- as.integer(input$staticSomCols)
    if (is.na(n_cols) || n_cols < 1L) n_cols <- 3L
    height <- paste0(
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

    row_starts <- seq(1L, length(panels), by = n_cols)
    rows <- lapply(row_starts, function(s) {
      do.call(
        shiny::fluidRow,
        panels[s:min(s + n_cols - 1L, length(panels))]
      )
    })
    do.call(htmltools::tagList, rows)
  })

  # Static SOM Dynamic Renderers ----
  shiny::observe({
    current <- dListRV()

    lapply(seq_along(current), function(idx) {
      local({
        i <- idx
        pair_i <- current[[i]]
        out_id <- paste0("staticSomDyn", i)

        output[[out_id]] <- shiny::renderPlot({
          rs <- rsUsed_d()
          shiny::req(rs)
          colorVar <- if (is.null(input$somColorVar)) "n" else input$somColorVar
          sizeVar <- if (is.null(input$somSizeVar)) "max" else input$somSizeVar
          pctl <- if (is.null(input$scatterPercentile)) 0.99 else input$scatterPercentile

          pp1 <- getBasePlot(pair_i[1L], pair_i[2L], colorVar, sizeVar)
          shiny::req(pp1)

          lim1 <- channelLimits[[pair_i[1L]]]
          lim2 <- channelLimits[[pair_i[2L]]]
          if (is.null(lim1)) lim1 <- c(min = 0, max = 1)
          if (is.null(lim2)) lim2 <- c(min = 0, max = 1)

          dimSel_i <- list(list(
            dims = pair_i,
            xlim = c(lim1["min"], lim1["max"]),
            ylim = c(lim2["min"], lim2["max"]),
            xzoom = c(NULL, NULL),
            yzoom = c(NULL, NULL)
          ))

          tailP <- (1 - pctl) / 2
          somCodes <- metaD[[somCodesName]]
          xlimP <- if (pair_i[1L] %in% colnames(somCodes)) {
            stats::quantile(somCodes[, pair_i[1L]],
              probs = c(tailP, 1 - tailP), na.rm = TRUE
            )
          } else {
            NULL
          }
          ylimP <- if (pair_i[2L] %in% colnames(somCodes)) {
            stats::quantile(somCodes[, pair_i[2L]],
              probs = c(tailP, 1 - tailP), na.rm = TRUE
            )
          } else {
            NULL
          }

          ggsomPlot(pp1, 1L, rs, dimSel_i,
            sce = sce,
            metaD = metaD,
            xlim = xlimP,
            ylim = ylimP
          )
        })

        shiny::outputOptions(output, out_id, suspendWhenHidden = TRUE)
      })
    })
  })

  # t-SNE, UMAP, PCA plots.
  output$tsne <- plotly::renderPlotly({
    p3 <- tsnePlot()
    shiny::req(p3)
    showLegend <- input$showlegend
    ret <- quiet_ggplotly(p3, source = "tsne", tooltip = "text")
    if (showLegend) {
      ret <- ret %>% plotly::layout(legend = list(
        x = 0, y = -0.2, xanchor = "left", yanchor = "bottom", orientation = "h"
      ))
    } else {
      ret <- ret %>% plotly::layout(showlegend = FALSE)
    }
    ret %>%
      plotly::layout(dragmode = "select") %>%
      plotly::event_register("plotly_selected") %>%
      plotly::event_register("plotly_relayout")
  })

  output$umap <- plotly::renderPlotly({
    p3 <- umapPlot()
    shiny::req(p3)
    showLegend <- input$showlegend
    ret <- quiet_ggplotly(p3, source = "umap", tooltip = "text")
    if (showLegend) {
      ret <- ret %>% plotly::layout(legend = list(
        x = 0, y = -0.2, xanchor = "left", yanchor = "bottom", orientation = "h"
      ))
    } else {
      ret <- ret %>% plotly::layout(showlegend = FALSE)
    }
    ret %>%
      plotly::layout(dragmode = "select") %>%
      plotly::event_register("plotly_selected") %>%
      plotly::event_register("plotly_relayout")
  })

  output$pca <- plotly::renderPlotly({
    p3 <- pcaPlot()
    shiny::req(p3)
    showLegend <- input$showlegend
    ret <- quiet_ggplotly(p3, source = "pca", tooltip = "text")
    if (showLegend) {
      ret <- ret %>% plotly::layout(legend = list(
        x = 0, y = -0.2, xanchor = "left", yanchor = "bottom", orientation = "h"
      ))
    } else {
      ret <- ret %>% plotly::layout(showlegend = FALSE)
    }
    ret %>%
      plotly::layout(dragmode = "select") %>%
      plotly::event_register("plotly_selected") %>%
      plotly::event_register("plotly_relayout")
  })

  # Dendrogram outputs.
  output$dend <- shiny::renderPlot({
    dendPlot() %>% plot(main = "dendrogram")
  })

  output$dendPlotly <- plotly::renderPlotly({
    dd <- dendPlotlyData()
    rs <- rsUsed()
    shiny::req(dd, rs)
    label_cols <- ifelse(
      dd$labels$label %in% as.character(rs), "red", "black"
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = dd$segments,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
      ) +
      ggplot2::geom_point(
        data = dd$labels,
        ggplot2::aes(x = x, y = y), # customdata removed
        color = label_cols
      ) +
      ggplot2::geom_text(
        data = dd$labels,
        ggplot2::aes(x = x, y = y, label = label),
        vjust = 1, size = 3
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text  = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )

    quiet_ggplotly(p, source = "dendPlotly", tooltip = "") %>%
      plotly::style( # inject after ggplotly()
        customdata = dd$labels$label,
        traces     = 2L # geom_point is trace 2
      ) %>%
      plotly::layout(dragmode = "select") %>%
      plotly::event_register("plotly_selected") %>%
      plotly::event_register("plotly_click")
  })

  # Stats outputs.
  output$somClusters <- shiny::renderText({
    rs <- rsUsed()
    shiny::req(rs)
    rs <- as.integer(intersect(rs, colnames(clusterPatientTable)))
    rowElement <- list()
    vals <- lapply(seq_len(40), function(x) {
      paste(rs[as.integer(rs / 40) + 1 == x], collapse = ", ")
    })
    paste(unlist(vals), collapse = "\n")
  })

  output$cellCounts <- DT::renderDT(options = list(
    lengthChange = FALSE,
    scrollX = TRUE
  ), {
    cst <- input$compareStatsTo
    sN <- as.integer(input$singleNode)
    rs <- rsUsed()
    shiny::req(rs)
    outputList <- rv$outputList
    shiny::req(clusterPatientTable)

    rSums <- rowSums(clusterPatientTable[, rs, drop = FALSE])
    names(rSums) <- rownames(clusterPatientTable)
    rSums <- data.table::as.data.table(t(rSums))

    if (!identical(cst, "none")) {
      if (cst %in% colnames(metaD$experiment_info)) {
        numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
        expInfo <- metaD$experiment_info[, numCols, drop = FALSE]
        eI <- expInfo[, !colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"), drop = FALSE]
        eI <- apply(eI, 2, as.numeric) %>% as.data.frame()
        rownames(eI) <- metaD$experiment_info$sample_id
        if (any(is.na(eI))) {
          stop("NAs produced when converting experiment_info to numeric")
        }
        ctStats <- data.table::as.data.table(t(eI[colnames(rSums), cst, drop = FALSE]))
        colnames(ctStats) <- colnames(rSums)
        rSums <- data.table::rbindlist(list(rSums, ctStats))
        rownames(rSums) <- c("selection", cst)
      } else if (cst %in% names(outputList)) {
        ctStats <- data.table::as.data.table(
          t(rowSums(clusterPatientTable[, outputList[[cst]], drop = FALSE]))
        )
        colnames(ctStats) <- rownames(clusterPatientTable)
        rSums <- data.table::rbindlist(list(rSums, ctStats))
        rownames(rSums) <- c("selection", cst)
      }
    }
    rn <- rownames(rSums)
    ctStats <- data.table::as.data.table(
      t(rowSums(clusterPatientTable[, sN, drop = FALSE]))
    )
    colnames(ctStats) <- rownames(clusterPatientTable)
    rSums <- data.table::rbindlist(list(rSums, ctStats))
    rownames(rSums) <- c(rn, paste("SOM node", sN))
    rSums
  })


  output$cellPercentages <- shiny::renderPrint({
    outputList <- rv$outputList
    rs <- rsUsed()
    shiny::req(rs)
    shiny::req(clusterPatientTable)
    relativeToCol <- input$relativeTo
    expInfo <- S4Vectors::metadata(sce)$experiment_info
    rownames(expInfo) <- expInfo$sample_id
    numCols <- unlist(lapply(expInfo, is.numeric), use.names = FALSE)
    expInfo <- expInfo[, numCols, drop = FALSE]
    eI <- apply(expInfo, 2, as.numeric)
    rSums <- compute_relative_counts(
      clusterPatientTable, rs, relativeToCol, eI, outputList
    )
    names(rSums) <- rownames(clusterPatientTable)
    noquote(formatC(signif(rSums, digits = 2), digits = 2, format = "fg", flag = "#"))
  })
  output$CountBar <- shiny::renderPlot({
    countBarPlot()
  })

  output$PercentBar <- shiny::renderPlot({
    PercentBarPlot()
  })

  output$ttestResult <- shiny::renderPrint({
      empty <- FALSE
      groupsVar <- input$groupsVar
      if (is.null(groupsVar) || length(groupsVar) == 0 ||
          !groupsVar %in% colnames(metaD$experiment_info)) {
          empty <- TRUE
      }
      grpInp <- groupsInput()
      if (purrr::is_empty(grpInp)) empty <- TRUE
      if (any(unlist(lapply(grpInp, purrr::is_empty)))) empty <- TRUE
      rs <- rsUsed()
      shiny::req(rs)
      relativeToCol <- input$relativeTo
      outputList <- rv$outputList
      if (length(rs) < 1) empty <- TRUE
      if (empty) {
          return("no data")
      }
      if (relativeToCol == "none") {
          rSums <- rep(1, nrow(clusterPatientTable))
      } else {
          numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
          expInfo <- metaD$experiment_info[, numCols, drop = FALSE]
          eI <- apply(expInfo, 2, as.numeric) %>% as.data.frame()
          rownames(eI) <- metaD$experiment_info$sample_id
          rSums <- .computeRelativeCounts(
              clusterPatientTable, rs, relativeToCol, eI, outputList
          ) / 100
          rSums <- rSums[rownames(clusterPatientTable)]
      }
      names(rSums) <- rownames(clusterPatientTable)
      cD <- SingleCellExperiment::colData(sce)
      cD <- cD[cD$cluster_id %in% rs, ]
      cellCounts <- cD %>%
          dplyr::as_tibble() %>%
          dplyr::group_by(.data$sample_id) %>%
          dplyr::count()
      x <- cellCounts %>%
          dplyr::filter(.data$sample_id %in% grpInp$group1) %>%
          dplyr::ungroup()
      x$rsums <- rSums[x$sample_id]
      x$val <- x$n / x$rsums
      x <- x$val
      y <- cellCounts %>%
          dplyr::filter(.data$sample_id %in% grpInp$group2) %>%
          dplyr::ungroup()
      y$rsums <- rSums[y$sample_id]
      y$val <- y$n / y$rsums
      y <- y$val
      shiny::req(x)
      shiny::req(y)
      if (stats::sd(x) == 0 && stats::sd(y) == 0) {
          return("not enough data")
      }
      stats::t.test(x, y)
  })

  # SOM raster outputs.
  output$somRaster <- shiny::renderPlot({
    xy <- somRasterPlot()
    shiny::req(xy)
    shiny::req(baseRasterGgplot)
    baseRasterGgplot +
      ggplot2::geom_point(
        data = xy, ggplot2::aes(x = x, y = y),
        color = "red", size = 1, inherit.aes = FALSE
      ) +
      ggplot2::geom_point(
        data = xy, ggplot2::aes(x = x, y = y),
        color = "red", shape = 3, size = 1, inherit.aes = FALSE
      )
  })

  output$somRasterSelect <- plotly::renderPlotly({
    rs <- rsUsed()
    if (is.null(rs)) {
      return(NULL)
    }
    p3 <- .buildSOMRasterSelectPlot(somRasterData, rs)
    shiny::req(inherits(p3, "ggplot"))
    quiet_ggplotly(p3, source = "somGrid", tooltip = "") %>%
      plotly::layout(showlegend = FALSE, dragmode = "select") %>%
      plotly::event_register("plotly_selected") %>%
      plotly::event_register("plotly_relayout")
  })

  # FlowSOM marker pies.
  output$flowSOMPieUI <- shiny::renderUI({
    rs <- rsUsed()
    shiny::req(rs)
    ncol_val <- input$flowSOMPieCols
    size_val <- input$flowSOMPieSize
    if (is.null(ncol_val) || is.null(size_val)) {
      return(NULL)
    }
    n <- length(rs)
    if (n < 1L) {
      return(NULL)
    }

    # Limit maximum dimensions to prevent PNG crash (max 50000px)
    ncol_val <- min(ncol_val, 20L)
    size_val <- min(size_val, 300L)

    rows <- ceiling(n / ncol_val)
    max_dim <- 8000L # Conservative limit to stay under 50000px
    total_width <- min(ncol_val * size_val, max_dim)
    total_height <- min(rows * size_val, max_dim)

    shiny::plotOutput(
      "flowSOMPie",
      width = paste0(total_width, "px"),
      height = paste0(total_height, "px")
    )
  })


  # FlowSOM star plot.
  if (!is.null(fsom)) {
    output$flowSOMStars <- plotly::renderPlotly({
      shiny::req(fsom)
      ret <- FlowSOM::PlotStars(
        fsom,
        backgroundValues = fsom$metaclustering,
        list_insteadof_ggarrange = TRUE
      )
      p <- ret$tree
      coords <- fsom$MST$l
      coord_df <- data.frame(
        x  = coords[, 1],
        y  = coords[, 2],
        id = seq_len(nrow(coords))
      )
      p_click <- p +
        ggplot2::geom_point(
          data = coord_df,
          ggplot2::aes(x = x, y = y), # ← customdata removed
          alpha = 0,
          size = 8,
          inherit.aes = FALSE
        )

      last_trace <- length(quiet_ggplotly(p_click)$x$data) # find invisible point trace

      quiet_ggplotly(p_click, source = "flowSOMStars", tooltip = "") %>%
        plotly::style( # ← customdata injected here
          customdata = coord_df$id,
          traces = last_trace
        ) %>%
        plotly::layout(dragmode = "select") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_click")
    })
  }

  output$flowSOMPie <- shiny::renderPlot({
    p <- flowSOMPiePlot()
    shiny::req(p)
    p
  })
  # Violin plots.
  output$VlnPlot <- shiny::renderPlot({
    vlnPlot()
  })

  output$VlnPlot2 <- shiny::renderPlot({
    VlnPlot2()
  })

  # UpSet plot.
  output$UpSet <- shiny::renderPlot({
    upSetPlot()
  })

  # Download handler.
  output$downloadPlots <- shiny::downloadHandler(
    filename = function() {
      nm <- input$clusterNameSelect
      if (is.null(nm) || !nzchar(nm)) nm <- "clusterSelector_plots"
      paste0(nm, ".pdf")
    },
    content = function(file) {
      rs <- rsUsed()
      shiny::req(rs)

      dlColorVar <- input$somColorVar %||% "n"
      dlSizeVar <- input$somSizeVar %||% "max"
      pctl <- input$scatterPercentile %||% 0.99

      upsetSel <- input$upsetSelection
      violinSel <- input$violinSelection %||% colsUsed
      dlOutputList <- rv$outputList
      if (length(upsetSel) < 3) upsetSel <- names(dlOutputList)

      plotRegistry <- list(
        countBar = countBarPlot,
        percentBar = PercentBarPlot,
        tsne = tsnePlot,
        umap = umapPlot,
        pca = pcaPlot,
        scatter = scatterPlot,
        flowSOMPie = flowSOMPiePlot,
        flowSOMStars = if (!is.null(fsom)) {
          function() FlowSOM::PlotStars(fsom, backgroundValues = fsom$metaclustering)$tree
        } else {
          NULL
        },
        violin1 = function() plotViolinFunc(sce, somCodesName, upsetSel, dlOutputList, violinSel),
        violin2 = function() plotViolin2Func(sce, somCodesName, violinSel, upsetSel, dlOutputList),
        upset = function() upsetPlotFunc(upsetSel, dlOutputList, sce)
      )

      .writeClusterSelectorPdf(
        file = file, rs = rs, sce = sce, sceRN = sceRN, sceCN = sceCN,
        metaD = metaD, somCodesName = somCodesName,
        dendPlotObj = dendPlot(),
        dimPairs = dListRV(), getBasePlotFn = getBasePlot,
        dlColorVar = dlColorVar, dlSizeVar = dlSizeVar, pctl = pctl,
        somRasterXy = somRasterPlot(), baseRasterGgplot = baseRasterGgplot,
        plotRegistry = plotRegistry
      )
    }
  )
}
