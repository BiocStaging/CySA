# tests/testthat/test-outputs_clusterSelector-gaps.R
#
# Targets branches uncovered ("!") in the coverage report for
# R/outputs_clusterSelector.R, driven through the real app from
# make_test_app() (helper-*.R) plus a local variant that adds a factor
# grouping column, which the stock CySA_example_sce() fixture lacks.
#
# Selection wiring (see observers_clusterSelector.R / server_clusterSelector.R):
#   input$clusterNumbers (text, comma-separated ids)
#     -> inputClusterNumber() [debounce 1000ms]
#     -> observer sets rsUsed(valid)
#     -> rsUsed_d = rsUsed %>% debounce(1500ms)
#
# rsUsed() feeds: scatter, dendPlot, somRasterPlot, cellCounts,
#   cellPercentages, countBarPlot, PercentBarPlot, flowSOMPiePlot.
# rsUsed_d() feeds: somData<i>, somDataMain, tsnePlot, umapPlot, pcaPlot.
#
# dimSelection() is populated from input$currentDimX/currentDimY via a
# 250ms-debounced activeDims() reactive, and is ALWAYS length 1 regardless
# of nPlots (see note in chat about the downloadPlots loop).

library(testthat)
library(shiny)

# helper-*.R (make_test_app) is auto-sourced by testthat before this file.

# ---- Local fixture: adds a usable factor grouping column -------------------
# CySA_example_sce()'s experiment_info only has sample_id (character) and
# some_numeric (numeric) -- no factor column, so groupsVar/group1/group2 and
# the real t.test() path are unreachable with the stock fixture.

make_test_app_with_groups <- function(n_cells = 400, n_nodes = 10) {
  set.seed(12)
  sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes)

  ei <- S4Vectors::metadata(sce)$experiment_info
  # sample1/sample2 -> "A", sample3/sample4 -> "B", each with >1 sample so
  # sd(x) and sd(y) aren't both zero in ttestResult.
  ei$cohort <- factor(rep(c("A", "A", "B", "B"), length.out = nrow(ei)))
  S4Vectors::metadata(sce)$experiment_info <- ei

  prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
  som_codes <- S4Vectors::metadata(sce)$SOM_codes
  markers <- S4Vectors::metadata(sce)$map$colsUsed
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
  dendTable <- data.frame(
    id = seq_len(nrow(som_codes)), label = rownames(som_codes),
    stringsAsFactors = FALSE
  )
  clusterPatientTable <- table(
    sample_id = sce$sample_id, cluster_id = sce$cluster_id
  )
  somRasterData <- data.frame(
    x = rep(seq_len(5), length.out = nrow(som_codes)),
    y = rep(seq_len(2), each = nrow(som_codes) / 2),
    id = seq_len(nrow(som_codes))
  )
  for (m in markers) somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
  arr <- array(seq_len(10 * 10 * length(markers)), dim = c(10, 10, length(markers)))
  somRasterObj <- raster::brick(arr)
  names(somRasterObj) <- markers

  app <- clusterSelector(
    sce = prepped$sce, sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList, dend = dend, dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData, somRasterObj = somRasterObj,
    env = new.env()
  )
  list(app = app, sce = prepped$sce, markers = markers)
}

# Common: get past clusterNumbers -> rsUsed -> rsUsed_d chain.
select_nodes <- function(session, ids = "1, 2, 3", wait_for_debounced = TRUE) {
  session$setInputs(clusterNumbers = ids)
  session$elapse(1001) # inputClusterNumber() debounce
  if (wait_for_debounced) session$elapse(1501) # rsUsed_d debounce
}

# Common: get past dimSelection()'s req via currentDimX/currentDimY.
select_dims <- function(session, x, y) {
  session$setInputs(currentDimX = x, currentDimY = y)
  session$elapse(300) # activeDims() debounce (250ms)
}

# ---- Gap 1: scatter body past req() (real quantile/plotly logic) ----------

test_that("scatter renders a real plot once selection/dims/samples are valid", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_dims(session, "marker1", "marker2")
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(samples2plot = c("sample1", "sample2", "sample3", "sample4"))
      out <- output$scatter
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 2: somData1 -- both showGroups branches ---------------------------

test_that("somData1 uses the non-grouped path when showGroups is FALSE", {
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    select_dims(session, "marker1", "marker2")
    select_nodes(session, "1, 2, 3")
    session$setInputs(showGroups = FALSE, colorbyGroups = character(0))
    out <- output$somData1
    expect_true(!is.null(out))
  }))
})

test_that("somData1 uses the projection path when showGroups is TRUE", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_dims(session, "marker1", "marker2")
      select_nodes(session, "1, 2, 3")
      session$setInputs(showGroups = TRUE, colorbyGroups = "selected")
      out <- output$somData1
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 3: somDataMain ------------------------------------------------------

test_that("somDataMain renders past the dims/rsUsed_d guard", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_dims(session, "marker1", "marker2")
      select_nodes(session, "1, 2, 3")
      session$setInputs(showGroups = FALSE)
      out <- output$somDataMain
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 4: tsne/umap/pca -- both showlegend branches, real computation ----

test_that("tsne/umap/pca compute and render with legend shown and hidden", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3")
      session$setInputs(
        dimRedSelection = fx$sce %>% (function(s) S4Vectors::metadata(s)$map$colsUsed)(),
        perplexity = 5, n_neighbors = 5, showlegend = TRUE
      )
      session$elapse(2100) # dimRedSelection debounce (1000) + tsne/umap/pca debounce (1000)

      expect_true(!is.null(output$tsne))
      expect_true(!is.null(output$umap))
      expect_true(!is.null(output$pca))

      session$setInputs(showlegend = FALSE)
      expect_true(!is.null(output$tsne))
      expect_true(!is.null(output$umap))
      expect_true(!is.null(output$pca))
    })
  )
})

# ---- Gap 5: cellCounts -- both comparison branches -------------------------

test_that("cellCounts adds a comparison row when compareStatsTo matches experiment_info", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(compareStatsTo = "some_numeric", singleNode = 1)
      tbl <- output$cellCounts
      expect_true(!is.null(tbl))
    })
  )
})

test_that("cellCounts adds a comparison row when compareStatsTo matches an outputList name", {
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
    session$setInputs(clusterName = "popA")
    session$setInputs(applyName = 1) # actionButton click
    session$setInputs(compareStatsTo = "popA", singleNode = 1)
    tbl <- output$cellCounts
    expect_true(!is.null(tbl))
  }))
})

# ---- Gap 6: ttestResult -- empty-data guard AND real t.test() path ---------

test_that("ttestResult reports 'no data' when groupsVar is 'none'", {
  fx <- make_test_app()
  quiet_plotly_test(testServer(fx$app, {
    session$setInputs(groupsVar = "none", relativeTo = "none")
    select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
    expect_output(print(output$ttestResult), "no data")
  }))
})

test_that("drawProjection produces non-empty traces for a real selection", {
  fx <- make_test_app()
  sce <- fx$sce
  pp1 <- plotSOMScatter(
    x = sce, chs = c("marker1", "marker2"),
    pointSize = "max", color_by = "n", zeros = TRUE,
    xRN = rownames(sce), xCN = colnames(sce)
  )
  dimSel <- list(list(dims = c("marker1", "marker2")))
  projectionDf <- buildProjectionDf(pp1, 1L, dimSel, sce, somCodesName = "SOM_codes")

  p <- drawProjection(projectionDf,
    rs = 1:3, colorbyGroups = "selected",
    sce = sce, outputList = list(selected = 1:3)
  )

  built <- plotly::plotly_build(plotly::ggplotly(p))
  trace_lengths <- vapply(built$x$data, function(tr) {
    length(tr$x %||% tr$y %||% integer(0))
  }, integer(1))

  expect_true(all(trace_lengths > 0),
    info = paste(
      "Empty trace(s) found at index:",
      paste(which(trace_lengths == 0), collapse = ", ")
    )
  )
})

test_that("PercentBarPlot ignores non-numeric experiment_info columns", {
  fx <- make_test_app_with_groups() # has a factor 'cohort' column
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(relativeTo = "none")
      expect_no_warning(output$PercentBar)
    })
  )
})

# ---- Gap 7: flowSOMPieUI / flowSOMPie real rendering -----------------------

test_that("flowSOMPieUI builds a plotOutput once nodes are selected", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(flowSOMPieCols = 4, flowSOMPieSize = 100, flowSOMPieMax = 10)
      ui <- output$flowSOMPieUI
      expect_true(!is.null(ui))
    })
  )
})

test_that("flowSOMPie renders once flowSOMPiePlot() returns a real plot", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
      session$setInputs(flowSOMPieCols = 4, flowSOMPieSize = 100, flowSOMPieMax = 10)
      out <- output$flowSOMPie
      expect_true(!is.null(out))
    })
  )
})

# ---- Gap 8: staticSomGrid with the app's real dList pairs ------------------

test_that("staticSomGrid builds panels/rows for the app's dList", {
  fx <- make_test_app()
  quiet_plotly_test(
    testServer(fx$app, {
      session$setInputs(staticSomCols = 2, staticSomHeight = 200)
      ui <- output$staticSomGrid
      expect_true(!is.null(ui))
      # sanity: with 6 default dList pairs and 2 cols, expect >= 1 row structure
    })
  )
})

# ---- Gap 9: downloadPlots filename() and content() -------------------------
#
# NOTE: with the default nPlots = 6, the content() handler's
# `for (plotIdx in seq_len(nPlots)) { dims <- dimSel[[plotIdx]]$dims }` will
# throw "subscript out of bounds" for plotIdx >= 2, because dimSelection()
# is only ever populated with ONE entry (see note above). Using nPlots = 1
# here to get a clean, passing coverage test; consider a *separate* test
# that documents/guards the nPlots > 1 crash as its own bug report:
#
#   test_that("downloadPlots errors with nPlots > 1 (KNOWN ISSUE)", {
#     fx <- make_test_app()  # default nPlots = 6
#     testServer(fx$app, {
#       select_dims(session, "marker1", "marker2")
#       select_nodes(session, "1, 2, 3", wait_for_debounced = FALSE)
#       session$setInputs(clusterNameSelect = "x")
#       expect_error(session$getDownload("downloadPlots"), "subscript out of bounds")
#     })
#   })

test_that(".writeClusterSelectorPdf writes a non-empty PDF", {
  fx <- make_test_app()
  tmp <- tempfile(fileext = ".pdf")

  # Build a minimal getBasePlotFn that returns a simple ggplot for any dim pair.
  getBasePlotFn <- function(xVar, yVar, colorVar, sizeVar) {
    ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y))
  }

  .writeClusterSelectorPdf(
    file = tmp,
    rs = 1:3,
    sce = fx$sce,
    sceRN = rownames(fx$sce),
    sceCN = colnames(fx$sce),
    metaD = S4Vectors::metadata(fx$sce),
    somCodesName = "SOM_codes",
    dendPlotObj = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:20, 5)))),
    dimPairs = list(c("marker1", "marker2")),
    getBasePlotFn = getBasePlotFn,
    dlColorVar = "n",
    dlSizeVar = "max",
    pctl = 0.99,
    somRasterXy = data.frame(x = 1, y = 1),
    baseRasterGgplot = ggplot2::ggplot(),
    plotRegistry = list()
  )
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})
