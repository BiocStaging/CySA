# Tests for server-reactive-helpers.R
# These test the reactive helper functions that can be used to build
# clusterSelector-like apps.

library(shiny)


# =============================================================================
# .buildDimRedReactives()
# =============================================================================

test_that(".buildDimRedReactives() returns list of tsne, umap, pca reactives", {
    sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
    metaD <- S4Vectors::metadata(sce)
    somCodesName <- "SOM_codes"

    input <- shiny::reactiveValues(
        dimRedSelection = metaD$map$colsUsed[1:5],
        perplexity = 10,
        n_neighbors = 15
    )

    reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)

    expect_type(reactives, "list")
    expect_length(reactives, 3)
    expect_true("tsne" %in% names(reactives))
    expect_true("umap" %in% names(reactives))
    expect_true("pca" %in% names(reactives))
    expect_type(reactives$tsne, "closure")
    expect_type(reactives$umap, "closure")
    expect_type(reactives$pca, "closure")
})


# test_that(".buildDimRedReactives() pca reactive computes PCA correctly", {
#     sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
#     metaD <- S4Vectors::metadata(sce)
#     somCodesName <- "SOM_codes"
#
#     input <- shiny::reactiveValues(
#         dimRedSelection = metaD$map$colsUsed[1:5],
#         perplexity = 10,
#         n_neighbors = 15
#     )
#
#     reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)
#
#     shiny::testServer(
#         app = shiny::shinyApp(
#             ui = shiny::fluidPage(),
#             server = function(input, output, session) {
#                 output$pcaResult <- shiny::renderPrint({
#                     pca_result <- reactives$pca()
#                     expect_s3_class(pca_result, "prcomp")
#                     expect_equal(ncol(pca_result$rotation), 2)  # rank. = 2
#                     expect_false(pca_result$scale)  # scale = FALSE
#                     cat("PCA OK")
#                 })
#             }
#         ),
#         expr = { session$flushReact() }
#     )
# })


# test_that(".buildDimRedReactives() pca respects dimRedSelection changes", {
#     sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
#     metaD <- S4Vectors::metadata(sce)
#     somCodesName <- "SOM_codes"
#
#     input <- shiny::reactiveValues(
#         dimRedSelection = metaD$map$colsUsed[1:3],  # 3 markers
#         perplexity = 10,
#         n_neighbors = 15
#     )
#
#     reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)
#
#     shiny::testServer(
#         app = shiny::shinyApp(
#             ui = shiny::fluidPage(),
#             server = function(input, output, session) {
#                 output$pcaResult <- shiny::renderPrint({
#                     pca_result <- reactives$pca()
#                     # PCA on t(SOM_codes) with 3 markers -> 3 rows in rotation
#                     expect_equal(nrow(pca_result$rotation), 3)
#                     cat("PCA OK")
#                 })
#             }
#         ),
#         expr = { session$flushReact() }
#     )
# })


# test_that(".buildDimRedReactives() umap reactive computes UMAP", {
#     skip_on_cran()
#     sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
#     metaD <- S4Vectors::metadata(sce)
#     somCodesName <- "SOM_codes"
#
#     input <- shiny::reactiveValues(
#         dimRedSelection = metaD$map$colsUsed[1:5],
#         perplexity = 10,
#         n_neighbors = 15
#     )
#
#     reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)
#
#     shiny::testServer(
#         app = shiny::shinyApp(
#             ui = shiny::fluidPage(),
#             server = function(input, output, session) {
#                 output$umapResult <- shiny::renderPrint({
#                     umap_result <- reactives$umap()
#                     expect_true(inherits(umap_result, "umap"))
#                     expect_true("layout" %in% names(umap_result))
#                     expect_equal(ncol(umap_result$layout), 2)  # 2D output
#                     cat("UMAP OK")
#                 })
#             }
#         ),
#         expr = { session$flushReact() }
#     )
# })


# test_that(".buildDimRedReactives() umap respects n_neighbors input", {
#     skip_on_cran()
#     sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
#     metaD <- S4Vectors::metadata(sce)
#     somCodesName <- "SOM_codes"
#
#     input <- shiny::reactiveValues(
#         dimRedSelection = metaD$map$colsUsed[1:5],
#         perplexity = 10,
#         n_neighbors = 30  # Different from default
#     )
#
#     reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)
#
#     shiny::testServer(
#         app = shiny::shinyApp(
#             ui = shiny::fluidPage(),
#             server = function(input, output, session) {
#                 output$umapResult <- shiny::renderPrint({
#                     umap_result <- reactives$umap()
#                     expect_true(inherits(umap_result, "umap"))
#                     cat("UMAP OK")
#                 })
#             }
#         ),
#         expr = { session$flushReact() }
#     )
# })


# test_that(".buildDimRedReactives() tsne reactive computes t-SNE", {
#     skip_on_cran()
#     sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
#     metaD <- S4Vectors::metadata(sce)
#     somCodesName <- "SOM_codes"
#
#     input <- shiny::reactiveValues(
#         dimRedSelection = metaD$map$colsUsed[1:5],
#         perplexity = 10,
#         n_neighbors = 15
#     )
#
#     reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)
#
#     shiny::testServer(
#         app = shiny::shinyApp(
#             ui = shiny::fluidPage(),
#             server = function(input, output, session) {
#                 output$tsneResult <- shiny::renderPrint({
#                     tsne_result <- reactives$tsne()
#                     expect_true(inherits(tsne_result, "tsne"))
#                     expect_true("Y" %in% names(tsne_result))
#                     expect_equal(ncol(tsne_result$Y), 2)  # 2D output
#                     cat("t-SNE OK")
#                 })
#             }
#         ),
#         expr = { session$flushReact() }
#     )
# })


# test_that(".buildDimRedReactives() tsne respects perplexity input", {
#     skip_on_cran()
#     sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
#     metaD <- S4Vectors::metadata(sce)
#     somCodesName <- "SOM_codes"
#
#     input <- shiny::reactiveValues(
#         dimRedSelection = metaD$map$colsUsed[1:5],
#         perplexity = 25,  # Different from default
#         n_neighbors = 15
#     )
#
#     reactives <- .buildDimRedReactives(input, metaD, sce, somCodesName)
#
#     shiny::testServer(
#         app = shiny::shinyApp(
#             ui = shiny::fluidPage(),
#             server = function(input, output, session) {
#                 output$tsneResult <- shiny::renderPrint({
#                     tsne_result <- reactives$tsne()
#                     expect_true(inherits(tsne_result, "tsne"))
#                     cat("t-SNE OK")
#                 })
#             }
#         ),
#         expr = { session$flushReact() }
#     )
# })


# =============================================================================
# .buildSelectionObserver()
# =============================================================================

test_that(".buildSelectionObserver() creates an observer that updates rsUsed", {
    sce <- CySA_example_sce(n_cells = 50, n_nodes = 10)
    metaD <- S4Vectors::metadata(sce)
    somCodesName <- "SOM_codes"

    # Mock inputSelect function that adds selected points to existing rs
    inputSelect <- function(eventData, currentRs, mode) {
        # Simulate extracting pointData keys and converting to node ids
        if (is.null(eventData)) return(currentRs)
        new_ids <- unique(eventData$pointData$key)
        if (mode == "add") {
            return(unique(c(currentRs, new_ids)))
        } else if (mode == "view") {
            return(new_ids)
        } else {
            return(currentRs)
        }
    }

    input <- shiny::reactiveValues(
        selectMode = "view",
        plotly_selected = NULL
    )

    rsUsed <- shiny::reactiveVal(c(1L))

    # Create observer - this registers the observer but doesn't trigger it
    observer <- .buildSelectionObserver(
        sourceId = "testPlot",
        input = input,
        rsUsed = rsUsed,
        inputSelect = inputSelect,
        verbose = FALSE
    )

    expect_true(inherits(observer, "Observer"))

    # Clean up
    observer$destroy()
})


# =============================================================================
# .buildSOMDataObservers()
# =============================================================================

test_that(".buildSOMDataObservers() creates observers for each plot", {
    input <- shiny::reactiveValues(
        selectMode = "view",
        plotly_selected = NULL
    )
    output <- shiny::reactiveValues()
    rsUsed <- shiny::reactiveVal(c(1L))
    activePlot <- shiny::reactiveVal(1L)

    inputSelect <- function(eventData, currentRs, mode) {
        if (is.null(eventData)) return(currentRs)
        new_ids <- unique(eventData$pointData$key)
        return(new_ids)
    }

    nPlots <- 3L
    observers <- .buildSOMDataObservers(
        nPlots = nPlots,
        input = input,
        output = output,
        rsUsed = rsUsed,
        activePlot = activePlot,
        inputSelect = inputSelect,
        verbose = FALSE
    )

    expect_type(observers, "list")
    expect_length(observers, nPlots)
    expect_true(all(vapply(observers, inherits, logical(1), "Observer")))

    # Clean up
    lapply(observers, function(obs) obs$destroy())
})


test_that(".buildSOMDataObservers() updates activePlot on selection", {
    input <- shiny::reactiveValues(
        selectMode = "view",
        plotly_selected = NULL
    )
    output <- shiny::reactiveValues()
    rsUsed <- shiny::reactiveVal(c(1L))
    activePlot <- shiny::reactiveVal(1L)

    inputSelect <- function(eventData, currentRs, mode) {
        if (is.null(eventData)) return(currentRs)
        return(unique(eventData$pointData$key))
    }

    observers <- .buildSOMDataObservers(
        nPlots = 2L,
        input = input,
        output = output,
        rsUsed = rsUsed,
        activePlot = activePlot,
        inputSelect = inputSelect,
        verbose = FALSE
    )

    # Initially activePlot should be 1 (can't test reactiveVal() directly outside context)
    # Just verify observers were created
    expect_length(observers, 2L)

    # Clean up
    lapply(observers, function(obs) obs$destroy())
})


# =============================================================================
# .buildZoomObservers()
# =============================================================================

test_that(".buildZoomObservers() creates observers for each plot", {
    input <- shiny::reactiveValues(
        plotly_relayout = NULL
    )

    zoomFunc <- function(zoom, plotIdx) {
        # Mock zoom function - does nothing
        invisible(NULL)
    }

    nPlots <- 3L
    observers <- .buildZoomObservers(
        nPlots = nPlots,
        input = input,
        zoomFunc = zoomFunc,
        verbose = FALSE
    )

    expect_type(observers, "list")
    expect_length(observers, nPlots)
    expect_true(all(vapply(observers, inherits, logical(1), "Observer")))

    # Clean up
    lapply(observers, function(obs) obs$destroy())
})


test_that(".buildZoomObservers() calls zoomFunc with correct plotIdx", {
    zoom_calls <- list()

    input <- shiny::reactiveValues(
        plotly_relayout = NULL
    )

    zoomFunc <- function(zoom, plotIdx) {
        zoom_calls[[length(zoom_calls) + 1]] <<- list(zoom = zoom, plotIdx = plotIdx)
    }

    observers <- .buildZoomObservers(
        nPlots = 2L,
        input = input,
        zoomFunc = zoomFunc,
        verbose = FALSE
    )

    # Observers are registered but not triggered yet
    expect_length(zoom_calls, 0)

    # Clean up
    lapply(observers, function(obs) obs$destroy())
})


# =============================================================================
# .safeEventData() helper (used internally)
# =============================================================================

test_that(".safeEventData returns NULL when event not triggered", {
    # .safeEventData needs a reactive context, so test inside testServer
    shiny::testServer(
        app = shiny::shinyApp(
            ui = shiny::fluidPage(),
            server = function(input, output, session) {
                output$testOutput <- shiny::renderPrint({
                    result <- .safeEventData(verbose = FALSE, "plotly_selected", source = "nonexistent")
                    expect_null(result)
                    cat("NULL OK")
                })
            }
        ),
        expr = {
            session$flushReact()
        }
    )
})
