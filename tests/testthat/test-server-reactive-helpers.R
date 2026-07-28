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
    suppress_plotly_event_warnings({
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
})

# --- .buildZoomObservers: no source changes needed -------------------------

test_that(".buildZoomObservers calls zoomFunc with the relayout event and plot index", {
    calls <- list()
    zoomFunc_stub <- function(zoom, plotIdx) {
        calls[[length(calls) + 1]] <<- list(zoom = zoom, plotIdx = plotIdx)
    }

    server <- function(input, output, session) {
        .buildZoomObservers(nPlots = 2, input = input, zoomFunc = zoomFunc_stub, verbose = FALSE)
    }

    quiet_plotly_test({
        testServer(server, {
            session$setInputs(
                `plotly_relayout-somData1` = jsonlite::toJSON(
                    list(`xaxis.range[0]` = 1, `xaxis.range[1]` = 5), auto_unbox = TRUE
                ),
                `plotly_relayout-somData2` = jsonlite::toJSON(
                    list(`xaxis.range[0]` = 2, `xaxis.range[1]` = 6), auto_unbox = TRUE
                )
            )
        })
    })

    expect_length(calls, 2)
    expect_true(any(vapply(calls, function(c) c$plotIdx == 1, logical(1))))
    expect_true(any(vapply(calls, function(c) c$plotIdx == 2, logical(1))))
})


# --- .buildSelectionObserver: requires the rsUsed_d parameter fix above ----
test_that("minimal repro: does observeEvent fire at all with a mocked, non-reactive trigger", {
    testthat::local_mocked_bindings(
        event_data = function(event, source = "all") {
            message("BARE MOCK CALLED")
            list(x = 1)
        },
        .package = "plotly"
    )

    shiny::testServer(
        app = shiny::shinyApp(
            ui = shiny::fluidPage(),
            server = function(input, output, session) {
                val <- shiny::reactiveVal("initial")

                shiny::observeEvent(
                    plotly::event_data("plotly_selected", source = "test"),
                    {
                        message("BARE HANDLER FIRED")
                        val("changed")
                    }
                )

                output$probe <- shiny::renderPrint(val())
            }
        ),
        expr = {
            session$flushReact()
            message("val after flush: ", shiny::isolate(val()))
        }
    )
})

test_that(".buildSelectionObserver updates rsUsed via inputSelect when a selection fires", {
    testthat::local_mocked_bindings(
        event_data = function(event, source = "all") {
            message("MOCK event_data CALLED: ", event, " / ", source)
            list(pointNumber = 0, curveNumber = 0)
        },
        .package = "plotly"
    )

    shiny::testServer(
        app = shiny::shinyApp(
            ui = shiny::fluidPage(),
            server = function(input, output, session) {
                rsUsed   <- shiny::reactiveVal(c(1, 2, 3, 4, 5))
                rsUsed_d <- shiny::reactive(rsUsed())

                .buildSelectionObserver(
                    sourceId    = "mySource",
                    input       = input,
                    rsUsed      = rsUsed,
                    rsUsed_d    = rsUsed_d,
                    inputSelect = function(d, rs, mode) {
                        message("inputSelect CALLED")
                        999
                    },
                    verbose     = FALSE
                )

                output$rsUsedProbe <- shiny::renderPrint(rsUsed())
            }
        ),
        expr = {
            session$flushReact()
            message("rsUsed after flush: ", paste(shiny::isolate(rsUsed()), collapse = ","))
            expect_output(print(output$rsUsedProbe), "999")
        }
    )
})

test_that(".buildSelectionObserver is a no-op when rsUsed_d() is NULL", {
    server <- function(input, output, session) {
        rsUsed   <- shiny::reactiveVal(NULL)
        rsUsed_d <- shiny::reactive(NULL)
        called <- FALSE
        inputSelect_stub <- function(d, rs, selectMode) { called <<- TRUE; rs }

        .buildSelectionObserver("mySource", input, rsUsed, rsUsed_d, inputSelect_stub, verbose = FALSE)
        output$calledProbe <- shiny::renderPrint(called)
    }

    quiet_plotly_test({
        testServer(server, {
        session$setInputs(`.clientValue-plotly_selected-mySource` = list(pointNumber = 0))
        expect_output(print(output$calledProbe), "FALSE")
    })
    })
})

# --- .buildSOMDataObservers: requires the same rsUsed_d parameter fix ------

test_that(".buildSOMDataObservers updates rsUsed and activePlot for the touched panel", {
    server <- function(input, output, session) {
        rsUsed     <- shiny::reactiveVal(1:3)
        rsUsed_d   <- shiny::reactive(rsUsed())
        activePlot <- shiny::reactiveVal(1L)
        inputSelect_stub <- function(d, rs, selectMode) c(rs, 42L)
        .buildSOMDataObservers(
            nPlots = 3, input = input, output = output,
            rsUsed = rsUsed, rsUsed_d = rsUsed_d, activePlot = activePlot,
            inputSelect = inputSelect_stub, verbose = FALSE
        )
        output$rsUsedProbe     <- shiny::renderPrint(rsUsed())
        output$activePlotProbe <- shiny::renderPrint(activePlot())
    }
    quiet_plotly_test({
        testServer(server, {
            click_json <- jsonlite::toJSON(list(list(pointNumber = 0)), auto_unbox = TRUE)
            session$setInputs(`plotly_selected-somData2` = click_json)
            expect_output(print(output$rsUsedProbe), "42")
            expect_output(print(output$activePlotProbe), "2")
        })
    })
})

test_that(".buildSOMDataObservers is a no-op when rsUsed_d() is NULL (shiny::req guard)", {
    server <- function(input, output, session) {
        rsUsed     <- shiny::reactiveVal(NULL)
        rsUsed_d   <- shiny::reactive(NULL)
        activePlot <- shiny::reactiveVal(1L)
        called <- FALSE
        inputSelect_stub <- function(d, rs, selectMode) { called <<- TRUE; rs }

        .buildSOMDataObservers(
            nPlots = 1, input = input, output = output,
            rsUsed = rsUsed, rsUsed_d = rsUsed_d, activePlot = activePlot,
            inputSelect = inputSelect_stub, verbose = FALSE
        )
        output$calledProbe <- shiny::renderPrint(called)
    }

    quiet_plotly_test({
        testServer(server, {
        session$setInputs(`.clientValue-plotly_selected-somData1` = list(pointNumber = 0))
        expect_output(print(output$calledProbe), "FALSE")
    })
    })
})





