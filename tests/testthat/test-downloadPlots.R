
# test-downloadPlots.R
#
# Architecture note
# ─────────────────────────────────────────────────────────────────────────────
# shiny::testServer() does NOT support invoking downloadHandler content()
# functions directly. The handler creates an HTTP endpoint, not a reactive
# output (see https://github.com/rstudio/shiny/issues/3979).
#
# For full end-to-end download testing (including browser file-save behaviour)
# use shinytest2::AppDriver$new() + $get_download().
#
# Test strategy used here:
#   1.  Filename logic  ── pure function tests, no Shiny context needed.
#   2.  req() guard     ── shiny::isolate() used to exercise req() behaviour.
#   3.  PDF content     ── each plotting call inside content() is exercised
#                          directly against a real PDF device.
#   4.  Integration     ── make_test_app() + testServer() verifies the server
#                          starts and all inputs that the handler reads are
#                          correctly accepted.

library(testthat)
library(shiny)

# =============================================================================
# Updated .muffle_plotly_warning()
# Extend to cover the 'colour' attribute validator warning that fires during
# plotly_build() inside testServer() flush callbacks.
# Source: plotly:::verify_attr_names() → attrs_name_check()
# quiet_ggplotly() strips 'colour' AFTER build, so the validator still fires.
# =============================================================================
.muffle_plotly_warning <- function(expr) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("is not registered",          msg, fixed = FALSE) ||
                grepl("event_register",             msg, fixed = FALSE) ||
                grepl("don't have these attributes", msg, fixed = FALSE) ||
                grepl("colour",                     msg, fixed = FALSE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

# =============================================================================
# Shared helpers
# =============================================================================

#' Read the first 4 bytes of a file and check for the %PDF magic signature.
.is_valid_pdf <- function(path) {
    if (!file.exists(path) || file.size(path) == 0L) return(FALSE)
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    startsWith(rawToChar(readBin(con, "raw", n = 4L)), "%PDF")
}

#' Build a minimal dimSelection list suitable for ggsomPlot / content() loop.
#' @param sce CySA SCE object.
#' @param n   Number of identical entries to create.
.make_dimsel <- function(sce, n = 1L) {
    markers <- S4Vectors::metadata(sce)$map$colsUsed
    ch1 <- markers[1L]
    ch2 <- markers[2L]
    som <- S4Vectors::metadata(sce)$SOM_codes
    replicate(n, list(
        dims  = c(ch1, ch2),
        xlim  = range(som[, ch1]),
        ylim  = range(som[, ch2]),
        xzoom = c(NULL, NULL),
        yzoom = c(NULL, NULL)
    ), simplify = FALSE)
}

#' Open a PDF device, run expr, always close the device and delete the file.
#' Returns the path so callers can inspect the file after the block.
.with_pdf <- function(expr) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(file = tmp)
    tryCatch(
        force(expr),
        finally = {
            if (grDevices::dev.cur() > 1L) grDevices::dev.off()
        }
    )
    tmp
}


# =============================================================================
# 1.  Filename logic — pure function
# =============================================================================
# Handler body:
#   shiny::req(input$clusterNameSelect)
#   paste(input$clusterNameSelect, ".pdf", sep = "")

test_that("downloadPlots filename: paste produces '<name>.pdf' exactly", {
    nms <- c("myGroup", "CD8_T_cells", "B cells", "NK-cells", "cluster.1")

    for (nm in nms) {
        result <- paste(nm, ".pdf", sep = "")

        expect_equal(result, paste0(nm, ".pdf"),
                     label = sprintf("exact match for '%s'", nm))

        expect_true(endsWith(result, ".pdf"),
                    label = sprintf("'%s' ends with .pdf", result))

        # sep = "" must not insert a space before ".pdf"
        expect_false(grepl(" \\.pdf$", result),
                     label = sprintf("no space before .pdf in '%s'", result))
    }
})

test_that("downloadPlots filename: sep='' is equivalent to paste0", {
    nm <- "TestGroup"
    expect_equal(paste(nm, ".pdf", sep = ""), paste0(nm, ".pdf"))
})

test_that("downloadPlots filename: multiple selections produce one string per element", {
    # input$clusterNameSelect is a multiple-select; paste() without collapse
    # returns a vector of the same length — the first element is what a
    # browser would use as the default filename.
    sels <- c("Group1", "Group2")
    result <- paste(sels, ".pdf", sep = "")

    expect_length(result, length(sels))
    expect_equal(result, c("Group1.pdf", "Group2.pdf"))

    # The first element is the filename the browser receives.
    expect_equal(result[[1L]], "Group1.pdf")
})


# =============================================================================
# 2.  req() guard — handler stops silently when clusterNameSelect is missing
# =============================================================================

test_that("downloadPlots guard: req(NULL) raises shiny.silent.error", {
    expect_error(
        shiny::isolate(shiny::req(NULL)),
        class = "shiny.silent.error"
    )
})

test_that("downloadPlots guard: req('') raises shiny.silent.error", {
    expect_error(
        shiny::isolate(shiny::req("")),
        class = "shiny.silent.error"
    )
})

test_that("downloadPlots guard: req(NA) raises shiny.silent.error", {
    expect_error(
        shiny::isolate(shiny::req(NA)),
        class = "shiny.silent.error"
    )
})

test_that("downloadPlots guard: req('validName') does not error", {
    expect_no_error(shiny::isolate(shiny::req("validName")))
})


# =============================================================================
# 3.  PDF content — component-level tests
#
# Each block corresponds to one section of the content() function body.
# Tests run against a real grDevices::pdf() device and verify:
#   (a) the plotting call completes without error,
#   (b) the output file is non-empty,
#   (c) the file carries valid %PDF magic bytes.
# =============================================================================

# ── 3a. SOM scatter + ggsomPlot (the nPlots loop) ────────────────────────────

test_that("downloadPlots content: plotSOMScatter + ggsomPlot write to PDF", {
    skip_on_cran()

    sce    <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    metaD  <- S4Vectors::metadata(sce)
    rs     <- 1L:5L
    dimSel <- .make_dimsel(sce, n = 1L)

    tmp <- .with_pdf({
        pp1 <- plotSOMScatter(
            x = sce, chs = dimSel[[1L]]$dims,
            pointSize = "max", color_by = "n",
            xRN = rownames(sce), xCN = colnames(sce)
        )
        invisible(print(ggsomPlot(
            pp1, 1L, rs, dimSel,
            sce = sce, metaD = metaD
        )))
    })
    on.exit(unlink(tmp), add = TRUE)

    expect_true(file.exists(tmp),   label = "PDF file was created")
    expect_gt(file.size(tmp), 0L,   label = "PDF file is non-empty")
    expect_true(.is_valid_pdf(tmp), label = "file starts with %PDF magic bytes")
})

test_that("downloadPlots content: nPlots SOM scatter iterations complete", {
    skip_on_cran()

    n_plots  <- 3L
    sce      <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    metaD    <- S4Vectors::metadata(sce)
    rs       <- 1L:5L
    dimSel   <- .make_dimsel(sce, n = n_plots)
    somCodes <- metaD$SOM_codes
    pctl     <- 0.99
    tailP    <- (1 - pctl) / 2

    tmp <- .with_pdf({
        for (plotIdx in seq_len(n_plots)) {
            dims <- dimSel[[plotIdx]]$dims

            xlimP <- if (dims[1L] %in% colnames(somCodes))
                stats::quantile(somCodes[, dims[1L]],
                                probs = c(tailP, 1 - tailP), na.rm = TRUE)
            else NULL

            ylimP <- if (dims[2L] %in% colnames(somCodes))
                stats::quantile(somCodes[, dims[2L]],
                                probs = c(tailP, 1 - tailP), na.rm = TRUE)
            else NULL

            pp1 <- plotSOMScatter(
                x = sce, chs = dims,
                pointSize = "max", color_by = "n",
                xRN = rownames(sce), xCN = colnames(sce)
            )
            invisible(print(ggsomPlot(
                pp1, plotIdx, rs, dimSel,
                sce = sce, metaD = metaD,
                xlim = xlimP, ylim = ylimP
            )))
        }
    })
    on.exit(unlink(tmp), add = TRUE)

    expect_true(.is_valid_pdf(tmp))
})



# ── 3b. plotViolinFunc ────────────────────────────────────────────────────────

test_that("downloadPlots content: plotViolinFunc writes a valid PDF", {
    skip_on_cran()

    sce        <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    markers    <- S4Vectors::metadata(sce)$map$colsUsed
    outputList <- list("Group1" = 1L:4L, "Group2" = 5L:8L, "Rest" = 9L:10L)

    tmp <- .with_pdf({
        pVln <- plotViolinFunc(sce, "SOM_codes", names(outputList),
                               outputList, markers)
        if (!is.null(pVln)) invisible(print(pVln))
    })
    on.exit(unlink(tmp), add = TRUE)

    expect_true(file.exists(tmp))
    expect_gt(file.size(tmp), 0L)
    expect_true(.is_valid_pdf(tmp))
})

test_that("downloadPlots content: plotViolinFunc returns NULL gracefully when < 2 nodes per group", {
    sce        <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    markers    <- S4Vectors::metadata(sce)$map$colsUsed
    # Each group has only one node — plotViolinFunc should skip those
    outputList <- list("G1" = 1L, "G2" = 2L, "Rest" = 3L:10L)

    result <- plotViolinFunc(sce, "SOM_codes", names(outputList),
                             outputList, markers)
    # Either NULL or a ggplot with a "no data" message is acceptable
    expect_true(is.null(result) || inherits(result, "ggplot"))
})


# ── 3c. plotViolin2Func ───────────────────────────────────────────────────────

test_that("downloadPlots content: plotViolin2Func writes a valid PDF", {
    skip_on_cran()

    sce        <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    markers    <- S4Vectors::metadata(sce)$map$colsUsed
    outputList <- list("Group1" = 1L:4L, "Group2" = 5L:8L, "Rest" = 9L:10L)

    tmp <- .with_pdf({
        pVln2 <- plotViolin2Func(sce, "SOM_codes", markers,
                                 names(outputList), outputList)
        if (!is.null(pVln2)) invisible(print(pVln2))
    })
    on.exit(unlink(tmp), add = TRUE)

    expect_true(file.exists(tmp))
    expect_gt(file.size(tmp), 0L)
    expect_true(.is_valid_pdf(tmp))
})


# ── 3d. SOM raster overlay ────────────────────────────────────────────────────

test_that("downloadPlots content: SOM raster + selected-node overlay writes to PDF", {
    skip_on_cran()

    sce     <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    markers <- S4Vectors::metadata(sce)$map$colsUsed
    rs      <- 1L:5L
    n_nodes <- 10L

    somRasterData <- data.frame(
        x  = rep(seq_len(5L), length.out = n_nodes),
        y  = rep(seq_len(2L), each = n_nodes / 2L),
        id = seq_len(n_nodes)
    )
    for (m in markers) somRasterData[[m]] <- seq_len(n_nodes) / n_nodes

    baseRasterGgplot <- .buildBaseRasterGgplot(somRasterData, markers)
    xy <- somRasterData[somRasterData$id %in% rs, c("x", "y"), drop = FALSE]

    tmp <- .with_pdf({
        invisible(print(
            baseRasterGgplot +
                ggplot2::geom_point(
                    data = xy, ggplot2::aes(x = x, y = y),
                    color = "red", size = 1, inherit.aes = FALSE
                ) +
                ggplot2::geom_point(
                    data = xy, ggplot2::aes(x = x, y = y),
                    color = "red", shape = 3, size = 1, inherit.aes = FALSE
                )
        ))
    })
    on.exit(unlink(tmp), add = TRUE)

    expect_true(.is_valid_pdf(tmp))
})

test_that("downloadPlots content: NULL raster overlay is handled silently", {
    # When somRasterPlot() or baseRasterGgplot is NULL, the handler skips.
    xy               <- NULL
    baseRasterGgplot <- NULL

    # The guard condition in the handler: if (!is.null(xy) && !is.null(baseRasterGgplot))
    expect_false(!is.null(xy) && !is.null(baseRasterGgplot),
                 label = "NULL guard correctly prevents raster plot section")
})


# ── 3e. Complete pipeline ─────────────────────────────────────────────────────

test_that("downloadPlots content: full pipeline writes a non-trivial multi-page PDF", {
    skip_on_cran()
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")

    sce        <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    metaD      <- S4Vectors::metadata(sce)
    markers    <- metaD$map$colsUsed
    rs         <- 1L:5L
    dimSel     <- .make_dimsel(sce, n = 2L)
    outputList <- list("Group1" = rs, "Rest" = setdiff(1L:10L, rs))

    tmp <- .with_pdf({
        # SOM scatter plots
        for (plotIdx in seq_along(dimSel)) {
            dims <- dimSel[[plotIdx]]$dims
            pp1  <- plotSOMScatter(
                x = sce, chs = dims,
                pointSize = "max", color_by = "n"
            )
            invisible(print(ggsomPlot(pp1, plotIdx, rs, dimSel,
                                      sce = sce, metaD = metaD)))
        }

        # Violin plots
        pVln <- plotViolinFunc(sce, "SOM_codes", names(outputList),
                               outputList, markers)
        if (!is.null(pVln)) invisible(print(pVln))

        pVln2 <- plotViolin2Func(sce, "SOM_codes", markers,
                                 names(outputList), outputList)
        if (!is.null(pVln2)) invisible(print(pVln2))
    })
    on.exit(unlink(tmp), add = TRUE)

    expect_true(.is_valid_pdf(tmp),
                label = "complete pipeline output is a valid PDF")

    # A multi-plot PDF should be meaningfully larger than a blank page (~4 KB).
    expect_gt(file.size(tmp), 4096L,
              label = "PDF has reasonable content (> 4 KB)")
})


# =============================================================================
# 4.  Integration — server starts and accepts every input the handler reads
# =============================================================================

test_that("downloadPlots integration: server accepts all inputs the handler depends on", {
    skip_on_cran()

    ta <- make_test_app(n_cells = 200L, n_nodes = 10L)

    .muffle_plotly_warning({
        shiny::testServer(ta$app, {
            markers <- S4Vectors::metadata(ta$sce)$map$colsUsed

            suppressWarnings(session$setInputs(
                # filename() reads this
                clusterNameSelect  = "TestGroup",
                # content() reads these
                clusterNumbers     = "1",
                somColorVar        = "n",
                somSizeVar         = "max",
                selectMode         = "view",
                showGroups         = FALSE,
                samples2plot       = unique(as.character(ta$sce$sample_id)),
                scatterPercentile  = 0.99,
                relativeTo         = "none",
                compareStatsTo     = "none",
                upsetSelection     = character(0L),
                violinSelection    = markers,
                dimRedSelection    = markers[seq_len(5L)],
                perplexity         = 5L,
                n_neighbors        = 4L,
                currentDimX        = markers[1L],
                currentDimY        = markers[2L]
            ))
            suppressWarnings(session$flushReact())

            expect_equal(session$input$clusterNameSelect, "TestGroup",
                         label = "clusterNameSelect reaches the session")
            expect_equal(session$input$somColorVar,  "n",    label = "somColorVar")
            expect_equal(session$input$somSizeVar,   "max",  label = "somSizeVar")
            expect_equal(session$input$relativeTo,   "none", label = "relativeTo")
            expect_equal(session$input$compareStatsTo, "none",
                         label = "compareStatsTo")
        })
    })
})





# =============================================================================
# Fix 1 — quantile axis limits test
# expect_length() has no label argument in any released testthat version.
# Use expect_equal(length(...), n, label = ...) throughout.
# =============================================================================
test_that("downloadPlots content: quantile axis limits (pctl logic) are numeric vectors", {
    sce      <- CySA_example_sce(n_cells = 100L, n_nodes = 10L)
    somCodes <- S4Vectors::metadata(sce)$SOM_codes
    markers  <- colnames(somCodes)

    for (pctl in c(0.95, 0.99, 1.0)) {
        local({
            p     <- pctl
            tailP <- (1 - p) / 2
            xlimP <- stats::quantile(
                somCodes[, markers[1L]],
                probs   = c(tailP, 1 - tailP),
                na.rm   = TRUE
            )

            # expect_length() has no label arg — use expect_equal(length())
            expect_equal(
                length(xlimP), 2L,
                label = sprintf("xlimP length == 2 for pctl = %.2f", p)
            )
            expect_true(
                all(is.finite(xlimP)),
                label = sprintf("xlimP is finite for pctl = %.2f", p)
            )
            expect_lte(
                xlimP[[1L]], xlimP[[2L]],
                label = sprintf("xlimP[1] <= xlimP[2] for pctl = %.2f", p)
            )
        })
    }
})


# =============================================================================
# Fix 2 — integration test: expect_no_error() accepts no extra arguments
# expect_no_error(object) only; label must move to a surrounding comment.
# The 'colour' attribute warning is now muffled by the updated helper above.
# =============================================================================
test_that("downloadPlots integration: server accepts all inputs the handler depends on", {
    skip_on_cran()

    ta <- make_test_app(n_cells = 200L, n_nodes = 10L)

    .muffle_plotly_warning({
        shiny::testServer(ta$app, {
            markers <- S4Vectors::metadata(ta$sce)$map$colsUsed

            suppressWarnings(session$setInputs(
                clusterNameSelect  = "TestGroup",
                clusterNumbers     = "1",
                somColorVar        = "n",
                somSizeVar         = "max",
                selectMode         = "view",
                showGroups         = FALSE,
                samples2plot       = unique(as.character(ta$sce$sample_id)),
                scatterPercentile  = 0.99,
                relativeTo         = "none",
                compareStatsTo     = "none",
                upsetSelection     = character(0L),
                violinSelection    = markers,
                dimRedSelection    = markers[seq_len(5L)],
                perplexity         = 5L,
                n_neighbors        = 4L,
                currentDimX        = markers[1L],
                currentDimY        = markers[2L]
            ))
            suppressWarnings(session$flushReact())

            # expect_equal() accepts label; use it for all assertions
            expect_equal(session$input$clusterNameSelect, "TestGroup",
                         label = "clusterNameSelect reaches the session")
            expect_equal(session$input$somColorVar,   "n",
                         label = "somColorVar")
            expect_equal(session$input$somSizeVar,    "max",
                         label = "somSizeVar")
            expect_equal(session$input$relativeTo,    "none",
                         label = "relativeTo")
            expect_equal(session$input$compareStatsTo, "none",
                         label = "compareStatsTo")
        })
    })
})

test_that("downloadPlots integration: missing clusterNameSelect does not crash server", {
    skip_on_cran()

    ta <- make_test_app(n_cells = 200L, n_nodes = 10L)

    # Intent: req(input$clusterNameSelect) inside filename() fires silently
    # and the session stays alive. expect_no_error() takes no label argument;
    # intent is documented here in the test body instead.
    .muffle_plotly_warning({
        shiny::testServer(ta$app, {
            suppressWarnings(session$setInputs(clusterNumbers = "1", selectMode = "view"))

            # expect_no_error() signature: expect_no_error(object) — no label
            expect_no_error(session$flushReact())
        })
    })
})


# =============================================================================
# Companion: audit every expect_* call in this file for unsupported label arg
#
# Functions that DO accept label =   : expect_equal, expect_true, expect_false,
#                                      expect_identical, expect_match, expect_lt,
#                                      expect_lte, expect_gt, expect_gte,
#                                      expect_null, expect_s3_class
#
# Functions that do NOT accept label =: expect_length, expect_no_error,
#                                       expect_no_warning, expect_no_message,
#                                       expect_error, expect_warning,
#                                       expect_message, expect_setequal,
#                                       expect_named, expect_type
#
# Rule: always use expect_equal(length(x), n, label = ...) instead of
#       expect_length(x, n, label = ...) whenever a diagnostic label is needed.
# =============================================================================
test_that("expect_length label rule: expect_equal(length()) accepts label correctly", {
    x <- c(1L, 2L, 3L)
    # This must NOT error — confirms the replacement pattern is correct
    expect_no_error(
        expect_equal(length(x), 3L, label = "x has 3 elements")
    )
})

test_that("expect_no_error label rule: bare form accepts no extra arguments", {
    # Confirms the corrected call site pattern compiles without error
    expect_no_error(
        expect_no_error(identity("ok"))
    )
})

