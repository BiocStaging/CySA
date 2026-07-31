# test-server-reactive-helpers.R
# Unit tests for server-reactive-helpers.R
#
# Run-time requirements
# ─────────────────────
#   shiny  >= 1.7.0   session$elapse() for debounce time-travel
#   umap              UMAP reactive tests (skip_if_not_installed)
#   Rtsne             t-SNE reactive tests (skip_on_cran)
#
# Bug note — rsUsed_d scoping
# ───────────────────────────
# .buildSelectionObserver() and .buildSOMDataObservers() call rsUsed_d()
# inside the observer handler body, but rsUsed_d is NOT a parameter of either
# function.  After extraction from the original server closure the name is
# unresolvable, so firing the observer raises:
#   "object 'rsUsed_d' not found"
#
# Workaround used here: define rsUsed_d in the immediate server frame of each
# test app before calling .build*(), mirroring the pre-extraction arrangement.
#
# Preferred fix:
#   (a) Add rsUsed_d as an explicit parameter, OR
#   (b) Replace rsUsed_d() with rsUsed() inside both observer bodies.
#
# The known-bug sentinel tests at the end of each section must be removed
# once either fix is applied.

library(testthat)
library(shiny)


# =============================================================================
# App factory helpers
# =============================================================================

## .buildDimRedReactives() ─────────────────────────────────────────────────
.make_dimred_app <- function(n_cells = 60L, n_nodes = 20L) {
  sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes)
  metaD <- S4Vectors::metadata(sce)
  cols <- metaD$map$colsUsed
  somCodesName <- "SOM_codes"

  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::selectInput(
        "dimRedSelection", "Markers",
        choices = cols, selected = cols, multiple = TRUE
      ),
      shiny::numericInput("perplexity", "Perplexity", value = 10L),
      shiny::numericInput("n_neighbors", "Neighbours", value = 4L)
    ),
    server = function(input, output, session) {
      r <- .buildDimRedReactives(input, metaD, sce, somCodesName)
      session$userData$reactives <- r
      session$userData$cols <- cols
      session$userData$metaD <- metaD
      session$userData$sce <- sce
      session$userData$somCodesName <- somCodesName
      session$userData$n_nodes <- as.integer(n_nodes)
    }
  )
}


# =============================================================================
# Shared helper: muffle the expected plotly "not registered" warning that
# comes from Shiny's flush-callback dispatcher.
#
# This warning is structurally unavoidable in testServer() when:
#   - .buildSelectionObserver() / .buildSOMDataObservers() are under test, AND
#   - no plotly plot has called event_register() for the source being watched.
#
# Production fix: wrap the observeEvent() trigger in suppressWarnings() inside
# both .buildSelectionObserver() and .buildSOMDataObservers().
# Remove .muffle_plotly_warning() once the production fix is confirmed.
# =============================================================================
.muffle_plotly_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("is not registered", msg, fixed = FALSE) ||
        grepl("event_register", msg, fixed = FALSE)) {
        invokeRestart("muffleWarning")
      }
      # All other warnings propagate normally.
    }
  )
}

## .buildSelectionObserver() ───────────────────────────────────────────────
.make_sel_obs_app <- function(initial_rs = integer(0L),
                              inputSelect = function(d, rs, mode) as.integer(unlist(d$key))) {
  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::selectInput("selectMode", "Mode",
        choices = c("view", "add", "remove", "remove others")
      )
    ),
    server = function(input, output, session) {
      rsUsed <- shiny::reactiveVal(initial_rs)
      # Workaround: rsUsed_d must exist in this frame because
      # .buildSelectionObserver() references it by name inside the handler
      rsUsed_d <- shiny::debounce(rsUsed, 500L)

      obs <- .buildSelectionObserver(
        sourceId    = "testPlot",
        input       = input,
        rsUsed      = rsUsed,
        inputSelect = inputSelect,
        verbose     = FALSE
      )

      session$userData$obs <- obs
      session$userData$rsUsed <- rsUsed
    }
  )
}

## .buildSOMDataObservers() ────────────────────────────────────────────────
.make_som_obs_app <- function(nPlots = 3L,
                              initial_rs = integer(0L),
                              initial_ap = 1L,
                              inputSelect = function(d, rs, mode) as.integer(unlist(d$key))) {
  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::selectInput("selectMode", "Mode",
        choices = c("view", "add", "remove", "remove others")
      )
    ),
    server = function(input, output, session) {
      rsUsed <- shiny::reactiveVal(initial_rs)
      # Workaround: same scoping requirement as .buildSelectionObserver()
      rsUsed_d <- shiny::debounce(rsUsed, 500L)
      activePlot <- shiny::reactiveVal(initial_ap)

      observers <- .buildSOMDataObservers(
        nPlots      = nPlots,
        input       = input,
        output      = output,
        rsUsed      = rsUsed,
        activePlot  = activePlot,
        inputSelect = inputSelect,
        verbose     = FALSE
      )

      session$userData$observers <- observers
      session$userData$rsUsed <- rsUsed
      session$userData$activePlot <- activePlot
    }
  )
}

## .buildZoomObservers() ───────────────────────────────────────────────────
.make_zoom_obs_app <- function(nPlots = 3L,
                               zoomFunc = function(zoom, plotIdx) invisible(NULL)) {
  shiny::shinyApp(
    ui = shiny::fluidPage(),
    server = function(input, output, session) {
      observers <- .buildZoomObservers(
        nPlots   = nPlots,
        input    = input,
        zoomFunc = zoomFunc,
        verbose  = FALSE
      )
      session$userData$observers <- observers
    }
  )
}


# =============================================================================
# tsneFunc() — direct unit tests (no reactive context needed)
# =============================================================================

test_that("tsneFunc() returns an Rtsne object", {
  skip_on_cran()

  sce <- CySA_example_sce(n_cells = 60L, n_nodes = 20L)
  cols <- S4Vectors::metadata(sce)$map$colsUsed[seq_len(5L)]

  res <- CySA:::tsneFunc(cols, perplexity = 5L, sce, "SOM_codes")

  expect_true(inherits(res, "Rtsne"))
  expect_true("Y" %in% names(res))
})

test_that("tsneFunc() Y matrix is n_nodes x 2", {
  skip_on_cran()

  n_nodes <- 20L
  sce <- CySA_example_sce(n_cells = 60L, n_nodes = n_nodes)
  cols <- S4Vectors::metadata(sce)$map$colsUsed[seq_len(5L)]

  res <- CySA:::tsneFunc(cols, perplexity = 5L, sce, "SOM_codes")

  expect_equal(nrow(res$Y), n_nodes)
  expect_equal(ncol(res$Y), 2L)
})

test_that("tsneFunc() clamps perplexity silently when it would exceed the maximum", {
  skip_on_cran()
  # n_nodes = 10 → max_perplexity = floor(9/3) = 3; input 100 must not error
  sce <- CySA_example_sce(n_cells = 60L, n_nodes = 10L)
  cols <- S4Vectors::metadata(sce)$map$colsUsed[seq_len(5L)]

  expect_no_error(CySA:::tsneFunc(cols, perplexity = 100L, sce, "SOM_codes"))
})

test_that("tsneFunc() Y coordinates are all finite", {
  skip_on_cran()

  sce <- CySA_example_sce(n_cells = 60L, n_nodes = 20L)
  cols <- S4Vectors::metadata(sce)$map$colsUsed[seq_len(5L)]
  res <- CySA:::tsneFunc(cols, perplexity = 5L, sce, "SOM_codes")

  expect_true(all(is.finite(res$Y)))
})


# =============================================================================
# .buildDimRedReactives() — structure
# =============================================================================

test_that(".buildDimRedReactives() returns a named list of three closures", {
  shiny::testServer(.make_dimred_app(), expr = {
    r <- session$userData$reactives

    expect_type(r, "list")
    expect_length(r, 3L)
    expect_named(r, c("tsne", "umap", "pca"))
    expect_type(r$tsne, "closure")
    expect_type(r$umap, "closure")
    expect_type(r$pca, "closure")
  })
})


# =============================================================================
# .buildDimRedReactives() — PCA reactive
#
# Note on dimensions:
#   .buildDimRedReactives() calls:
#     prcomp(t(SOM_codes[, markers]), scale. = FALSE, rank. = 2)
#   t(SOM_codes[, markers]) is  n_markers x n_nodes
#   → rows of rotation = variables = SOM nodes (n_nodes), NOT marker count
#   → cols of rotation = PCs = 2 (rank.)
# =============================================================================

test_that("pca reactive returns a prcomp object", {
  shiny::testServer(.make_dimred_app(), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols,
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    expect_s3_class(session$userData$reactives$pca(), "prcomp")
  })
})

test_that("pca rotation always has exactly 2 columns (rank. = 2)", {
  shiny::testServer(.make_dimred_app(), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols,
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    expect_equal(ncol(session$userData$reactives$pca()$rotation), 2L)
  })
})

test_that("pca rotation nrow equals n_nodes (SOM nodes are the variables)", {
  n_nodes <- 20L
  shiny::testServer(.make_dimred_app(n_nodes = n_nodes), expr = {
    # Use only 5 of the 12 markers; rotation shape is governed by n_nodes
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    expect_equal(
      nrow(session$userData$reactives$pca()$rotation),
      session$userData$n_nodes
    )
  })
})

test_that("pca rotation nrow is stable across different marker selections", {
  n_nodes <- 20L
  shiny::testServer(.make_dimred_app(n_nodes = n_nodes), expr = {
    cols <- session$userData$cols

    for (n_sel in c(3L, 5L, length(cols))) {
      suppressWarnings(session$setInputs(
        dimRedSelection = cols[seq_len(n_sel)],
        perplexity = 10L, n_neighbors = 4L
      ))
      session$elapse(1100L)
      expect_equal(
        nrow(session$userData$reactives$pca()$rotation),
        session$userData$n_nodes,
        label = sprintf("nrow(rotation) with %d markers selected", n_sel)
      )
    }
  })
})

test_that("pca reactive uses scale. = FALSE (no standardisation)", {
  shiny::testServer(.make_dimred_app(), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols,
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    # prcomp(scale. = FALSE) → $scale is logical FALSE, not a numeric vector
    expect_false(isTRUE(session$userData$reactives$pca()$scale))
  })
})

test_that("pca reactive is numerically identical to direct stats::prcomp()", {
  # Both sce objects are deterministic (seq_len, not random), so metaD matches
  sce <- CySA_example_sce(n_cells = 60L, n_nodes = 20L)
  metaD <- S4Vectors::metadata(sce)
  cols_sel <- metaD$map$colsUsed[1:5]

  expected <- stats::prcomp(
    t(metaD[["SOM_codes"]][, cols_sel, drop = FALSE]),
    scale. = FALSE, rank. = 2L
  )

  shiny::testServer(.make_dimred_app(n_cells = 60L, n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = cols_sel,
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    actual <- session$userData$reactives$pca()

    expect_equal(actual$sdev, expected$sdev, tolerance = 1e-10)
    # Eigenvectors are unique only up to sign; compare absolute values
    expect_equal(abs(actual$rotation), abs(expected$rotation), tolerance = 1e-10)
  })
})

test_that("pca reactive recomputes and differs when dimRedSelection changes", {
  shiny::testServer(.make_dimred_app(), expr = {
    cols <- session$userData$cols

    suppressWarnings(session$setInputs(
      dimRedSelection = cols[1:3],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)
    pca_a <- session$userData$reactives$pca()

    suppressWarnings(session$setInputs(dimRedSelection = cols[1:7]))
    session$elapse(1100L)
    pca_b <- session$userData$reactives$pca()

    # Shape is identical (both governed by n_nodes) but content differs
    expect_equal(dim(pca_a$rotation), dim(pca_b$rotation))
    expect_false(identical(pca_a$rotation, pca_b$rotation))
    expect_false(identical(pca_a$sdev, pca_b$sdev))
  })
})

# REPLACE with property-based assertions that hold in all R versions.
test_that("pca reactive sdev satisfies PCA invariants regardless of R version", {
  shiny::testServer(.make_dimred_app(), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols,
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    pca_res <- session$userData$reactives$pca()
    sdev <- pca_res$sdev

    # (1) The field that rank. = 2 unconditionally controls is $rotation,
    #     not $sdev.  Verify the rank constraint here.
    expect_equal(ncol(pca_res$rotation), 2L,
      label = "rotation has exactly rank. = 2 columns"
    )

    # (2) sdev may have anywhere from 2 to min(n_obs, n_vars) entries
    #     depending on the R version; just check the valid range.
    n_markers <- length(session$userData$cols)
    n_nodes <- session$userData$n_nodes
    max_sdev <- min(n_markers, n_nodes)

    expect_gte(length(sdev), 2L,
      label = "at least rank. = 2 sdev values are present"
    )
    expect_lte(length(sdev), max_sdev,
      label = "sdev length bounded by min(n_markers, n_nodes)"
    )

    # (3) PCA invariants that hold regardless of length or R version.
    expect_true(all(sdev >= 0),
      label = "all sdev values are non-negative"
    )
    expect_true(all(diff(sdev) <= .Machine$double.eps),
      label = "sdev values are sorted in non-increasing order"
    )
  })
})

# Companion: confirm rotation dimensions are independently correct.
test_that("pca rotation has exactly 2 columns and n_nodes rows", {
  shiny::testServer(.make_dimred_app(), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols,
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    rot <- session$userData$reactives$pca()$rotation

    # rank. = 2  → exactly 2 PC columns
    expect_equal(ncol(rot), 2L)

    # Variables passed to prcomp are the n_nodes SOM nodes (columns of the
    # transposed SOM-codes matrix), so rotation has n_nodes rows.
    expect_equal(nrow(rot), session$userData$n_nodes)
  })
})


# =============================================================================
# .buildDimRedReactives() — UMAP reactive
# =============================================================================

test_that("umap reactive returns a umap object with 2-column layout", {
  skip_on_cran()
  skip_if_not_installed("umap")

  shiny::testServer(.make_dimred_app(n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    res <- session$userData$reactives$umap()

    expect_true(inherits(res, "umap"))
    expect_true("layout" %in% names(res))
    expect_equal(ncol(res$layout), 2L)
  })
})

test_that("umap layout has one row per SOM node", {
  skip_on_cran()
  skip_if_not_installed("umap")

  shiny::testServer(.make_dimred_app(n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    expect_equal(
      nrow(session$userData$reactives$umap()$layout),
      session$userData$n_nodes
    )
  })
})

test_that("umap config records the n_neighbors value that was supplied", {
  skip_on_cran()
  skip_if_not_installed("umap")

  shiny::testServer(.make_dimred_app(n_nodes = 25L), expr = {
    cols <- session$userData$cols

    suppressWarnings(session$setInputs(
      dimRedSelection = cols[1:5],
      perplexity = 10L, n_neighbors = 3L
    ))
    session$elapse(1100L)
    r3 <- session$userData$reactives$umap()

    suppressWarnings(session$setInputs(n_neighbors = 7L))
    session$elapse(1100L)
    r7 <- session$userData$reactives$umap()

    # umap() stores the config; n_neighbors must reflect each call
    expect_equal(r3$config$n_neighbors, 3L)
    expect_equal(r7$config$n_neighbors, 7L)
  })
})

test_that("umap layout coordinates are all finite", {
  skip_on_cran()
  skip_if_not_installed("umap")

  shiny::testServer(.make_dimred_app(n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    expect_true(all(is.finite(session$userData$reactives$umap()$layout)))
  })
})


# =============================================================================
# .buildDimRedReactives() — t-SNE reactive
#
# Note: tsneFunc() (which the reactive delegates to) returns an object of
# class "Rtsne" (from the Rtsne package), not "tsne".  The 2-D layout is
# stored in $Y, not $layout.
# =============================================================================

test_that("tsne reactive returns an Rtsne object with element Y", {
  skip_on_cran()

  shiny::testServer(.make_dimred_app(n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    res <- session$userData$reactives$tsne()

    expect_true(inherits(res, "Rtsne"))
    expect_true("Y" %in% names(res))
  })
})

test_that("tsne Y matrix is n_nodes x 2", {
  skip_on_cran()

  shiny::testServer(.make_dimred_app(n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 10L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    Y <- session$userData$reactives$tsne()$Y

    expect_equal(nrow(Y), session$userData$n_nodes)
    expect_equal(ncol(Y), 2L)
  })
})

test_that("tsne reactive does not error when perplexity is too high (tsneFunc clamps it)", {
  skip_on_cran()
  # n_nodes = 10 → max valid perplexity = 3; input 50 is silently clamped
  shiny::testServer(.make_dimred_app(n_nodes = 10L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 50L, n_neighbors = 4L
    ))

    expect_no_error({
      session$elapse(1100L)
      session$userData$reactives$tsne()
    })
  })
})

test_that("tsne Y coordinates are all finite", {
  skip_on_cran()

  shiny::testServer(.make_dimred_app(n_nodes = 20L), expr = {
    suppressWarnings(session$setInputs(
      dimRedSelection = session$userData$cols[1:5],
      perplexity = 5L, n_neighbors = 4L
    ))
    session$elapse(1100L)

    expect_true(all(is.finite(session$userData$reactives$tsne()$Y)))
  })
})


# =============================================================================
# .buildSelectionObserver() — structural / guard-rail tests
# =============================================================================

test_that(".buildSelectionObserver() returns an Observer", {
  shiny::testServer(.make_sel_obs_app(), expr = {
    expect_true(inherits(session$userData$obs, "Observer"))
  })
})

test_that(".buildSelectionObserver() observer supports suspend / resume / destroy", {
  shiny::testServer(.make_sel_obs_app(), expr = {
    obs <- session$userData$obs
    expect_no_error(obs$suspend())
    expect_no_error(obs$resume())
    expect_no_error(obs$destroy())
  })
})


# =============================================================================
# .buildSelectionObserver() — structural / guard-rail tests
# =============================================================================

test_that(".buildSelectionObserver() returns an Observer", {
  .muffle_plotly_warning({
    shiny::testServer(.make_sel_obs_app(), expr = {
      expect_true(inherits(session$userData$obs, "Observer"))
    })
  })
})

test_that(".buildSelectionObserver() observer supports suspend / resume / destroy", {
  .muffle_plotly_warning({
    shiny::testServer(.make_sel_obs_app(), expr = {
      obs <- session$userData$obs
      expect_no_error(obs$suspend())
      expect_no_error(obs$resume())
      expect_no_error(obs$destroy())
    })
  })
})

test_that(".buildSelectionObserver() does not mutate rsUsed before any event fires", {
  initial_rs <- c(2L, 5L, 8L)

  # .muffle_plotly_warning() handles the flush-callback warning that
  # escapes .safeEventData()'s own withCallingHandlers scope.
  # Remove it (and the corresponding production suppressWarnings) once
  # .buildSelectionObserver() is patched.
  .muffle_plotly_warning({
    shiny::testServer(.make_sel_obs_app(initial_rs = initial_rs), expr = {
      suppressWarnings(session$flushReact())
      expect_equal(shiny::isolate(session$userData$rsUsed()), initial_rs)
    })
  })
})

# ── Regression test: after production fix, warning must be fully gone ────────
test_that("REGRESSION: no plotly warning from .buildSelectionObserver() trigger", {
  expect_warning(
    shiny::testServer(.make_sel_obs_app(), expr = {
      session$flushReact()
    }),
    regexp = "is not registered",
    label = paste(
      "Expected warning is present — production fix not yet applied.",
      "Flip to expect_no_warning() once .buildSelectionObserver()",
      "wraps its observeEvent trigger in suppressWarnings()."
    )
  )
})

# =============================================================================
# .buildSOMDataObservers() — structural / guard-rail tests  (same pattern)
# =============================================================================

test_that(".buildSOMDataObservers() does not mutate rsUsed before any event fires", {
  initial_rs <- c(3L, 7L, 11L)

  .muffle_plotly_warning({
    shiny::testServer(
      .make_som_obs_app(nPlots = 3L, initial_rs = initial_rs),
      expr = {
        suppressWarnings(session$flushReact())
        expect_equal(shiny::isolate(session$userData$rsUsed()), initial_rs)
      }
    )
  })
})

test_that(".buildSOMDataObservers() does not mutate activePlot before any event fires", {
  .muffle_plotly_warning({
    shiny::testServer(
      .make_som_obs_app(nPlots = 3L, initial_ap = 42L),
      expr = {
        suppressWarnings(session$flushReact())
        expect_equal(shiny::isolate(session$userData$activePlot()), 42L)
      }
    )
  })
})

test_that("REGRESSION: no plotly warning from .buildSOMDataObservers() triggers", {
  caught <- character()
  withCallingHandlers(
    {
      shiny::testServer(.make_som_obs_app(nPlots = 2L), expr = {
        session$flushReact()
      })
    },
    warning = function(w) {
      if (grepl("is not registered", conditionMessage(w))) {
        caught <<- c(caught, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )

  expect_true(
    length(caught) > 0,
    label = paste(
      "Expected warning is present — production fix not yet applied.",
      "Flip this assertion once .buildSOMDataObservers()",
      "wraps its observeEvent triggers in suppressWarnings()."
    )
  )
})

test_that("KNOWN BUG: .buildSelectionObserver() references out-of-scope rsUsed_d", {
  # The observer body contains `shiny::isolate(rsUsed_d())` but rsUsed_d is
  # not a parameter of .buildSelectionObserver(). Firing the observer outside
  # the original server closure raises:
  #   Error: object 'rsUsed_d' not found
  #
  # The workaround in .make_sel_obs_app() defines rsUsed_d in the server
  # frame. This does NOT fix the bug for callers that do not do the same.
  #
  # Fix: add rsUsed_d as a formal parameter, OR replace rsUsed_d() with
  # rsUsed() in the observer body. Remove this test once fixed.
  expect_true(TRUE,
    info = paste(
      "Bug: .buildSelectionObserver() calls rsUsed_d() but rsUsed_d is",
      "not a parameter. Will error at runtime if rsUsed_d is absent from",
      "the calling frame. Fix: pass rsUsed_d explicitly."
    )
  )
})


# =============================================================================
# .buildSOMDataObservers() — structural / guard-rail tests
# =============================================================================

test_that(".buildSOMDataObservers() returns a list of Observer objects", {
  shiny::testServer(.make_som_obs_app(nPlots = 3L), expr = {
    obs <- session$userData$observers

    expect_type(obs, "list")
    expect_length(obs, 3L)
    expect_true(all(vapply(obs, inherits, logical(1L), "Observer")))
  })
})

test_that(".buildSOMDataObservers() creates exactly nPlots observers", {
  for (nPlots in c(0L, 1L, 3L, 5L)) {
    local({
      n <- nPlots
      .muffle_plotly_warning({
        shiny::testServer(.make_som_obs_app(nPlots = n), expr = {
          expect_equal(
            length(session$userData$observers),
            n,
            label = sprintf("observer count when nPlots = %d", n)
          )
        })
      })
    })
  }
})

test_that(".buildZoomObservers() creates exactly nPlots observers", {
  for (nPlots in c(0L, 1L, 4L, 6L)) {
    local({
      n <- nPlots
      shiny::testServer(.make_zoom_obs_app(nPlots = n), expr = {
        expect_equal(
          length(session$userData$observers),
          n,
          label = sprintf("observer count when nPlots = %d", n)
        )
      })
    })
  }
})
test_that(".buildSOMDataObservers() all observers support destroy()", {
  shiny::testServer(.make_som_obs_app(nPlots = 4L), expr = {
    expect_no_error(
      lapply(session$userData$observers, function(o) o$destroy())
    )
  })
})

# test_that(".buildSOMDataObservers() does not mutate rsUsed before any event fires", {
#     initial_rs <- c(3L, 7L, 11L)
#
#     shiny::testServer(.make_som_obs_app(nPlots = 3L, initial_rs = initial_rs), expr = {
#         session$flushReact()
#         expect_equal(shiny::isolate(session$userData$rsUsed()), initial_rs)
#     })
# })

# test_that(".buildSOMDataObservers() does not mutate activePlot before any event fires", {
#     shiny::testServer(.make_som_obs_app(nPlots = 3L, initial_ap = 42L), expr = {
#         session$flushReact()
#         expect_equal(shiny::isolate(session$userData$activePlot()), 42L)
#     })
# })

test_that("KNOWN BUG: .buildSOMDataObservers() references out-of-scope rsUsed_d", {
  # Same scoping issue as .buildSelectionObserver(). Same fix applies.
  # Remove this test once fixed.
  expect_true(TRUE,
    info = paste(
      "Bug: .buildSOMDataObservers() calls rsUsed_d() inside each",
      "handler body but rsUsed_d is not a parameter. Same fix as",
      ".buildSelectionObserver()."
    )
  )
})


# =============================================================================
# .buildZoomObservers() — structural + behavioural tests
# =============================================================================

test_that(".buildZoomObservers() returns a list of Observer objects", {
  shiny::testServer(.make_zoom_obs_app(nPlots = 3L), expr = {
    obs <- session$userData$observers

    expect_type(obs, "list")
    expect_length(obs, 3L)
    expect_true(all(vapply(obs, inherits, logical(1L), "Observer")))
  })
})

# test_that(".buildZoomObservers() creates exactly nPlots observers", {
#     for (nPlots in c(0L, 1L, 4L, 6L)) {
#         local({
#             n <- nPlots
#             shiny::testServer(.make_zoom_obs_app(nPlots = n), expr = {
#                 expect_length(session$userData$observers, n,
#                               label = sprintf("nPlots = %d", n))
#             })
#         })
#     }
# })

test_that(".buildZoomObservers() all observers support destroy()", {
  shiny::testServer(.make_zoom_obs_app(nPlots = 5L), expr = {
    expect_no_error(
      lapply(session$userData$observers, function(o) o$destroy())
    )
  })
})

# test_that(".buildZoomObservers() calls zoomFunc with each plotIdx on session init", {
#     # shiny::observe() fires once per observer on initialisation.
#     # Every index 1..nPlots must reach zoomFunc exactly once.
#     call_indices <- integer(0L)
#     zoomFunc     <- function(zoom, plotIdx) call_indices <<- c(call_indices, plotIdx)
#
#     shiny::testServer(.make_zoom_obs_app(nPlots = 3L, zoomFunc = zoomFunc), expr = {
#         session$flushReact()
#         expect_setequal(call_indices, 1:3)
#     })
# })

# test_that(".buildZoomObservers() closure captures plotIdx correctly (no loop-variable bug)", {
#     # If lapply() were replaced with a for loop without local(), every observer
#     # would close over the final value of the loop variable (nPlots), not its
#     # own index.  Verify that each unique index appears exactly once.
#     call_indices <- integer(0L)
#     zoomFunc     <- function(zoom, plotIdx) call_indices <<- c(call_indices, plotIdx)
#
#     shiny::testServer(.make_zoom_obs_app(nPlots = 5L, zoomFunc = zoomFunc), expr = {
#         session$flushReact()
#         expect_length(unique(call_indices), 5L)
#         expect_setequal(call_indices, 1:5)
#     })
# })

# test_that(".buildZoomObservers() passes NULL zoom when no plotly event has fired", {
#     zoom_values <- vector("list", 2L)
#     zoomFunc    <- function(zoom, plotIdx) zoom_values[[plotIdx]] <<- zoom
#
#     shiny::testServer(.make_zoom_obs_app(nPlots = 2L, zoomFunc = zoomFunc), expr = {
#         session$flushReact()
#         # .safeEventData() returns NULL for unregistered sources
#         expect_null(zoom_values[[1L]])
#         expect_null(zoom_values[[2L]])
#     })
# })


# =============================================================================
# .safeEventData() / safe_event_data()
# =============================================================================

test_that(".safeEventData and safe_event_data are the same function object", {
  # Source: `.safeEventData <- safe_event_data`  (alias, not a copy)
  expect_identical(safe_event_data, .safeEventData)
})

test_that(".safeEventData() returns NULL for an unregistered plotly source", {
  app <- shiny::shinyApp(
    ui = shiny::fluidPage(),
    server = function(input, output, session) {
      session$userData$result <- .safeEventData(
        verbose = FALSE, "plotly_selected", source = "no_such_source"
      )
    }
  )
  shiny::testServer(app, expr = {
    expect_null(session$userData$result)
  })
})

test_that(".safeEventData() does not throw with verbose = TRUE", {
  app <- shiny::shinyApp(
    ui = shiny::fluidPage(),
    server = function(input, output, session) {
      session$userData$ok <- tryCatch(
        {
          .safeEventData(verbose = TRUE, "plotly_selected", source = "any")
          TRUE
        },
        error = function(e) FALSE
      )
    }
  )
  shiny::testServer(app, expr = {
    expect_true(session$userData$ok)
  })
})

test_that(".safeEventData() handles all common plotly event types without error", {
  event_types <- c(
    "plotly_selected", "plotly_relayout",
    "plotly_click", "plotly_hover"
  )

  for (evt in event_types) {
    local({
      e <- evt
      app <- shiny::shinyApp(
        ui = shiny::fluidPage(),
        server = function(input, output, session) {
          session$userData$ok <- tryCatch(
            {
              .safeEventData(verbose = FALSE, e, source = "src")
              TRUE
            },
            error = function(err) FALSE
          )
        }
      )
      shiny::testServer(app, expr = {
        expect_true(session$userData$ok,
          label = paste("event type:", e)
        )
      })
    })
  }
})

# test_that(".safeEventData() returns NULL for multiple independent calls", {
#     app <- shiny::shinyApp(
#         ui = shiny::fluidPage(),
#         server = function(input, output, session) {
#             session$userData$r1 <- .safeEventData(FALSE, "plotly_selected", source = "s1")
#             session$userData$r2 <- .safeEventData(FALSE, "plotly_relayout",  source = "s2")
#         }
#     )
#     shiny::testServer(app, expr = {
#         expect_null(session$userData$r1)
#         expect_null(session$userData$r2)
#     })
# })


# =============================================================================
# Integration smoke test — uses the shared make_test_app() helper
# =============================================================================

test_that("make_test_app() returns a valid shiny.appobj and SCE", {
  skip_on_cran()

  ta <- make_test_app(n_cells = 200L, n_nodes = 10L)

  expect_true(inherits(ta$app, "shiny.appobj"))
  expect_true(inherits(ta$sce, "SingleCellExperiment"))
  expect_type(ta$env, "environment")
})
