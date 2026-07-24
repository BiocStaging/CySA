# tests/testthat/test-CySA.R
library(testthat)
library(CySA)
library(SingleCellExperiment)
library(S4Vectors)
library(SummarizedExperiment)

# ── shared fixture ─────────────────────────────────────────────────────────────
small_sce <- function(n_cells = 300, n_nodes = 10) {
    CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes)
}

# ══════════════════════════════════════════════════════════════════════════════
# CySA_example_sce()
# ══════════════════════════════════════════════════════════════════════════════

test_that("CySA_example_sce() returns a SingleCellExperiment", {
    sce <- small_sce()
    expect_s4_class(sce, "SingleCellExperiment")
})

test_that("CySA_example_sce() colData has cluster_id and sample_id", {
    sce <- small_sce()
    expect_true("cluster_id" %in% names(colData(sce)))
    expect_true("sample_id"  %in% names(colData(sce)))
})

test_that("CySA_example_sce() metadata has SOM_codes, SOM_stats, map$colsUsed", {
    sce <- small_sce()
    md  <- S4Vectors::metadata(sce)
    expect_true("SOM_codes"  %in% names(md))
    expect_true("SOM_stats"  %in% names(md))
    expect_true("map"        %in% names(md))
    expect_true("colsUsed"   %in% names(md$map))
})

test_that("CySA_example_sce() n_nodes controls SOM grid rows", {
    sce <- small_sce(n_nodes = 8)
    expect_equal(nrow(S4Vectors::metadata(sce)$SOM_codes), 8L)
})

test_that("CySA_example_sce() has at least one assay", {
    sce <- small_sce()
    expect_gte(length(SummarizedExperiment::assays(sce)), 1L)
})

# ══════════════════════════════════════════════════════════════════════════════
# CySA_default_cluster_cols()
# ══════════════════════════════════════════════════════════════════════════════

test_that("CySA_default_cluster_cols() returns 20 valid hex colours", {
    cols <- CySA_default_cluster_cols()
    expect_length(cols, 20L)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
})

# ══════════════════════════════════════════════════════════════════════════════
# prepClusterSelectorData()
# ══════════════════════════════════════════════════════════════════════════════

test_that("prepClusterSelectorData() returns a named list", {
    sce    <- small_sce()
    result <- prepClusterSelectorData(sce, total_cells_to_sample = 80)
    expect_named(result, c("sce", "sce_subsampled", "dList"), ignore.order = TRUE)
})

test_that("prepClusterSelectorData() sce_subsampled is a SingleCellExperiment", {
    sce    <- small_sce()
    result <- prepClusterSelectorData(sce, total_cells_to_sample = 80)
    expect_s4_class(result$sce_subsampled, "SingleCellExperiment")
})

test_that("prepClusterSelectorData() subsamples to fewer cells", {
    sce    <- small_sce(n_cells = 400)
    result <- prepClusterSelectorData(sce, total_cells_to_sample = 100)
    expect_lte(ncol(result$sce_subsampled), ncol(result$sce))
})

test_that("prepClusterSelectorData() dList contains pairs of markers", {
    sce    <- small_sce()
    result <- prepClusterSelectorData(sce, total_cells_to_sample = 80)
    dl     <- result$dList
    expect_type(dl, "list")
    expect_gte(length(dl), 1L)
    expect_true(all(lengths(dl) == 2L))
})

test_that("prepClusterSelectorData() is reproducible with the same seed", {
    sce <- small_sce(n_cells = 400)
    r1  <- prepClusterSelectorData(sce, total_cells_to_sample = 100, seed = 7)
    r2  <- prepClusterSelectorData(sce, total_cells_to_sample = 100, seed = 7)
    expect_equal(colnames(r1$sce_subsampled), colnames(r2$sce_subsampled))
})

# ══════════════════════════════════════════════════════════════════════════════
# plotSOMScatter()
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotSOMScatter() returns a ggplot for valid channels", {
    sce  <- small_sce()
    chs  <- colnames(S4Vectors::metadata(sce)$SOM_codes)[seq_len(2)]
    p    <- plotSOMScatter(sce, chs = chs)
    expect_s3_class(p, "gg")
})

test_that("plotSOMScatter() errors when metaSlot is absent", {
    sce <- small_sce()
    expect_error(
        plotSOMScatter(sce, chs = c("marker1", "marker2"),
                       metaSlot = "ABSENT_SLOT"),
        regexp = "Need 'ABSENT_SLOT' in metadata"
    )
})

test_that("plotSOMScatter() warns for unknown channels", {
    sce <- small_sce()
    expect_warning(
        p <- plotSOMScatter(sce, chs = c("NO_SUCH_MARKER", "marker1", "marker2")),
        regexp = "Unknown channels/markers dropped"
    )
    expect_s3_class(p, "gg")
})

test_that("plotSOMScatter() returns NULL when < 2 valid channels remain", {
    sce    <- small_sce()
    result <- suppressWarnings(
        plotSOMScatter(sce, chs = c("MISS1", "MISS2"))
    )
    expect_null(result)
})

# ══════════════════════════════════════════════════════════════════════════════
# plotCytoScatter()
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotCytoScatter() returns a ggplot for valid channels", {
    skip_if_not_installed("KernSmooth")
    sce <- small_sce()
    chs <- rownames(sce)[seq_len(2)]
    p   <- plotCytoScatter(sce, chs = chs)
    expect_s3_class(p, "gg")
})

test_that("plotCytoScatter() errors on unknown channels", {
    sce <- small_sce()
    expect_error(
        plotCytoScatter(sce, chs = c("NO_SUCH_MARKER", "marker1")),
        regexp = "Unknown channels"
    )
})

test_that("plotCytoScatter() warns when gate is non-NULL", {
    skip_if_not_installed("KernSmooth")
    sce <- small_sce()
    chs <- rownames(sce)[seq_len(2)]
    expect_warning(
        plotCytoScatter(sce, chs = chs, gate = list()),
        regexp = "'gate' is not currently implemented"
    )
})

test_that("plotCytoScatter() works when color_by is NULL (density mode)", {
    skip_if_not_installed("KernSmooth")
    sce <- small_sce()
    chs <- rownames(sce)[seq_len(2)]
    p   <- plotCytoScatter(sce, chs = chs, color_by = NULL)
    expect_s3_class(p, "gg")
})

# ══════════════════════════════════════════════════════════════════════════════
# Internal helpers (via :::)
# ══════════════════════════════════════════════════════════════════════════════

test_that(".initializeOutputList() assigns unassigned clusters to Rest", {
    ol     <- list(Group1 = c(1L, 2L))
    result <- CySA:::.initializeOutputList(ol, clusterLevels = as.character(1:5))
    expect_true("Rest" %in% names(result))
    expect_setequal(result$Rest, 3:5)
    expect_false(any(result$Rest %in% c(1L, 2L)))
})

test_that(".initializeOutputList() removes empty groups", {
    ol     <- list(Group1 = c(1L, 2L), EmptyGrp = integer(0))
    result <- CySA:::.initializeOutputList(ol, clusterLevels = as.character(1:5))
    expect_false("EmptyGrp" %in% names(result))
})

test_that(".rebuildOutputList() recalculates Rest after removal", {
    ol     <- list(Group1 = c(1L, 2L), Group2 = c(3L, 4L))
    result <- CySA:::.rebuildOutputList(ol, sceLevels = as.character(1:6))
    expect_true("Rest" %in% names(result))
    expect_setequal(result$Rest, 5:6)
})

test_that(".computeScatterLimits() returns named xlim/ylim lists", {
    set.seed(1)
    lims <- CySA:::.computeScatterLimits(rnorm(200), rnorm(200), pctl = 0.95)
    expect_named(lims, c("xlim", "ylim"))
    expect_length(lims$xlim, 2L)
    expect_length(lims$ylim, 2L)
    expect_lt(lims$xlim[1], lims$xlim[2])
})

test_that(".computeRelativeCounts() returns vector of percentages (none mode)", {
    cpt <- matrix(c(10L, 20L, 5L, 15L), nrow = 2,
                  dimnames = list(c("s1", "s2"), c("1", "2")))
    cpt <- as.table(cpt)
    pct <- CySA:::.computeRelativeCounts(cpt, rs = "1",
                                         relativeToCol = "none",
                                         expInfo = NULL,
                                         outputList = list())
    # s1: 10/(10+20) * 100 = 33.3; s2: 20/(20+15) * 100 = 57.1
    expect_true(all(pct >= 0 & pct <= 100))
})

test_that(".inputSelect() applies 'add' mode correctly", {
    d  <- data.frame(customdata = c(3L, 5L), curveNumber = 0)
    rs <- c(1L, 2L, 3L)
    expect_setequal(CySA:::.inputSelect(d, rs, "add"), c(1L, 2L, 3L, 5L))
})

test_that(".inputSelect() applies 'remove' mode correctly", {
    d  <- data.frame(customdata = c(2L), curveNumber = 0)
    rs <- c(1L, 2L, 3L)
    expect_setequal(CySA:::.inputSelect(d, rs, "remove"), c(1L, 3L))
})

test_that(".inputSelect() returns rs unchanged for NULL customdata", {
    rs <- c(1L, 2L)
    expect_equal(CySA:::.inputSelect(NULL, rs, "add"), rs)
})

test_that(".buildFlowSOMPiePlot() returns a ggplot for valid rs", {
    sce      <- small_sce()
    codes    <- S4Vectors::metadata(sce)$SOM_codes
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed
    rs       <- seq_len(min(3L, nrow(codes)))
    p        <- CySA:::.buildFlowSOMPiePlot(codes, rs, colsUsed)
    expect_s3_class(p, "gg")
})

test_that(".buildFlowSOMPiePlot() returns NULL for empty rs", {
    sce      <- small_sce()
    codes    <- S4Vectors::metadata(sce)$SOM_codes
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed
    expect_null(CySA:::.buildFlowSOMPiePlot(codes, integer(0), colsUsed))
})

test_that(".buildBaseRasterGgplot() returns NULL for NULL input", {
    expect_null(CySA:::.buildBaseRasterGgplot(NULL, "marker1"))
})

test_that(".buildBaseRasterGgplot() returns a ggplot for valid data", {
    # Build a minimal somRasterData manually (no square-grid requirement)
    rd <- data.frame(
        x = rep(1:3, 2), y = rep(1:2, each = 3),
        marker1 = runif(6), marker2 = runif(6)
    )
    p <- CySA:::.buildBaseRasterGgplot(rd, c("marker1", "marker2"))
    expect_s3_class(p, "gg")
})

# ══════════════════════════════════════════════════════════════════════════════
# clusterSelector() — non-interactive structural test
# ══════════════════════════════════════════════════════════════════════════════

.make_app_inputs <- function() {
    sce     <- small_sce(n_cells = 300, n_nodes = 10)
    prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)

    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    dend      <- stats::as.dendrogram(
        stats::hclust(stats::dist(som_codes))
    )
    dendTable <- data.frame(
        id    = seq_len(nrow(som_codes)),
        label = rownames(som_codes),
        stringsAsFactors = FALSE
    )
    cpt <- table(
        sample_id  = sce$sample_id,
        cluster_id = sce$cluster_id
    )
    markers <- S4Vectors::metadata(sce)$map$colsUsed
    # Manual rasterData to avoid the perfect-square requirement
    n <- nrow(som_codes)
    rd <- data.frame(
        x  = rep(seq_len(5), length.out = n),
        y  = rep(seq_len(ceiling(n / 5)), each = 5)[seq_len(n)],
        id = seq_len(n)
    )
    for (m in markers) rd[[m]] <- seq_len(n) / n

    list(
        sce                = prepped$sce,
        sce_subsampled     = prepped$sce_subsampled,
        dList              = prepped$dList,
        dend               = dend,
        dendTable          = dendTable,
        clusterPatientTable = cpt,
        somRasterData      = rd
    )
}

test_that("clusterSelector() returns a shiny.appobj without error", {
    skip_if_not_installed("shiny")
    inp <- .make_app_inputs()
    app <- clusterSelector(
        sce                 = inp$sce,
        sce_subsampled      = inp$sce_subsampled,
        dList               = inp$dList,
        dend                = inp$dend,
        dendTable           = inp$dendTable,
        clusterPatientTable = inp$clusterPatientTable,
        somRasterData       = inp$somRasterData,
        somRasterObj        = NULL
    )
    expect_s3_class(app, "shiny.appobj")
})

test_that(".validateClusterSelectorInputs() errors when SOM_codes is missing", {
    inp  <- .make_app_inputs()
    sce2 <- inp$sce
    S4Vectors::metadata(sce2)$SOM_codes <- NULL

    expect_error(
        CySA:::.validateClusterSelectorInputs(
            sce2, inp$sce_subsampled, list(),
            inp$dList, inp$dend, inp$dendTable,
            inp$clusterPatientTable, inp$somRasterData,
            "SOM_codes", 6L, NULL
        ),
        regexp = "SOM_codes"
    )
})

test_that(".validateClusterSelectorInputs() errors when SOM_stats is missing", {
    inp  <- .make_app_inputs()
    sce2 <- inp$sce
    S4Vectors::metadata(sce2)$SOM_stats <- NULL

    expect_error(
        CySA:::.validateClusterSelectorInputs(
            sce2, inp$sce_subsampled, list(),
            inp$dList, inp$dend, inp$dendTable,
            inp$clusterPatientTable, inp$somRasterData,
            "SOM_codes", 6L, NULL
        ),
        regexp = "SOM_stats"
    )
})

test_that(".validateClusterSelectorInputs() errors when dList is too short", {
    inp <- .make_app_inputs()
    expect_error(
        CySA:::.validateClusterSelectorInputs(
            inp$sce, inp$sce_subsampled, list(),
            list(c("m1","m2")),   # only 1 pair, need >= nPlots = 6
            inp$dend, inp$dendTable,
            inp$clusterPatientTable, inp$somRasterData,
            "SOM_codes", 6L, NULL
        ),
        regexp = "at least 6"
    )
})
