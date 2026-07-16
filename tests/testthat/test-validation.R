test_that(".initializeOutputList() adds Rest and removes empties", {
    ol <- list(A = c(1, 2), B = c(3))
    levels <- c(1, 2, 3, 4, 5)
    res <- .initializeOutputList(ol, levels)

    expect_true("Rest" %in% names(res))
    expect_equal(sort(res$Rest), c(4L, 5L))
})

test_that(".computeRelativeCounts() returns percentages", {
    cpt <- matrix(
        c(10, 20, 30, 40),
        nrow = 2,
        dimnames = list(c("s1", "s2"), c("1", "2"))
    )
    class(cpt) <- "table"

    res <- .computeRelativeCounts(cpt, rs = c("1", "2"), relativeToCol = "none",
                                  expInfo = data.frame(), outputList = list())
    expect_type(res, "double")
    expect_equal(unname(res), c(100, 100))
})

test_that(".validateClusterSelectorInputs() checks required inputs", {
    sce <- CySA_example_sce()
    sce_bad <- sce
    S4Vectors::metadata(sce_bad)$map <- NULL

    expect_error(
        clusterSelector(
            sce = sce_bad,
            sce_subsampled = sce_bad,
            dList = list(c("marker1", "marker2")),
            dend = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
            dendTable = data.frame(id = 1),
            clusterPatientTable = table(c(1)),
            somRasterData = data.frame(x = 1, y = 1, marker1 = 1)
        ),
        "colsUsed is required"
    )
})

test_that(".validateClusterSelectorInputs() requires cluster_id and sample_id", {
    sce <- CySA_example_sce()
    sce_bad <- sce
    SummarizedExperiment::colData(sce_bad)$cluster_id <- NULL

    expect_error(
        clusterSelector(
            sce = sce_bad,
            sce_subsampled = sce_bad,
            dList = list(c("marker1", "marker2")),
            dend = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
            dendTable = data.frame(id = 1),
            clusterPatientTable = table(c(1)),
            somRasterData = data.frame(x = 1, y = 1, marker1 = 1)
        ),
        "cluster_id.*colData"
    )
})

test_that(".validateClusterSelectorInputs() requires dList with enough pairs", {
    sce <- CySA_example_sce()
    expect_error(
        clusterSelector(
            sce = sce,
            sce_subsampled = sce,
            dList = list(c("marker1", "marker2")),
            nPlots = 6,
            dend = stats::as.dendrogram(stats::hclust(stats::dist(matrix(1:4, nrow = 2)))),
            dendTable = data.frame(id = 1),
            clusterPatientTable = table(c(1)),
            somRasterData = data.frame(x = 1, y = 1, marker1 = 1)
        ),
        "dList.*at least"
    )
})

test_that(".buildSOMCodes() adds missing markers to SOM_codes", {
    sce <- CySA_example_sce()
    S4Vectors::metadata(sce)$SOM_codes <- S4Vectors::metadata(sce)$SOM_codes[, 1:6, drop = FALSE]
    sce <- .buildSOMCodes(sce)

    expect_true(all(rownames(sce) %in% colnames(S4Vectors::metadata(sce)$SOM_codes)))
})
