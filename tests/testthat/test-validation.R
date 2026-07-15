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

test_that(".buildFlowSOMPiePlot() returns NULL for invalid rs", {
    sce <- CySA_example_sce()
    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed

    expect_null(.buildFlowSOMPiePlot(som_codes, rs = integer(0), colsUsed))
    expect_null(.buildFlowSOMPiePlot(som_codes, rs = 9999, colsUsed))
})
