test_that("plotSOMScatter() returns a ggplot object", {
    sce <- CySA_example_sce()
    p <- plotSOMScatter(sce, chs = c("marker1", "marker2"))

    expect_s3_class(p, "ggplot")
})

test_that("plotSOMScatter() handles missing channels gracefully", {
    sce <- CySA_example_sce()

    expect_warning(
        p <- plotSOMScatter(sce, chs = c("missing_marker", "marker1", "marker2")),
        "Unknown channels/markers"
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotSOMScatter() works with more than two channels", {
    sce <- CySA_example_sce()
    p <- plotSOMScatter(sce, chs = c("marker1", "marker2", "marker3", "marker4"))
    expect_s3_class(p, "ggplot")
})

test_that("plotCytoScatter() returns a ggplot object", {
    skip_if_not_installed("KernSmooth")
    sce <- CySA_example_sce()
    p <- plotCytoScatter(x = sce, chs = c("marker1", "marker2"))

    expect_s3_class(p, "ggplot")
})

test_that("plotCytoScatter() errors on unknown channels", {
    sce <- CySA_example_sce()
    expect_error(
        plotCytoScatter(x = sce, chs = c("not_a_marker", "marker1")),
        "Unknown channels/markers"
    )
})

test_that("plotCytoScatter() respects facet_by when supplied", {
    skip_if_not_installed("KernSmooth")
    sce <- CySA_example_sce()
    sce$group <- sample(c("A", "B"), ncol(sce), replace = TRUE)
    p <- plotCytoScatter(x = sce, chs = c("marker1", "marker2"), facet_by = "group")
    expect_s3_class(p, "ggplot")
})

test_that("plot-helpers build expected plot objects", {
    sce <- CySA_example_sce()
    metaD <- S4Vectors::metadata(sce)
    som_codes <- metaD$SOM_codes
    colsUsed <- metaD$map$colsUsed

    p_pie <- .buildFlowSOMPiePlot(som_codes, rs = c(1, 2, 3), colsUsed)
    expect_s3_class(p_pie, "ggplot")

    somRasterData <- data.frame(
        x = rep(seq_len(5), length.out = nrow(som_codes)),
        y = rep(seq_len(2), each = nrow(som_codes) / 2),
        id = seq_len(nrow(som_codes))
    )
    for (m in colsUsed) {
        somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
    }

    p_raster <- .buildBaseRasterGgplot(somRasterData, colsUsed)
    expect_s3_class(p_raster, "ggplot")

    p_select <- .buildSOMRasterSelectPlot(somRasterData, rs = c(1, 2))
    expect_s3_class(p_select, "ggplot")
})

test_that(".buildFlowSOMPiePlot() returns NULL for invalid rs", {
    sce <- CySA_example_sce()
    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    colsUsed <- S4Vectors::metadata(sce)$map$colsUsed

    expect_null(.buildFlowSOMPiePlot(som_codes, rs = integer(0), colsUsed))
    expect_null(.buildFlowSOMPiePlot(som_codes, rs = 9999, colsUsed))
})
