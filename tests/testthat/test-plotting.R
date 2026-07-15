test_that("plotSOMScatter() returns a ggplot object", {
    sce <- CySA_example_sce()
    p <- plotSOMScatter(sce, chs = c("marker1", "marker2"))

    expect_s3_class(p, "ggplot")
})

test_that("plotSOMScatter() handles missing channels gracefully", {
    sce <- CySA_example_sce()
    p <- plotSOMScatter(sce, chs = c("missing_marker", "marker1"))

    expect_true(is.null(p) || inherits(p, "ggplot"))
})

test_that("plotScatterBJ() returns a ggplot object", {
    sce <- CySA_example_sce()
    p <- plotScatterBJ(sce, chs = c("marker1", "marker2"))

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
