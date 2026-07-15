test_that("prepClusterSelectorData() returns required components", {
    sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
    prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)

    expect_type(prepped, "list")
    expect_named(prepped, c("sce", "sce_subsampled", "dList"))
    expect_s4_class(prepped$sce, "SingleCellExperiment")
    expect_s4_class(prepped$sce_subsampled, "SingleCellExperiment")
    expect_type(prepped$dList, "list")
})

test_that("prepClusterSelectorData() subsampling reduces cell count", {
    sce <- CySA_example_sce(n_cells = 500, n_nodes = 10)
    prepped <- prepClusterSelectorData(
        sce,
        total_cells_to_sample = 100,
        seed = 42
    )
    expect_lt(ncol(prepped$sce_subsampled), ncol(sce))
    expect_true(ncol(prepped$sce_subsampled) >= 100)
})

test_that("prepClusterSelectorData() works with user-provided dList", {
    sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
    dList <- list(
        d1 = c("marker1", "marker2"),
        d2 = c("marker3", "marker4")
    )
    prepped <- prepClusterSelectorData(
        sce,
        dList = dList,
        total_cells_to_sample = 100
    )
    expect_equal(prepped$dList, dList)
})

test_that("prepClusterSelectorData() is reproducible with seed", {
    sce <- CySA_example_sce(n_cells = 200, n_nodes = 10)
    prepped1 <- prepClusterSelectorData(
        sce,
        total_cells_to_sample = 100,
        seed = 123
    )
    prepped2 <- prepClusterSelectorData(
        sce,
        total_cells_to_sample = 100,
        seed = 123
    )
    expect_equal(
        SingleCellExperiment::colData(prepped1$sce_subsampled)$id,
        SingleCellExperiment::colData(prepped2$sce_subsampled)$id
    )
})
