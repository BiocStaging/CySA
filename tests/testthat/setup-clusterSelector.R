# Setup helpers for clusterSelector tests.
#
# This file creates a reusable small example app and shared inputs. It is sourced
# once by testthat before any test file runs, so objects created here are
# available in all tests (but not exported from the package).

library(shiny)

# Build a small, deterministic app object that tests can share.
make_test_app <- function(n_cells = 200, n_nodes = 10) {
    sce <- CySA_example_sce(n_cells = n_cells, n_nodes = n_nodes)
    prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 100)

    som_codes <- S4Vectors::metadata(sce)$SOM_codes
    markers <- S4Vectors::metadata(sce)$map$colsUsed

    dend <- stats::as.dendrogram(stats::hclust(stats::dist(som_codes)))
    dendTable <- data.frame(
        id = seq_len(nrow(som_codes)),
        label = rownames(som_codes),
        stringsAsFactors = FALSE
    )

    clusterPatientTable <- table(
        sample_id = sce$sample_id,
        cluster_id = sce$cluster_id
    )

    somRasterData <- data.frame(
        x = rep(seq_len(5), length.out = nrow(som_codes)),
        y = rep(seq_len(2), each = nrow(som_codes) / 2),
        id = seq_len(nrow(som_codes))
    )
    for (m in markers) {
        somRasterData[[m]] <- seq_len(nrow(som_codes)) / nrow(som_codes)
    }

    arr <- array(
        data = seq_len(10 * 10 * length(markers)),
        dim = c(10, 10, length(markers))
    )
    somRasterObj <- raster::brick(arr)
    names(somRasterObj) <- markers

    env <- new.env()
    app <- clusterSelector(
        sce = prepped$sce,
        sce_subsampled = prepped$sce_subsampled,
        dList = prepped$dList,
        dend = dend,
        dendTable = dendTable,
        clusterPatientTable = clusterPatientTable,
        somRasterData = somRasterData,
        somRasterObj = somRasterObj,
        env = env
    )

    list(app = app, env = env, sce = prepped$sce)
}
