#!/usr/bin/env Rscript

# CySA standalone Shiny app launcher.
#
# This file is a runnable example for local development. Bioconductor packages
# should keep the Shiny app factory inside R/ and expose it via
# clusterSelector(); app.R is optional and only needed when you want to launch
# the app directly from the package root (for example with shiny::runApp()).

pkg_path <- if (file.exists("DESCRIPTION")) {
    getwd()
} else {
    system.file(package = "CySA")
}

if (pkg_path == "" || !file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    stop("Cannot locate CySA package root. Run this script from the package root.")
}

# During development, load the local package. Once installed, library(CySA) is
# enough.
if (file.exists(file.path(pkg_path, "R"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
} else {
    library(CySA)
}

# Build a minimal example dataset and the inputs required by clusterSelector().
sce <- CySA_example_sce(n_cells = 500, n_nodes = 25)
prepped <- prepClusterSelectorData(sce, total_cells_to_sample = 300)

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
    y = rep(seq_len(5), each = nrow(som_codes) / 5),
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

# Create the Shiny application object.
app <- clusterSelector(
    sce = prepped$sce,
    sce_subsampled = prepped$sce_subsampled,
    dList = prepped$dList,
    dend = dend,
    dendTable = dendTable,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = somRasterObj
)

# Launch the app when this script is run directly (e.g. Rscript app.R).
# When sourced by shiny::runApp() or during tests, just return the app object.
if (sys.nframe() == 0L) {
    shiny::runApp(app)
}

app
