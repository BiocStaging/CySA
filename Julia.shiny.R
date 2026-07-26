setwd("/pasteur/helix/projects/scBiomarkers/bernd/cytometry/CySA")
roxygen2::roxygenise()
devtools::document()
devtools::test()
devtools::check()
devtools::check_built()
cov <- covr::codecov()
BiocCheck::BiocCheck()

# After using any workaround above, confirm coverage
result <- BiocCheck::BiocCheck(
    "/pasteur/helix/projects/scBiomarkers/bernd/cytometry/CyFj11",
    `quit-with-status` = FALSE   # collect all results, don't stop on ERROR
)

# Inspect what was checked
names(result)           # "error", "warning", "note"
result$error            # must be length 0 before Bioconductor submission
result$warning
result$note


system("R CMD INSTALL .")
detach("package:CySA", unload = TRUE)
library(CySA)
devtools::load_all()
cp = load(file.path(
    "/pasteur/helix/projects/scBiomarkers/projects/julia",
    "channel_rename",
    "keepChannels",
    "train",
    "indexed",
    "rmMargins",
    "noDoublets",
    "cyEt",
    "transformed",
    "PeacoQC_results",
    "cyCombine",
    "scaleP",
    "sce.som.dim.40.rlen.200.mst.1.radius.40.seed.47",
    "_data",
    "sce.shiny.RData"
))


app <- clusterSelector(
    sce = sce,
    sce_subsampled = sce_subsampled,
    dList = dList,
    dend = dend_obj,
    dendTable = dend_table,
    clusterPatientTable = clusterPatientTable,
    somRasterData = somRasterData,
    somRasterObj = somRasterObj
)



shiny::runApp(app)
