library(CySA)

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
