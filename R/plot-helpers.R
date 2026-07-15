# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# plot-helpers.R ----
# Plot building helpers extracted from the clusterSelector() server closure.
# These functions produce static ggplot2 objects that are used both in the app
# and can be tested independently.

#' Build a Base SOM Raster Heatmap
#'
#' Creates a faceted ggplot2 raster heatmap from \code{somRasterData}. The
#' resulting object is reused as the background for selected-node overlays.
#'
#' @param somRasterData Data frame with \code{x}, \code{y}, and marker columns.
#' @param colsUsed Marker columns to include.
#'
#' @return A \code{ggplot} object or \code{NULL} if \code{somRasterData} is
#'   \code{NULL}.
#'
#' @keywords internal
.buildBaseRasterGgplot <- function(somRasterData, colsUsed) {
    if (is.null(somRasterData)) return(NULL)

    somRasterData <- somRasterData[, c("x", "y", colsUsed), drop = FALSE]
    raster_long <- tidyr::pivot_longer(
        somRasterData,
        cols = c(-"x", -"y"),
        names_to = "marker",
        values_to = "value"
    )
    marker_cols <- setdiff(names(somRasterData), c("x", "y"))
    raster_long$marker <- factor(raster_long$marker, levels = marker_cols)

    ggplot(raster_long, aes(x = x, y = y, fill = value)) +
        geom_raster() +
        facet_wrap(~marker) +
        scale_fill_viridis() +
        coord_fixed() +
        theme_minimal()
}


#' Build a FlowSOM Marker Pie Plot
#'
#' Creates a faceted pie chart showing mean marker expression per selected SOM
#' node.
#'
#' @param somCodes SOM codes matrix.
#' @param rs Selected SOM node ids.
#' @param colsUsed Marker names to include.
#'
#' @return A \code{ggplot} object or \code{NULL} when no markers are available.
#'
#' @keywords internal
.buildFlowSOMPiePlot <- function(somCodes, rs, colsUsed) {
    rs <- as.integer(rs)
    if (length(rs) < 1 || any(rs <= 0) || any(rs > nrow(somCodes))) {
        return(NULL)
    }

    markers <- intersect(colsUsed, colnames(somCodes))
    if (length(markers) < 1) return(NULL)

    expr <- as.data.frame(somCodes[rs, markers, drop = FALSE])
    expr$cluster_id <- factor(as.character(rs), levels = unique(as.character(rs)))

    long <- tidyr::pivot_longer(
        expr,
        cols = tidyr::all_of(markers),
        names_to = "marker",
        values_to = "expression"
    )
    long$marker <- factor(long$marker, levels = markers)

    n_markers <- length(markers)
    marker_cols <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Set2")
    )(n_markers)

    ggplot(long, aes(x = "", y = expression, fill = marker)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y") +
        facet_wrap(~cluster_id) +
        scale_fill_manual(values = marker_cols) +
        theme_minimal() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        )
}


#' Build a SOM Raster Selection Plot
#'
#' Creates a plotly-friendly ggplot showing the full SOM grid with selected
#' nodes highlighted.
#'
#' @param somRasterData Data frame containing \code{x} and \code{y} grid columns.
#' @param rs Selected SOM node ids.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
.buildSOMRasterSelectPlot <- function(somRasterData, rs) {
    n_rows <- max(somRasterData$x)
    n_cols <- max(somRasterData$y)
    data.points <- expand.grid(seq_len(n_rows), seq_len(n_cols))
    colnames(data.points) <- c("x", "y")

    ggplot(data.points, aes(x = x, y = y, customdata = seq_len(nrow(data.points)))) +
        geom_point() +
        geom_point(
            data = data.points[rs, , drop = FALSE],
            aes(x = x, y = y),
            color = "red",
            size = 0.9
        )
}


#' Compute Percentile Axis Limits for Scatter Plots
#'
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param pctl Quantile between 0.5 and 1.
#'
#' @return A named list with \code{xlim} and \code{ylim}.
#'
#' @keywords internal
.computeScatterLimits <- function(x, y, pctl) {
    if (is.null(pctl)) pctl <- 0.99
    tailP <- (1 - pctl) / 2
    list(
        xlim = stats::quantile(x, probs = c(tailP, 1 - tailP), na.rm = TRUE),
        ylim = stats::quantile(y, probs = c(tailP, 1 - tailP), na.rm = TRUE)
    )
}
