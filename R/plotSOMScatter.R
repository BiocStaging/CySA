#' @rdname plotSOMScatter
#' @title Scatter plot
#'
#' @description Bivariate scatter plots including visualization of
#' (group-specific) gates, their boundaries and percentage of selected cells.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character vector specifying which channels to plot.
#' @param metaSlot name of the metadata slot containing SOM codes.
#' @param pointSize column in SOM stats used to size points.
#' @param color_by column to color points by.
#' @param bins number of bins when coloring by density.
#' @param assay name of the assay to use.
#' @param statsSlot name of the metadata slot containing SOM stats.
#' @param label how axis labels should be constructed.
#' @param zeros logical specifying whether to include 0 values.
#' @param k_pal optional cluster color palette.
#' @param xRN optional row names.
#' @param xCN optional channel names.
#'
#' @return a \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#' @importFrom ggplot2 scale_color_manual facet_wrap facet_grid vars
#' @importFrom ggplot2 guides guide_legend ylab theme_bw theme
#' @importFrom ggplot2 element_blank element_text element_rect
#' @importFrom grid unit
#' @importFrom rlang sym .data
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom SingleCellExperiment int_colData
#' @importFrom CATALYST channels cluster_codes cluster_ids
#' @importFrom S4Vectors metadata
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' sce <- CySA_example_sce()
#' plotSOMScatter(sce, chs = c("marker1", "marker2"))
#'
#' @export
plotSOMScatter <- function(x, chs, metaSlot = "SOM_codes", pointSize = "n",
                           color_by = "n",
                           bins = 100, assay = "exprs", statsSlot = "SOM_stats",
                           label = c("target", "channel", "both"),
                           zeros = FALSE, k_pal = NULL,
                           xRN = NULL, xCN = NULL) {
    label <- match.arg(label)
    pointSize <- pointSize[1]
    color_by <- color_by[1]

    # ---- validate metadata slots ----
    md <- S4Vectors::metadata(x)
    if (!metaSlot %in% names(md)) {
        stop("Need '", metaSlot, "' in metadata of x", call. = FALSE)
    }
    som_codes <- md[[metaSlot]]

    if (!is.null(statsSlot) && !statsSlot %in% names(md)) {
        warning("'", statsSlot, "' not found in metadata of x; proceeding without stats")
        statsSlot <- NULL
    }
    if (!is.null(statsSlot)) {
        stats <- md[[statsSlot]]
    }

    # ---- resolve names ----
    if (is.null(xRN)) xRN <- rownames(x)
    if (is.null(xCN)) xCN <- CATALYST::channels(x)

    # ---- validate channels ----
    valid <- unique(c(colnames(som_codes), xRN, xCN))
    missing <- setdiff(chs, valid)
    if (length(missing) == length(chs)) {
        warning(
            "Unknown channels/markers: ", paste(missing, collapse = ", "),
            "\nAvailable in SOM_codes: ", paste(colnames(som_codes), collapse = ", "),
            call. = FALSE
        )
        return(NULL)
    }
    if (length(missing)) {
        warning(
            "Unknown channels/markers dropped: ", paste(missing, collapse = ", "),
            call. = FALSE
        )
        chs <- setdiff(chs, missing)
    }

    # ---- ensure requested channels exist in SOM_codes (compute locally, do not mutate x) ----
    y <- SummarizedExperiment::assay(x, assay)
    assay_rows <- rownames(y)

    missing_in_som <- setdiff(chs, colnames(som_codes))
    if (length(missing_in_som)) {
        for (ch in missing_in_som) {
            if (!ch %in% assay_rows) {
                stop(
                    "Cannot compute SOM stats for '", ch,
                    "': not found in assay rows.",
                    call. = FALSE
                )
            }
            warning("computing SOM stats for '", ch, "'")
            cl_ids <- as.integer(unique(SummarizedExperiment::colData(x)$cluster_id))
            vals <- vapply(cl_ids, function(cl) {
                cells <- which(SummarizedExperiment::colData(x)$cluster_id == cl)
                mean(y[ch, cells, drop = TRUE])
            }, numeric(1))
            som_codes <- cbind(som_codes, stats::setNames(vals, NULL))
            colnames(som_codes)[ncol(som_codes)] <- ch
        }
    }

    yy <- som_codes[, chs, drop = FALSE]

    # ---- construct label names ----
    nms <- switch(label,
                  target = xRN,
                  channel = xCN,
                  both = ifelse(xCN == xRN, xCN, paste(xCN, xRN, sep = "-"))
    )
    # rename SOM-code columns for plotting
    # (only rename those that correspond to rownames(x)/channels)
    rename_map <- stats::setNames(nms, xRN)
    found <- chs %in% names(rename_map)
    colnames(yy)[found] <- rename_map[chs[found]]

    # ---- add cluster color column if requested ----
    if (isTRUE(color_by %in% names(CATALYST::cluster_codes(x)))) {
        x[[color_by]] <- CATALYST::cluster_ids(x, color_by)
    }

    # ---- build data.frame ----
    if (!is.null(statsSlot)) {
        df <- cbind(yy, stats)
    } else {
        df <- as.data.frame(yy)
        df$n <- 1
    }

    if (is.null(chs) || length(chs) < 2) {
        stop("At least two channels are required for a scatter plot.", call. = FALSE)
    }

    # filter zero rows
    if (!zeros) {
        df <- df[rowSums(df[, chs[c(1, 2)], drop = FALSE] == 0) == 0, , drop = FALSE]
    }

    # ---- melt if >2 channels ----
    if (length(chs) > 2) {
        id_vars <- if (!is.null(statsSlot)) unique(c(chs[1], colnames(stats))) else chs[1]
        df <- reshape2::melt(df, id.vars = id_vars)
        facet <- "variable"
        ylab <- ggplot2::ylab(NULL)
        chs[2] <- "value"
    } else {
        facet <- NULL
        ylab <- NULL
    }

    # ---- color/size setup ----
    col_var <- if (color_by %in% colnames(df)) color_by else NULL

    if (!pointSize %in% colnames(df)) pointSize <- NULL
    if (!is.null(pointSize) && pointSize == "n" && !"n" %in% colnames(df)) {
        pointSize <- NULL
    }

    geom <- ggplot2::geom_point(alpha = 0.2, na.rm = TRUE)
    guides <- NULL
    scales <- NULL

    if (!is.null(col_var) && is.numeric(df[[col_var]])) {
        scales <- ggplot2::scale_color_gradientn(
            colors = c("navy", rev(RColorBrewer::brewer.pal(11, "Spectral"))),
            guide = "colorbar"
        )
    } else if (!is.null(col_var) && col_var %in% names(CATALYST::cluster_codes(x))) {
        if (is.null(k_pal)) k_pal <- CySA_default_cluster_cols()
        scales <- ggplot2::scale_color_manual(values = k_pal)
        guides <- ggplot2::guides(col = ggplot2::guide_legend(
            override.aes = list(alpha = 1, size = 3)
        ))
    }

    # ---- facets ----
    if (!is.null(facet)) {
        if (length(facet) == 1) {
            facet <- ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet)))
        } else {
            facet <- ggplot2::facet_grid(
                cols = ggplot2::vars(!!rlang::sym(facet[1])),
                rows = ggplot2::vars(!!rlang::sym(facet[2]))
            )
        }
    }

    # ---- filter zeros again after melting ----
    if (!zeros) {
        df <- df[rowSums(df[, chs[c(1, 2)], drop = FALSE] == 0) == 0, , drop = FALSE]
    }

    # ---- aesthetics ----
    aes_args <- list(
        x = rlang::sym(chs[1]),
        y = rlang::sym(chs[2])
    )
    if (!is.null(pointSize)) aes_args$size <- rlang::sym(pointSize)
    if (!is.null(col_var)) aes_args$color <- rlang::sym(col_var)

    p1 <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) + geom

    if (!is.null(scales)) p1 <- p1 + scales
    if (!is.null(guides)) p1 <- p1 + guides
    if (!is.null(ylab)) p1 <- p1 + ylab
    if (!is.null(facet)) p1 <- p1 + facet

    p1 + ggplot2::theme_bw() + ggplot2::theme(
        aspect.ratio = 1,
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(color = "black"),
        strip.background = ggplot2::element_rect(fill = "white"),
        legend.key.height = grid::unit(0.8, "lines")
    )
}
