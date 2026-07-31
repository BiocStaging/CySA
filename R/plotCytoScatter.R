#' Plot Scatter for a SingleCellExperiment
#'
#' Generates a 2-D scatter / density / colored point plot for a
#' \code{SingleCellExperiment}.
#'
#' @param x A \code{SingleCellExperiment}.
#' @param chs Character vector of one or two markers/channels to plot. More than
#'   two channels will be melted and faceted.
#' @param gate Currently unused.
#' @param color_by Optional column name used to color points.
#' @param facet_by Optional column name used for faceting.
#' @param bins Number of bins for the 2-D histogram (default 100).
#' @param assay Assay name to use (default \code{"exprs"}).
#' @param label Axis label style: \code{"target"}, \code{"channel"} or
#'   \code{"both"}.
#' @param zeros Logical: if \code{TRUE}, remove rows where both plotted
#'   dimensions are zero.
#' @param k_pal Color palette used when \code{color_by} is a clustering.
#' @param rownms Optional character vector used instead of \code{rownames(x)} for
#'   matching \code{chs}.
#'
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_gradientn
#' @importFrom ggplot2 scale_color_gradientn scale_color_manual facet_wrap
#' @importFrom ggplot2 facet_grid vars ylab guides guide_legend theme_bw theme
#' @importFrom ggplot2 element_blank element_text element_rect
#' @importFrom rlang .data sym
#' @importFrom grid unit
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment int_colData
#' @importFrom CATALYST channels cluster_codes cluster_ids
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' sce <- CySA_example_sce()
#' plotCytoScatter(sce, chs = c("marker1", "marker2"))
#'
#' @export
plotCytoScatter <- function(x, chs, gate = NULL, color_by = NULL, facet_by = NULL,
                            bins = 100, assay = "exprs",
                            label = c("target", "channel", "both"),
                            zeros = FALSE, k_pal = NULL, rownms = NULL) {
  label <- match.arg(label)

  if (!is.null(gate)) {
    warning("'gate' is not currently implemented and will be ignored.")
  }

  # match channels / markers
  m <- if (is.null(rownms)) rownames(x) else rownms
  chx <- CATALYST::channels(x)

  i <- match(chs, m, nomatch = 0L)
  if (all(i == 0L)) {
    i <- match(chs, chx, nomatch = 0L)
  }

  miss <- chs[i == 0L]
  if (length(miss)) {
    stop("Unknown channels/markers: ", paste(miss, collapse = ", "))
  }

  # extract assay
  y <- SummarizedExperiment::assay(x[i, , drop = FALSE], assay)

  # label axes
  nms <- switch(label,
    target = m,
    channel = chx,
    both = ifelse(chx == m, chx, paste(chx, m, sep = "-"))
  )
  chs[i != 0L] <- rownames(y) <- nms[i]

  # add clustering column if requested
  if (!is.null(color_by) &&
    color_by %in% names(CATALYST::cluster_codes(x))) {
    x[[color_by]] <- CATALYST::cluster_ids(x, color_by)
  }

  # build data.frame
  cd <- cbind(
    SummarizedExperiment::colData(x),
    SingleCellExperiment::int_colData(x)
  )
  df <- data.frame(
    t(as.matrix(y)),
    cd,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  cd_vars <- intersect(names(cd), names(df))

  facet <- NULL
  ylab <- NULL

  if (length(chs) > 2L) {
    df <- reshape2::melt(
      df,
      id.vars = unique(c(chs[1L], cd_vars))
    )
    facet <- "variable"
    ylab <- ggplot2::ylab(NULL)
    chs[2L] <- "value"
  }

  # geometry + scales
  if (is.null(color_by)) {
    dc <- if (nrow(df) >= 5L) {
      grDevices::densCols(
        x = df[[chs[1L]]],
        y = df[[chs[2L]]],
        nbin = min(as.integer(bins), 128L),
        colramp = grDevices::colorRampPalette(
          c("navy", rev(RColorBrewer::brewer.pal(11L, "Spectral")))
        )
      )
    } else {
      rep("#000080", nrow(df))
    }
    df$.dens_col <- dc
    geom <- ggplot2::geom_point(
      ggplot2::aes(color = .dens_col),
      size = 0.5, na.rm = TRUE, show.legend = FALSE
    )
    scales <- ggplot2::scale_color_identity(guide = "none")
    guides <- NULL
  } else {
    geom <- ggplot2::geom_point(
      alpha = 0.2, size = 0.8, na.rm = TRUE
    )
    if (is.numeric(df[[color_by]])) {
      scales <- ggplot2::scale_color_gradientn(
        colors = c(
          "navy",
          rev(RColorBrewer::brewer.pal(11L, "Spectral"))
        )
      )
      guides <- NULL
    } else {
      if (color_by %in% names(CATALYST::cluster_codes(x))) {
        if (is.null(k_pal)) k_pal <- CySA_default_cluster_cols()
        scales <- ggplot2::scale_color_manual(values = k_pal)
      } else {
        scales <- NULL
      }
      guides <- ggplot2::guides(
        col = ggplot2::guide_legend(
          override.aes = list(alpha = 1, size = 3)
        )
      )
    }
  }

  # facets
  facet <- c(facet, facet_by)
  if (!is.null(facet)) {
    if (length(facet) == 1L) {
      facet <- ggplot2::facet_wrap(
        ggplot2::vars(!!rlang::sym(facet))
      )
    } else {
      facet <- ggplot2::facet_grid(
        cols = ggplot2::vars(!!rlang::sym(facet[1L])),
        rows = ggplot2::vars(!!rlang::sym(facet[2L]))
      )
    }
  }

  # optional zero filtering
  if (isTRUE(zeros)) {
    df <- df[rowSums(
      df[, chs[c(1L, 2L)], drop = FALSE] == 0
    ) == 0, , drop = FALSE]
  }

  # plot
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[chs[1L]]],
      y = .data[[chs[2L]]],
      col = if (is.null(color_by)) NULL else .data[[color_by]]
    )
  ) +
    geom +
    scales +
    guides +
    facet +
    ylab +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "white"),
      legend.key.height = grid::unit(0.8, "lines")
    )
}
