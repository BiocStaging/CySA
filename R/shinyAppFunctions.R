# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# functions used in shiny app
#' @importFrom tidyselect all_of
# safe_event_data ----
# Wrapper around plotly::event_data that suppresses the "not registered"
# warnings unless verbose mode is enabled.
safe_event_data <- function(event, source = "all", verbose = FALSE) {
    if (verbose) {
        return(plotly::event_data(event, source = source))
    }
    withCallingHandlers(
        plotly::event_data(event, source = source),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("is not registered", msg) || grepl("event_register", msg)) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

# Alias with dot prefix used by extracted server helpers.
.safeEventData <- safe_event_data

# quiet_ggplotly ----
# Wrapper around ggplotly that suppresses the harmless ggplot2 warning about
# the customdata aesthetic, which is consumed by plotly after conversion.
quiet_ggplotly <- function(p, ...) {
    result <- withCallingHandlers(
        suppressMessages(ggplotly(p, ...)),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("unknown aesthetics.*(customdata|key)", msg, ignore.case = TRUE) ||
                grepl("Aspect ratios aren't yet implemented", msg, ignore.case = TRUE) ||
                grepl("has yet to be implemented in plotly",  msg, ignore.case = TRUE) ||
                grepl("don't have these attributes",          msg, ignore.case = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
    # Remove ggplot2-internal 'colour' attribute from every trace.
    result$x$data <- lapply(result$x$data, function(tr) {
        tr[["colour"]] <- NULL
        tr
    })
    result <- result %>%
        plotly::layout(autosize = TRUE) %>%
        plotly::config(responsive = TRUE)
    # Force early evaluation of the lazy plotly object. plotly_build() is what
    # actually resolves trace type/mode and emits the "No trace type
    # specified"/"No scatter mode specified" messages -- doing it here, inside
    # suppressMessages(), means Shiny's later renderPlotly()/toJSON step (or a
    # test accessing output$xxx) just serializes an already-resolved object
    # and has nothing left to emit.
    suppressMessages(plotly::plotly_build(result))
}

# highlight_df = function ----
highlight_df <- function(x, y, rs, somCodesName = "SOM_codes", metaD = NULL) {
    if (is.null(metaD)) stop("'metaD' must be provided")
    if (is.null(rs) || length(rs) == 0L) return(data.frame(x = numeric(0), y = numeric(0), id = integer(0)))

    # Filter out NA and out-of-bounds indices
    somCodes <- metaD[[somCodesName]]
    rs_valid <- rs[!is.na(rs) & rs > 0L & rs <= nrow(somCodes)]

    if (length(rs_valid) == 0L) return(data.frame(x = numeric(0), y = numeric(0), id = integer(0)))

    data.frame(
        x = somCodes[rs_valid, x],
        y = somCodes[rs_valid, y],
        id = rs_valid
    )
}


countBarPlotFunc <- function(rs, clusterPatientTable, cst, sce, outputList, groupsInput) {
    rs <- as.integer(intersect(rs, colnames(clusterPatientTable)))
    # selected counts
    rSums <- data.frame(counts = rowSums(clusterPatientTable[, rs, drop = FALSE]))
    rSums$id <- rownames(clusterPatientTable)

    # counts to compare to
    if (!cst == "none") {
        if (cst %in% colnames(S4Vectors::metadata(sce)$experiment_info)) {
            # expInfo = S4Vectors::metadata(sce)$experiment_info
            numCols <- unlist(lapply(S4Vectors::metadata(sce)$experiment_info, is.numeric), use.names = FALSE)
            expInfo <- S4Vectors::metadata(sce)$experiment_info[, numCols, drop = FALSE]
            # eI = expInfo[,!colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"),drop=FALSE]
            eI <- apply(expInfo, 2, as.numeric) %>% as.data.frame()
            rownames(eI) <- S4Vectors::metadata(sce)$experiment_info[, "sample_id"]
            if (is.na(eI) %>% any()) {
                stop("NAs produced when converting experiment_info to numeric")
            }
            ctStats <- eI[, cst]
            rSums <- cbind(rSums, eI[rSums$id, cst])
            colnames(rSums) <- c("counts", "id", cst)
            dfm <- reshape2::melt(rSums[, c("id", cst, "counts")], id.vars = 1)
            colnames(dfm) <- c("id", "variable", "counts")
        } else if (cst %in% names(outputList)) {
            rSums <- cbind(rSums, rowSums(clusterPatientTable[, outputList[[cst]], drop = FALSE]))
            colnames(rSums) <- c("counts", "id", cst)
            dfm <- reshape2::melt(rSums[, c("id", cst, "counts")], id.vars = 1)
            colnames(dfm) <- c("id", "variable", "counts")
        } else {
            stop("Comparison column '", cst, "' not found in experiment_info or outputList")
        }
    } else {
        dfm <- rSums
        dfm$variable <- "counts"
    }
    if ("group" %in% colnames(S4Vectors::metadata(sce)$experiment_info)) {
        dfm <- merge(dfm, S4Vectors::metadata(sce)$experiment_info, by.x = "id", by.y = "sample_id")
        lvs <- S4Vectors::metadata(sce)$experiment_info[order(S4Vectors::metadata(sce)$experiment_info$group), "sample_id"]
        # dfm = dfm %>% group_by(group)
        dfm$id <- factor(dfm$id, levels = lvs)
    }
    dfm$groups <- "other"
    if (!is_empty(groupsInput)) {
        if (!is_empty(groupsInput$group1)) {
            dfm$groups[dfm$id %in% groupsInput$group1] <- "group1"
        }
        if (!is_empty(groupsInput$group2)) {
            dfm$groups[dfm$id %in% groupsInput$group2] <- "group2"
        }
    }
    dfm$groups <- as.factor(dfm$groups)
    dfm <- as_tibble(dfm)
    # save(file = "/pasteur/appa/scratch/bernd/Rtest2.Rdata",
    #      list = c("clusterPatientTable", "rs", "expInfo", "eI", "cst", "rSums", "dfm", "groupsInput"))
    # cp =load(file = "/pasteur/appa/scratch/bernd/Rtest2.Rdata")
    ggplot(dfm, aes(x = reorder(id, counts), y = counts)) +
        geom_bar(aes(fill = interaction(variable, groups)), stat = "identity", position = "dodge") +
        coord_flip() +
        guides(fill = guide_legend(title = "groups"))
    # ggplot(dfm, aes(x=id, y=counts)) +
    #   geom_bar(stat = "identity") + coord_flip()
}


# Deprecated alias kept for backward compatibility.
# Prefer .computeRelativeCounts() in new code.
compute_relative_counts <- function(clusterPatientTable, rs, relativeToCol, expInfo, outputList) {
    switch(relativeToCol,
           "none" = rowSums(clusterPatientTable[, rs, drop = FALSE]) /
               rowSums(clusterPatientTable[, , drop = FALSE]) * 100,
           {
               if (relativeToCol %in% colnames(expInfo)) {
                   rowSums(clusterPatientTable[, rs, drop = FALSE]) /
                       expInfo[rownames(clusterPatientTable), relativeToCol] * 100
               } else if (relativeToCol %in% names(outputList)) {
                   rowSums(clusterPatientTable[, rs, drop = FALSE]) /
                       rowSums(clusterPatientTable[, outputList[[relativeToCol]], drop = FALSE]) * 100
               } else {
                   stop("relativeToCol '", relativeToCol, "' not found")
               }
           }
    )
}

PercentBarPlotFunc <- function(sce, relativeToCol, clusterPatientTable, rs, outputList, group, groupsInput) {
    expInfo <- metadata(sce)$experiment_info
    rownames(expInfo) <- expInfo$sample_id
    expInfo <- expInfo[, !colnames(expInfo) %in% c("sample_nr", "sample_id", "sample"), drop = FALSE]

    rSums <- data.frame(
        Percent = compute_relative_counts(
            clusterPatientTable, rs, relativeToCol, expInfo, outputList
        )
    )
    rSums$id <- rownames(clusterPatientTable)

    if ("group" %in% colnames(S4Vectors::metadata(sce)$experiment_info)) {
        rSums <- merge(rSums, S4Vectors::metadata(sce)$experiment_info, by.x = "id", by.y = "sample_id")
        rSums <- rSums %>% group_by(group)
        rSums$id <- factor(rSums$id, levels = rSums[order(rSums$group), "id"][[1]])
    }
    rSums$groups <- "other"
    if (!is_empty(groupsInput)) {
        if (!is_empty(groupsInput$group1)) {
            rSums$groups[rSums$id %in% groupsInput$group1] <- "group1"
        }
        if (!is_empty(groupsInput$group2)) {
            rSums$groups[rSums$id %in% groupsInput$group2] <- "group2"
        }
    }
    rSums$groups <- as.factor(rSums$groups)
    rSums <- as_tibble(rSums)
    # save(file = "/pasteur/appa/scratch/bernd/Rtest.Rdata",
    #      list = c("clusterPatientTable", "rs", "expInfo", "eI", "relativeToCol", "rSums", "groupsInput"))
    # cp = load(file = "/pasteur/appa/scratch/bernd/Rtest.Rdata")
    ggplot(rSums, aes(x = reorder(id, Percent), y = Percent, fill = groups)) +
        # geom_bar(stat = "identity") + coord_flip()
        geom_col() +
        coord_flip()
}


ggsomPlot <- function(pp1, plotIdx, rs, dimSelection, somCodesName = "SOM_codes", sce, metaD = S4Vectors::metadata(sce), xlim = NULL, ylim = NULL) {
    newData <- highlight_df(dimSelection[[plotIdx]]$dims[1], dimSelection[[plotIdx]]$dims[2], rs, somCodesName, metaD = metaD)
    p3 <- pp1 + geom_point(
        data = newData,
        aes(x = `x`, y = `y`),
        color = "red",
        size = 0.3
    )
    if (!is.null(xlim) && !is.null(ylim)) {
        p3 <- p3 + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)  # clips, not drops

    } else if (is.null(dimSelection[[plotIdx]]$xzoom[1])) {
        p3 <- p3 + ggplot2::coord_cartesian(
            xlim = c(
                as.numeric(dimSelection[[plotIdx]]$xlim[1]),
                as.numeric(dimSelection[[plotIdx]]$xlim[2])
            ),
            ylim = c(
                as.numeric(dimSelection[[plotIdx]]$ylim[1]),
                as.numeric(dimSelection[[plotIdx]]$ylim[2])
            )
        )
    } else {
        p3 <- p3 + ggplot2::coord_cartesian(
            xlim = c(
                as.numeric(dimSelection[[plotIdx]]$xzoom[1]),
                as.numeric(dimSelection[[plotIdx]]$xzoom[2])
            ),
            ylim = c(
                as.numeric(dimSelection[[plotIdx]]$yzoom[1]),
                as.numeric(dimSelection[[plotIdx]]$yzoom[2])
            )
        )
    }
    return(p3)
}

somPlot <- function(pp1, plotIdx, rs, colorbyGroups, showGroups,
                    dimSelection = NULL, somCodesName = "SOM_codes",
                    sce, metaD = S4Vectors::metadata(sce),
                    outputList = list(), projectionDf = NULL,
                    xlim = NULL, ylim = NULL,
                    source = NULL) {
    if (is.null(pp1)) return(NULL)

    # Always use index 1 for dimSelection - the caller (somDataMain) passes a single-element list
    dsIdx <- 1L

    if (showGroups) {
        req(projectionDf)
        # drawProjection returns a ggplot object
        p3 <- drawProjection(projectionDf, rs,
                             colorbyGroups = colorbyGroups,
                             sce = sce, outputList = outputList)
        # Apply zoom/limits to ggplot before conversion
        if (!is.null(dimSelection[[dsIdx]]$xzoom[1])) {
            p3 <- p3 + ggplot2::coord_cartesian(
                xlim = c(as.numeric(dimSelection[[dsIdx]]$xzoom[1]),
                       as.numeric(dimSelection[[dsIdx]]$xzoom[2])),
                ylim = c(as.numeric(dimSelection[[dsIdx]]$yzoom[1]),
                       as.numeric(dimSelection[[dsIdx]]$yzoom[2]))
            )
        } else if (!is.null(xlim) && !is.null(ylim)) {
            p3 <- p3 + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
        } else if (!is.null(dimSelection[[dsIdx]]$xlim) && !is.null(dimSelection[[dsIdx]]$ylim)) {
            p3 <- p3 + ggplot2::coord_cartesian(
                xlim = c(as.numeric(dimSelection[[dsIdx]]$xlim[1]),
                       as.numeric(dimSelection[[dsIdx]]$xlim[2])),
                ylim = c(as.numeric(dimSelection[[dsIdx]]$ylim[1]),
                       as.numeric(dimSelection[[dsIdx]]$ylim[2]))
            )
        }
    } else {
        # ggplot path - apply limits before converting to plotly
        p3 <- pp1 + geom_point(
            data = highlight_df(
                dimSelection[[dsIdx]]$dims[1],
                dimSelection[[dsIdx]]$dims[2],
                rs, somCodesName, metaD = metaD
            ),
            aes(x = x, y = y), color = "red", size = 0.3
        )

        # Apply zoom limits if set by user, otherwise use provided limits or defaults
        if (!is.null(dimSelection[[dsIdx]]$xzoom[1])) {
            p3 <- p3 + ggplot2::coord_cartesian(
                xlim = c(as.numeric(dimSelection[[dsIdx]]$xzoom[1]),
                       as.numeric(dimSelection[[dsIdx]]$xzoom[2])),
                ylim = c(as.numeric(dimSelection[[dsIdx]]$yzoom[1]),
                       as.numeric(dimSelection[[dsIdx]]$yzoom[2]))
            )
        } else if (!is.null(xlim) && !is.null(ylim)) {
            p3 <- p3 + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
        } else if (!is.null(dimSelection[[dsIdx]]$xlim) && !is.null(dimSelection[[dsIdx]]$ylim)) {
            p3 <- p3 + ggplot2::coord_cartesian(
                xlim = c(as.numeric(dimSelection[[dsIdx]]$xlim[1]),
                       as.numeric(dimSelection[[dsIdx]]$xlim[2])),
                ylim = c(as.numeric(dimSelection[[dsIdx]]$ylim[1]),
                       as.numeric(dimSelection[[dsIdx]]$ylim[2]))
            )
        }
    }

    if (is.null(p3)) return(NULL)

    # Use explicit source when supplied; otherwise derive from plot index
    if (plotIdx > 5L) plotIdx <- 6L
    plot_source <- if (!is.null(source)) source else paste0("somData", plotIdx)

    quiet_ggplotly(p3, source = plot_source, tooltip = "") %>%
        layout(showlegend = FALSE, dragmode = "select") %>%
        event_register("plotly_selected") %>%
        event_register("plotly_relayout")
}

tsneFunc <- function(dimRedSelection, perplexity, sce, somCodesName = "SOM_codes") {
    dimRedCols <- dimRedSelection
    mat <- S4Vectors::metadata(sce)[[somCodesName]][, dimRedCols, drop = FALSE]
    max_perplexity <- max(1, floor((nrow(mat) - 1) / 3))
    perplexity <- min(as.integer(perplexity), max_perplexity)
    Rtsne::Rtsne(mat, perplexity = perplexity)
}


plotViolinFunc <- function(sce, somCodesName = "SOM_codes", upsetSelection, outputList, violinSelection) {
    somCodes <- sce@metadata[[somCodesName]]
    markers <- colnames(somCodes)
    if (length(upsetSelection) < 3) {
        upsetSelection <- names(outputList)
    }
    markers <- intersect(markers, violinSelection)
    nNodes <- nrow(somCodes)
    parts <- lapply(upsetSelection, function(na) {
        rows <- outputList[[na]]
        rows <- rows[rows != 0]
        if (length(rows) < 2) {
            return(NULL)
        }
        wide <- as.data.frame(somCodes[rows, markers, drop = FALSE])
        wide$somNode <- factor(rows, levels = seq_len(nNodes))
        long <- tidyr::pivot_longer(wide,
                                    cols = tidyselect::all_of(markers),
                                    names_to = "marker", values_to = "expr"
        )
        long$grpName <- na
        long
    })
    data <- data.table::rbindlist(parts)
    if (nrow(data) == 0) {
        return(
            ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::labs(title = "No groups with \u2265 2 SOM nodes selected")
        )
    }

    if (nrow(data) > 0) {
        data$marker <- factor(data$marker, levels = markers)
        data$grpName <- factor(data$grpName, levels = upsetSelection)
    }
    if (nrow(data) > 10) {
        nb.cols <- length(unique(data$marker))
        mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

        p <- ggplot(data, aes(factor(marker), expr, fill = marker)) +
            geom_violin(scale = "width", adjust = 1, trim = TRUE) +
            scale_y_continuous(expand = c(0, 0), position = "right", labels = function(x) {
                c(rep(x = "", times = length(x) - 2), x[length(x) - 1], "")
            }) +
            facet_grid(rows = vars(grpName), scales = "free", switch = "y") +
            theme_cowplot(font_size = 12) +
            theme(
                legend.position = "none", panel.spacing = unit(0, "lines"),
                plot.title = element_text(hjust = 0.5),
                panel.background = element_rect(fill = NA, color = "black"),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold"),
                strip.text.y.left = element_text(angle = 0),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            ) +
            scale_fill_manual(values = mycolors) +
            ggtitle("Marker on x-axis") +
            xlab("Marker") +
            ylab("Expression Level")
        return(p)
    } else {
        message("please check that all upsetSelection are in outputList: no violin data to plot")
        return(NULL)
    }
}


#' Build a Projection Data Frame for drawProjection
#'
#' Constructs the per-SOM-node data frame that \code{drawProjection} expects:
#' the two requested channel columns come first, followed by all SOM_stats
#' columns.  Reads coordinates directly from \code{metadata(sce)[[somCodesName]]}
#' instead of parsing \code{ggplot_build()} output, which is both faster and
#' independent of which aesthetics the base plot happens to expose.
#'
#' @param pp1 Base ggplot object (kept in signature for call-site compatibility;
#'   no longer used internally).
#' @param plotIdx Index into \code{dimSelection}.
#' @param dimSelection List of per-plot dimension specs.
#' @param sce Full \code{SingleCellExperiment}.
#' @param somCodesName Name of the SOM codes metadata slot
#'   (default \code{"SOM_codes"}).
#'
#' @return A data frame or \code{NULL} on failure.
#'
#' @keywords internal
buildProjectionDf <- function(pp1, plotIdx, dimSelection, sce,
                              somCodesName = "SOM_codes") {
    shiny::req(pp1, dimSelection)

    ch_names  <- dimSelection[[plotIdx]]$dims
    md        <- S4Vectors::metadata(sce)
    som_codes <- md[[somCodesName]]
    som_stats <- md[["SOM_stats"]]

    ## ---- guard: required metadata slots ------------------------------------
    if (is.null(som_codes)) {
        warning("buildProjectionDf: '", somCodesName,
                "' not found in metadata(sce).")
        return(NULL)
    }
    if (is.null(som_stats)) {
        warning("buildProjectionDf: 'SOM_stats' not found in metadata(sce).")
        return(NULL)
    }

    ## ---- ensure SOM_stats has an id column ---------------------------------
    if (!"id" %in% colnames(som_stats)) {
        som_stats$id <- seq_len(nrow(som_stats))
    }

    ## ---- resolve channels --------------------------------------------------
    available  <- colnames(som_codes)
    missing_ch <- setdiff(ch_names, available)
    if (length(missing_ch) > 0L) {
        warning(
            "buildProjectionDf: channel(s) not in SOM codes - falling back: ",
            paste(missing_ch, collapse = ", ")
        )
        ## substitute with the first available channels
        fallback <- head(available, 2L)
        ch_names[!ch_names %in% available] <-
            fallback[seq_along(missing_ch)]
    }

    ## ---- build coordinate data frame directly from SOM codes ---------------
    ## Row i of som_codes corresponds to SOM node i.
    coord_df <- data.frame(
        id = seq_len(nrow(som_codes)),
        ch1 = som_codes[, ch_names[1L]],
        ch2 = som_codes[, ch_names[2L]],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    colnames(coord_df)[colnames(coord_df) == "ch1"] <- ch_names[1L]
    colnames(coord_df)[colnames(coord_df) == "ch2"] <- ch_names[2L]

    ## ---- join with SOM stats -----------------------------------------------
    df2 <- dplyr::left_join(som_stats, coord_df, by = "id")
    df2 <- df2[order(df2$id), ]
    df2[is.na(df2)] <- 0

    ## ---- channel columns first (drawProjection reads names(df)[1:2]) -------
    coord_idx <- match(ch_names, names(df2))
    other_idx <- setdiff(seq_along(df2), coord_idx)
    df2 <- df2[, c(coord_idx, other_idx), drop = FALSE]
    df2
}
drawProjection <- function(df, rs, colorbyGroups, sce, outputList = list()) {
    colN <- names(df)[seq_len(2)]
    if (!"id" %in% names(df)) {
        df$id <- seq_len(nrow(df))
        df <- dplyr::left_join(df, S4Vectors::metadata(sce)$SOM_stats, by = "id")
    }
    df <- df[order(df$id), ]   # row i == SOM node i  ->  pointNumber+1 fallback works

    if ("n"       %in% names(df) && !"N"      %in% names(df)) df$N      <- df$n
    if ("rdQu"    %in% names(df) && !"thrdQu" %in% names(df)) df$thrdQu <- df$rdQu
    if ("thirdQu" %in% names(df) && !"thrdQu" %in% names(df)) df$thrdQu <- df$thirdQu
    df$cluster <- df$id

    nGrps  <- 1L
    sl     <- FALSE
    mycolors <- NULL
    if (length(colorbyGroups) >= 1L) {
        df$colGrp <- "other"
        for (cg in colorbyGroups) {
            mask <- df$cluster %in% outputList[[cg]]
            df$colGrp[mask] <- ifelse(
                df$colGrp[mask] == "other",
                cg,
                paste(df$colGrp[mask], cg)   # multi-group nodes get combined label
            )
        }
        df$colGrp <- factor(df$colGrp)
        sl    <- TRUE
        nGrps <- length(levels(df$colGrp))
        mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nGrps + 1L)
    }

    rs_int      <- as.integer(rs)
    df_selected <- df[df$id %in% rs_int, , drop = FALSE]

    # --- Build incrementally: fixes the if/else operator-precedence bug ---
    # For t-SNE/UMAP/PCA: use aes(color = colGrp) when grouping is requested.
    # Note: This creates multiple traces in plotly, but selection works because
    # the observer collects points from all traces via the 'key' aesthetic.
    if (nGrps == 1L) {
        p3 <- ggplot2::ggplot(
            data = df,
            ggplot2::aes(
                x    = .data[[colN[1]]],
                y    = .data[[colN[2]]],
                key  = id,              # sufficient for selection
                text = paste(
                    "cluster:", cluster, "<br>N cells:", N, "<br>mean:",
                    format(mean, digits = 2), "<br>3rd Q:",
                    format(thrdQu, digits = 2), "<br>max:",
                    format(max, digits = 2)
                )
            )
        ) +
            ggplot2::geom_point(
                color = "lightblue", size = 1.5,
                show.legend = FALSE, na.rm = TRUE
            )
    } else {
        p3 <- ggplot2::ggplot(
            data = df,
            ggplot2::aes(
                x     = .data[[colN[1]]],
                y     = .data[[colN[2]]],
                color = colGrp,
                key   = id,             # sufficient for selection
                text  = paste(
                    "cluster:", cluster, "<br>N cells:", N, "<br>mean:",
                    format(mean, digits = 2), "<br>3rd Q:",
                    format(thrdQu, digits = 2), "<br>max:",
                    format(max, digits = 2)
                )
            )
        ) +
            ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
            ggplot2::scale_color_manual(values = mycolors, drop = FALSE)
    }

    # Red overlay — no customdata needed, purely visual
    if (nrow(df_selected) > 0L) {
        p3 <- p3 + ggplot2::geom_point(
            data        = df_selected,
            mapping     = ggplot2::aes(
                x = .data[[colN[1]]],
                y = .data[[colN[2]]]   # customdata removed
            ),
            color       = "red",
            size        = 0.8,
            inherit.aes = FALSE,
            na.rm       = TRUE
        )
    }
    p3
}
plotViolin2Func <- function(sce, somCodesName = "SOM_codes", violinSelection, upsetSelection, outputList) {
    somCodes <- sce@metadata[[somCodesName]]
    markers <- colnames(somCodes)
    markers <- intersect(markers, violinSelection)
    if (length(upsetSelection) < 3) {
        upsetSelection <- names(outputList)
    }
    nNodes <- nrow(somCodes)
    parts <- lapply(upsetSelection, function(na) {
        rows <- outputList[[na]]
        rows <- rows[rows != 0]
        if (length(rows) < 2) {
            return(NULL)
        }
        wide <- as.data.frame(somCodes[rows, markers, drop = FALSE])
        wide$somNode <- factor(rows, levels = seq_len(nNodes))
        long <- tidyr::pivot_longer(wide,
                                    cols = tidyselect::all_of(markers),
                                    names_to = "marker", values_to = "expr"
        )
        long$grpName <- na
        long
    })
    data <- data.table::rbindlist(parts)
    if (nrow(data) == 0) {
        return(
            ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::labs(title = "No groups with \u2265 2 SOM nodes selected")
        )
    }

    if (nrow(data) > 0) {
        data$marker <- factor(data$marker, levels = markers)
        data$grpName <- factor(data$grpName, levels = upsetSelection)
    }
    nb.cols <- length(unique(data$grpName))
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

    p <- ggplot(data, aes(factor(grpName), expr, fill = grpName)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_y_continuous(expand = c(0, 0), position = "right", labels = function(x) {
            c(rep(x = "", times = length(x) - 2), x[length(x) - 1], "")
        }) +
        facet_grid(rows = vars(marker), scales = "free", switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(
            legend.position = "none", panel.spacing = unit(0, "lines"),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = NA, color = "black"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.y.left = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) +
        scale_fill_manual(values = mycolors) +
        ggtitle("grp name on x-axis") +
        xlab("Groups") +
        ylab("Expression Level")
    p
}


upsetPlotFunc <- function(upsetSelection, outputList, sce, maxCombs = 100) {
    if (length(upsetSelection) < 2) {
        return(NULL)
    }
    if (is.null(grDevices::dev.list())) {
        grDevices::pdf(NULL)
        on.exit(grDevices::dev.off(), add = TRUE)
    }
    if (length(upsetSelection) > 31) upsetSelection <- head(upsetSelection, 31)
    cm <- ComplexHeatmap::make_comb_mat(outputList[upsetSelection])

    nCombs <- length(ComplexHeatmap::comb_name(cm))
    if (nCombs > maxCombs) {
        stop(
            "Too many distinct group overlap combinations (", nCombs,
            ") to render an UpSet plot; select fewer/less-overlapping groups ",
            "(currently allows up to ", maxCombs, ")."
        )
    }
    ncells <- rep(0, length(ComplexHeatmap::comb_name(cm)))
    grpCells <- rep(0, length(upsetSelection))
    for (cidx in seq_along(ComplexHeatmap::comb_name(cm))) {
        ncells[cidx] <- S4Vectors::metadata(sce)$SOM_stats[as.integer(ComplexHeatmap::extract_comb(cm, ComplexHeatmap::comb_name(cm)[cidx])), "n"] %>% sum()
    }
    for (olIdx in seq_along(upsetSelection)) {
        grpCells[olIdx] <- S4Vectors::metadata(sce)$SOM_stats[as.integer(outputList[[upsetSelection[olIdx]]]), "n"] %>% sum()
    }
    names(ncells) <- ComplexHeatmap::comb_name(cm)
    names(grpCells) <- upsetSelection
    grpSoms <- lapply(outputList[upsetSelection], length) %>% unlist()
    top_ha <- ComplexHeatmap::HeatmapAnnotation(
        "cell #" = ComplexHeatmap::anno_barplot(ncells,
                                                add_numbers = TRUE,
                                                gp = grid::gpar(fill = "black"), width = grid::unit(4, "cm")
        ),
        "som #" = ComplexHeatmap::anno_barplot(ComplexHeatmap::comb_size(cm),
                                               add_numbers = TRUE,
                                               gp = grid::gpar(fill = "black"), width = grid::unit(4, "cm")
        ),
        gap = grid::unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0
    )

    side_ha <- ComplexHeatmap::rowAnnotation(
        "som #" = ComplexHeatmap::anno_barplot(grpSoms,
                                               add_numbers = TRUE,
                                               gp = grid::gpar(fill = NULL), width = grid::unit(3, "cm")
        ),
        "cell #" = ComplexHeatmap::anno_barplot(grpCells,
                                                add_numbers = TRUE,
                                                gp = grid::gpar(fill = NULL), width = grid::unit(3, "cm")
        ),
        gap = grid::unit(2, "mm"), annotation_name_rot = 0
    )
    ComplexHeatmap::UpSet(cm,
                          comb_order = order(
                              ComplexHeatmap::comb_degree(cm),
                              -ComplexHeatmap::comb_size(cm)
                          ),
                          top_annotation = top_ha,
                          right_annotation = side_ha
    )
}
