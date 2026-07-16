# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# server-helpers.R ----
# Non-reactive helper functions used inside the clusterSelector server. Moving
# these out of the server closure reduces the size of clusterSelector() and
# makes each piece independently testable.

#' Synchronize Select Input Choices with outputList
#'
#' Updates the UI select inputs that depend on the names of \code{outputList}.
#'
#' @param session Shiny session object.
#' @param input Shiny input object.
#' @param outputList Named list of cluster groupings.
#' @param metaD \code{metadata(sce)} list.
#'
#' @keywords internal
.updateOutputListInputs <- function(session, input, outputList, metaD) {
    ol_names <- names(outputList)
    oldVal <- isolate(input$clusterNameSelect)
    shiny::updateSelectInput(session = session, "clusterNameSelect", choices = ol_names, selected = oldVal)

    oldval <- isolate(input$compareStatsTo)
    numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
    expInfo <- metaD$experiment_info[, numCols, drop = FALSE]
    eI <- apply(expInfo, 2, as.numeric)
    if (any(is.na(eI))) {
        stop("NAs produced when converting experiment_info to numeric")
    }
    choices <- c("none", colnames(eI), ol_names)
    shiny::updateSelectInput(session = session, "compareStatsTo", choices = choices, selected = oldval)

    oldval <- isolate(input$relativeTo)
    shiny::updateSelectInput(session = session, "relativeTo", choices = choices, selected = oldval)

    oldval <- isolate(input$upsetSelection)
    new_sel <- union(oldval, intersect(c("Rest", "selected"), ol_names))
    if (length(new_sel) == 0) new_sel <- ol_names
    shiny::updateSelectInput(session = session, "upsetSelection", choices = ol_names, selected = new_sel)

    oldval <- isolate(input$colorbyGroups)
    shiny::updateSelectInput(session = session, "colorbyGroups", choices = ol_names, selected = oldval)

    oldval <- isolate(input$groupRM)
    shiny::updateSelectInput(session = session, "groupRM", choices = ol_names, selected = oldval)

    shiny::updateSelectInput(session = session, "clusterNameRM", choices = ol_names)
}


#' Compute Relative Cell Counts
#'
#' Computes the percentage of selected cells relative to a reference column,
#' another cluster group, or the total population.
#'
#' @param clusterPatientTable Table of sample by cluster counts.
#' @param rs Selected cluster ids.
#' @param relativeToCol Reference column or group name.
#' @param expInfo Sample-level metadata data.frame.
#' @param outputList Named list of cluster groupings.
#'
#' @return Numeric vector of percentages.
#'
#' @keywords internal
.computeRelativeCounts <- function(clusterPatientTable, rs, relativeToCol, expInfo, outputList) {
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


#' Update outputList After Cluster Assignment Change
#'
#' Rebuilds \code{outputList} after a group is removed or renamed, ensuring the
#' \code{"Rest"} group is always correct.
#'
#' @param outputList Current named list of cluster groupings.
#' @param sceLevels Character or integer vector of all cluster levels.
#' @param removed Optional vector of cluster ids that should be excluded.
#'
#' @return Updated named list.
#'
#' @keywords internal
.rebuildOutputList <- function(outputList, sceLevels, removed = NULL) {
    if (!is.null(removed)) {
        for (nm in names(outputList)) {
            outputList[[nm]] <- setdiff(outputList[[nm]], removed)
        }
    }

    used <- unique(unlist(outputList))
    outputList[["Rest"]] <- as.integer(sceLevels[!sceLevels %in% used])

    for (na in names(outputList)) {
        if (length(outputList[[na]]) == 0) outputList[[na]] <- NULL
    }
    outputList
}

# inputSelect = function ----
.inputSelect <- function(d, rs, mode) {
    if (is.null(rs)) return(integer(0))
    if (is.null(d))  return(as.integer(rs))

    node_ids <- NULL

    if (is.data.frame(d)) {
        if ("key" %in% names(d) && length(d$key) > 0L && !all(is.na(d$key))) {
            # ggplotly key aesthetic — works across all color-group traces
            node_ids <- as.integer(unlist(d$key))

        } else if ("customdata" %in% names(d) && !all(is.na(d$customdata))) {
            # injected via plotly::style() on SOM/dend plots
            node_ids <- as.integer(unlist(d$customdata))

        } else if ("pointNumber" %in% names(d)) {
            # last resort: row-order encoding, only valid when data is node-ordered
            # filter to curveNumber == 0 (the all-nodes base layer)
            if ("curveNumber" %in% names(d)) {
                d <- d[d$curveNumber == 0L, , drop = FALSE]
                if (nrow(d) == 0L) return(as.integer(rs))
            }
            node_ids <- as.integer(unlist(d$pointNumber)) + 1L
        }
    } else {
        node_ids <- suppressWarnings(as.integer(unlist(d)))
    }

    if (is.null(node_ids)) return(as.integer(rs))
    node_ids <- node_ids[!is.na(node_ids) & node_ids > 0L]
    if (length(node_ids) == 0L) return(as.integer(rs))

    as.integer(switch(EXPR = mode,
                      "remove others" = intersect(node_ids, rs),
                      "add"           = unique(c(rs, node_ids)),
                      "remove"        = setdiff(rs, node_ids),
                      node_ids
    ))
}


