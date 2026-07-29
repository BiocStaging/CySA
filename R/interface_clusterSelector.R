# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# interface_clusterSelector.R ----
# UI construction helpers for clusterSelector(). Splitting the dashboard into
# small builder functions makes the main factory easier to read and test.


# Central tooltip texts ---------------------------------------------------------
# Named character vector used to keep hover-over help consistent across all UI
# builder functions.  Names correspond to Shiny inputIds so that adding a new
# control only requires adding one entry here.
.clusterSelectorTooltips <- c(
    selectMode = paste(
        "Choose how plot selections change the currently selected SOM nodes:",
        "view = replace; add = union; remove = subtract; remove others = intersect."
    ),
    samples2plot = "Restrict scatter plots and statistics to the selected samples.",
    scatterPercentile = paste(
        "Percentile range used to auto-zoom the scatter plot axes around the",
        "selected cells."
    ),
    somColorVar = "Statistic used to colour SOM nodes in the 2D plots.",
    somSizeVar = "Statistic used to scale SOM node points in the 2D plots.",
    clusterNumbers = "Comma-separated SOM node ids to add to the current selection.",
    applyclusterNumbers = "Create or update the selection from the node ids above.",
    clusterName = "Name for the currently selected group of nodes.",
    applyName = "Save the current selection under the name above.",
    groupRM = "Existing named groups to delete.",
    rmGroups = "Delete the selected named groups.",
    clusterNameSelect = "Existing groups to activate or view together.",
    clusterNameRM = "Single existing group to remove.",
    rmGrp = "Remove the named group selected above.",
    downloadPlots = "Download the current plots as a PDF.",
    close = "Close the clusterSelector window.",
    currentDimX = "Marker/channel for the X axis of the main SOM 2D plot.",
    currentDimY = "Marker/channel for the Y axis of the main SOM 2D plot.",
    dimPairSelect = "Preset marker pair; choosing it updates both axes.",
    showGroups = "Colour nodes by saved group membership instead of by statistic.",
    colorbyGroups = "Saved groups used for node colouring.",
    dimRedSelection = "Markers used for t-SNE, UMAP and PCA calculations.",
    perplexity = "t-SNE perplexity; lower values focus on local structure.",
    n_neighbors = "UMAP number of neighbours; affects local/global balance.",
    showlegend = "Show a legend on the dimension-reduction plots.",
    staticSomX = "Marker for the X axis of a new static SOM pair view.",
    staticSomY = "Marker for the Y axis of a new static SOM pair view.",
    staticSomAdd = "Add the chosen pair to the static grid.",
    staticSomRemove = "Existing preset pair to remove from the grid.",
    staticSomRemoveBtn = "Delete the selected preset pair.",
    staticSomCols = "Number of plots per row in the static grid.",
    staticSomHeight = "Height in pixels for each static SOM plot.",
    flowSOMPieCols = "Number of pie-chart columns in the FlowSOM view.",
    flowSOMPieSize = "Pixel size of each FlowSOM pie.",
    flowSOMPieMax = "Maximum number of pies to display in the FlowSOM view.",
    compareStatsTo = "Numerical column or saved group to compare cell counts against.",
    relativeTo = "Reference population used for percentage calculations.",
    singleNode = "Show counts for one specific SOM node alongside the selection.",
    groupsVar = "Sample-level factor used to split samples for the t-test.",
    group1 = "Samples belonging to group 1 for the t-test.",
    group2 = "Samples belonging to group 2 for the t-test.",
    violinSelection = "Markers displayed in the violin plots.",
    upsetSelection = "Saved groups compared in the UpSet plot."
)


#' Build clusterSelector Sidebar UI
#'
#' Constructs the left-hand control panel: selection mode, sample filter,
#' scatter-plot auto-zoom, SOM 2D plot aesthetics, and group management.
#'
#' Axis controls have moved to the SOM 2D plots body box to avoid the
#' selectize-dropdown overflow that occurs inside the narrow sidebar.
#'
#' @param sce Full \code{SingleCellExperiment} for sample-id choices.
#' @param outputList Named list of cluster groupings.
#' @param nPlots Kept for API compatibility; no longer controls dim-UI layout.
#'
#' @return A \code{shinydashboard::dashboardSidebar} UI object.
#'
#' @keywords internal
.buildClusterSelectorSidebar <- function(sce, outputList, nPlots) {
    shinydashboard::dashboardSidebar(
        shinyjs::useShinyjs(),
        shinyjs::extendShinyjs(
            text      = "shinyjs.closeWindow = function() { window.close(); }",
            functions = c("closeWindow")
        ),
        shiny::tags$script(shiny::HTML("
            (function() {
                var resizeTimer = null;
                function resetAndResizePlots() {
                    $('.ui-resizable').css({ width: '100%', height: '100%' });
                    $('.js-plotly-plot').each(function() {
                        Plotly.Plots.resize(this);
                    });
                }
                $(window).on('resize', function() {
                    clearTimeout(resizeTimer);
                    resizeTimer = setTimeout(resetAndResizePlots, 200);
                });
            })();
        ")),
        shinydashboard::sidebarMenu(
            id = "tabs",
            shinydashboard::menuItem(
                "Parameters", tabName = "parameters",
                icon = shiny::icon("sliders")
            )
        ),
        shiny::hr(),
        shiny::conditionalPanel(
            "input.tabs == 'parameters'",
            shiny::div(
                style = "padding: 0 12px;",

                # Selection mode and sample filter ----
                shiny::div(
                    title = .clusterSelectorTooltips[["selectMode"]],
                    shiny::radioButtons(
                        "selectMode", "Select mode",
                        choices  = c("view", "remove others", "add", "remove"),
                        selected = "view",
                        inline   = FALSE
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["samples2plot"]],
                    shiny::selectInput(
                        "samples2plot", "Samples to plot",
                        choices  = unique(as.character(sce$sample_id)),
                        selected = unique(as.character(sce$sample_id)),
                        multiple = TRUE, selectize = TRUE
                    )
                ),

                # Scatter plot auto-zoom ----
                shiny::div(
                    title = .clusterSelectorTooltips[["scatterPercentile"]],
                    shiny::sliderInput(
                        "scatterPercentile", "Scatter auto-zoom percentile",
                        min = 0.5, max = 1, value = 0.99, step = 0.01
                    )
                ),
                shiny::helpText(
                    "Lower percentiles zoom in further on dense regions."
                ),

                # SOM 2D plot aesthetics ----
                shiny::div(
                    title = .clusterSelectorTooltips[["somColorVar"]],
                    shiny::selectInput(
                        "somColorVar", "SOM 2D color by",
                        choices  = c("n", "mean", "median", "rdQu", "max"),
                        selected = "n", multiple = FALSE
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["somSizeVar"]],
                    shiny::selectInput(
                        "somSizeVar", "SOM 2D size by",
                        choices  = c("n", "mean", "median", "rdQu", "max"),
                        selected = "max", multiple = FALSE
                    )
                ),
                shiny::helpText(
                    "Color and size summarise marker expression per SOM node."
                ),

                shiny::hr(),

                # Selection by node ids ----
                shiny::div(
                    title = .clusterSelectorTooltips[["clusterNumbers"]],
                    shiny::textInput(
                        "clusterNumbers", "Cluster numbers",
                        value = "1"
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["applyclusterNumbers"]],
                    shiny::actionButton("applyclusterNumbers", "Apply cluster numbers")
                ),

                shiny::hr(),

                # Save current selection as named group ----
                shiny::div(
                    title = .clusterSelectorTooltips[["clusterName"]],
                    shiny::textInput(
                        "clusterName", "Name selection",
                        value = ""
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["applyName"]],
                    shiny::actionButton("applyName", "Apply name")
                ),

                shiny::hr(),

                # Remove named groups ----
                shiny::div(
                    title = .clusterSelectorTooltips[["groupRM"]],
                    shiny::selectInput(
                        "groupRM", "Select to remove",
                        choices  = names(outputList),
                        multiple = TRUE, selectize = TRUE
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["rmGroups"]],
                    shiny::actionButton("rmGroups", "Remove groups")
                ),

                shiny::hr(),

                # Remove a single named group ----
                shiny::div(
                    title = .clusterSelectorTooltips[["clusterNameSelect"]],
                    shiny::selectInput(
                        "clusterNameSelect", "Select named",
                        choices  = names(outputList),
                        multiple = TRUE, selectize = TRUE
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["clusterNameRM"]],
                    shiny::selectInput(
                        "clusterNameRM", "Remove by name",
                        choices  = names(outputList),
                        multiple = FALSE, selectize = TRUE
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["rmGrp"]],
                    shiny::actionButton("rmGrp", "Remove")
                ),

                shiny::hr(),

                # Window controls ----
                shiny::div(
                    title = .clusterSelectorTooltips[["downloadPlots"]],
                    shiny::downloadButton("downloadPlots", "Download Plots")
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["close"]],
                    shiny::actionButton("close", "Close window")
                )
            )
        )
    )
}


#' Build the First Dashboard Body Row
#'
#' Top row of the dashboard: collapsible sample tree (optional), the main 2D
#' scatter plot, and the interactive dendrogram. All three outputs support
#' node/point selection.
#'
#' @param colTree Optional collapsible tree object.
#' @param nPlots Number of 2D SOM plots (kept for signature consistency).
#'
#' @return A \code{shiny::fluidRow} object.
#'
#' @keywords internal
.buildFirstBodyRow <- function(colTree, nPlots) {
    selectionHint <- htmltools::tags$p(
        "Tip: click or lasso points/nodes to select clusters. Selection mode is",
        "controlled from the sidebar.",
        style = "color:#888; font-style:italic; font-size:12px; margin:4px 0;"
    )

    shiny::fluidRow(
        if (!is.null(colTree)) {
            shiny::column(
                width = 3,
                shinydashboardPlus::box(
                    title = "interactive tree", solidHeader = TRUE,
                    width = 12, status = "primary",
                    collapsible = TRUE, collapsed = TRUE,
                    collapsibleTree::collapsibleTreeOutput("plot")
                )
            )
        },
        shiny::column(
            width = 6,
            shinydashboardPlus::box(
                title = "2D plot", solidHeader = TRUE, width = 12, status = "primary",
                collapsible = TRUE, collapsed = TRUE,
                selectionHint,
                plotly::plotlyOutput("scatter") %>% shinyjqui::jqui_resizable()
            )
        ),
        shiny::column(
            width = 6,
            shinydashboardPlus::box(
                title = "interactive dendrogram", solidHeader = TRUE,
                width = 12, status = "primary",
                collapsible = TRUE, collapsed = TRUE,
                selectionHint,
                plotly::plotlyOutput("dendPlotly") %>% shinyjqui::jqui_resizable()
            )
        )
    )
}


#' Build SOM 2D Plot Box
#'
#' Main SOM visualisation panel: interactive axis controls, the primary 2D SOM
#' plot, dimension-reduction controls, and t-SNE/UMAP/PCA projection plots.
#' All plot outputs are selection-enabled.
#'
#' @param nPlots Kept for API compatibility.
#' @param colsUsed Character vector of SOM column names for dim-red selector.
#' @param outputList Named list of cluster groupings.
#' @param markers Character vector of all available markers (\code{rownames(sce)}).
#' @param dList List of marker pairs used to populate the preset-pair picker.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildSOM2DPlotsBox <- function(nPlots, colsUsed, outputList, markers, dList) {
    cn     <- setdiff(colsUsed, c("label", "clusterid"))
    init_x <- if (length(dList) > 0L) dList[[1L]][1L] else markers[1L]
    init_y <- if (length(dList) > 0L) dList[[1L]][2L] else markers[min(2L, length(markers))]

    selectionHint <- htmltools::tags$p(
        "Tip: click or lasso points to select SOM nodes. Selection mode is",
        "controlled from the sidebar.",
        style = "color:#888; font-style:italic; font-size:12px; margin:4px 0;"
    )

    shinydashboardPlus::box(
        title = "SOM 2D plots", solidHeader = TRUE,
        width = 12L, status = "primary",
        collapsible = TRUE, collapsed = TRUE,

        ## -- axis controls (full body width) ----------------------------------
        shiny::fluidRow(
            shiny::column(
                width = 4L,
                shiny::div(
                    title = .clusterSelectorTooltips[["currentDimX"]],
                    shiny::selectInput(
                        "currentDimX", "X axis",
                        choices = markers, selected = init_x,
                        multiple = FALSE, selectize = TRUE, width = "100%"
                    )
                )
            ),
            shiny::column(
                width = 4L,
                shiny::div(
                    title = .clusterSelectorTooltips[["currentDimY"]],
                    shiny::selectInput(
                        "currentDimY", "Y axis",
                        choices = markers, selected = init_y,
                        multiple = FALSE, selectize = TRUE, width = "100%"
                    )
                )
            ),
            shiny::column(
                width = 4L,
                ## Rendered by server so it stays in sync with dListRV
                shiny::div(
                    title = .clusterSelectorTooltips[["dimPairSelect"]],
                    shiny::uiOutput("dimPairSelectUI")
                )
            )
        ),

        ## -- show-groups + single interactive plot -----------------------------
        shiny::fluidRow(shiny::column(
            width = 2L,
            shiny::div(
                title = .clusterSelectorTooltips[["showGroups"]],
                shiny::checkboxInput(
                    "showGroups", "Colour groups", value = FALSE
                )
            )
        )),
        shiny::fluidRow(shiny::column(
            width = 12L,
            selectionHint,
            plotly::plotlyOutput("somDataMain") %>% shinyjqui::jqui_resizable()
        )),

        ## -- dim-red controls -------------------------------------------------
        shiny::fluidRow(
            shiny::column(
                width = 8L,
                shiny::div(
                    title = .clusterSelectorTooltips[["colorbyGroups"]],
                    shiny::selectInput(
                        "colorbyGroups", "Colour by groups",
                        choices = names(outputList), multiple = TRUE, selectize = TRUE
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["dimRedSelection"]],
                    shiny::selectInput(
                        "dimRedSelection", "Markers for dim. reduction",
                        choices = cn, selected = cn, multiple = TRUE, selectize = TRUE
                    )
                ),
                shiny::helpText(
                    "Dimension-reduction runs on SOM code vectors using the selected markers."
                )
            ),
            shiny::column(
                width = 2L,
                shiny::div(
                    title = .clusterSelectorTooltips[["perplexity"]],
                    shiny::numericInput(
                        "perplexity", "Perplexity",
                        value = 30, min = 1, max = 500
                    )
                ),
                shiny::div(
                    title = .clusterSelectorTooltips[["showlegend"]],
                    shiny::checkboxInput(
                        "showlegend", "Show legend", value = FALSE
                    )
                )
            ),
            shiny::column(
                width = 2L,
                shiny::div(
                    title = .clusterSelectorTooltips[["n_neighbors"]],
                    shiny::numericInput(
                        "n_neighbors", "n_neighbors",
                        value = 4, min = 2, max = 500
                    )
                )
            )
        ),

        ## -- dim-red projection plots ------------------------------------------
        shiny::fluidRow(
            shiny::column(width = 4L,
                          selectionHint,
                          plotly::plotlyOutput("tsne") %>% shinyjqui::jqui_resizable()),
            shiny::column(width = 4L,
                          selectionHint,
                          plotly::plotlyOutput("umap") %>% shinyjqui::jqui_resizable()),
            shiny::column(width = 4L,
                          selectionHint,
                          plotly::plotlyOutput("pca")  %>% shinyjqui::jqui_resizable())
        )
    )
}


#' Build Six Static SOM Views Box
#'
#' Grid of static ggplot SOM pair views. Users can add/remove marker pairs,
#' change the grid layout, and adjust per-plot height.
#'
#' @param markers Character vector of all available markers (\code{rownames(sce)}).
#' @param dList Initial list of marker pairs.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#' @keywords internal
.buildSixStaticSOMBox <- function(markers, dList) {
    ## Initial remove-picker choices
    if (length(dList) > 0L) {
        pair_labels    <- vapply(dList,
                                 function(p) paste(p[1L], "-", p[2L]), character(1L))
        remove_choices <- stats::setNames(seq_along(dList), pair_labels)
    } else {
        remove_choices <- character(0L)
    }

    init_x <- if (length(dList) > 0L) dList[[1L]][1L] else markers[1L]
    init_y <- if (length(dList) > 0L) dList[[1L]][2L] else markers[min(2L, length(markers))]

    shinydashboardPlus::box(
        title       = "All SOM pair views (static)",
        solidHeader = TRUE, width = 12L, status = "primary",
        collapsible = TRUE, collapsed = TRUE,

        ## -- pair management row ----------------------------------------------
        shiny::fluidRow(
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["staticSomX"]],
                    shiny::selectInput(
                        "staticSomX", "X axis (add)",
                        choices = markers, selected = init_x,
                        multiple = FALSE, selectize = TRUE, width = "100%"
                    )
                )
            ),
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["staticSomY"]],
                    shiny::selectInput(
                        "staticSomY", "Y axis (add)",
                        choices = markers, selected = init_y,
                        multiple = FALSE, selectize = TRUE, width = "100%"
                    )
                )
            ),
            shiny::column(
                width = 2L,
                shiny::div(
                    style = "margin-top:25px;",
                    title = .clusterSelectorTooltips[["staticSomAdd"]],
                    shiny::actionButton(
                        "staticSomAdd", label = "Add",
                        icon  = shiny::icon("plus"),
                        class = "btn-success btn-sm"
                    )
                )
            ),
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["staticSomRemove"]],
                    shiny::selectInput(
                        "staticSomRemove", "Remove pair",
                        choices   = remove_choices,
                        multiple  = FALSE,
                        selectize = FALSE,   # native <select>: never overflows
                        width     = "100%"
                    )
                )
            ),
            shiny::column(
                width = 1L,
                shiny::div(
                    style = "margin-top:25px;",
                    title = .clusterSelectorTooltips[["staticSomRemoveBtn"]],
                    shiny::actionButton(
                        "staticSomRemoveBtn", label = NULL,
                        icon  = shiny::icon("trash"),
                        class = "btn-danger btn-sm"
                    )
                )
            )
        ),

        ## -- layout controls --------------------------------------------------
        shiny::fluidRow(
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["staticSomCols"]],
                    shiny::selectInput(
                        "staticSomCols", "Columns per row",
                        choices  = c("1" = 1L, "2" = 2L, "3" = 3L,
                                     "4" = 4L, "6" = 6L),
                        selected = 3L,
                        multiple = FALSE, selectize = FALSE, width = "100%"
                    )
                )
            ),
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["staticSomHeight"]],
                    shiny::sliderInput(
                        "staticSomHeight", "Plot height (px)",
                        min = 150L, max = 600L, value = 240L, step = 10L
                    )
                )
            )
        ),

        shiny::hr(),
        htmltools::tags$p(
            "Red dots = selected SOM nodes. Expand box to render.",
            style = paste(
                "color:#888; font-style:italic;",
                "font-size:12px; margin-bottom:8px;"
            )
        ),

        ## -- dynamic plot grid (filled by renderUI in server) -----------------
        shiny::uiOutput("staticSomGrid")
    )
}

#' Build FlowSOM Marker Pies UI
#'
#' Displays each selected SOM node as a pie chart coloured by mean marker
#' expression. Users control the grid layout and pie size.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#' @keywords internal
.buildFlowSOMPieBox <- function() {
    shinydashboardPlus::box(
        title       = "FlowSOM marker pies",
        solidHeader = TRUE, width = 12L, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["flowSOMPieCols"]],
                    shiny::numericInput(
                        inputId = "flowSOMPieCols",
                        label   = "Columns",
                        value   = 5L, min = 1L, max = 20L, step = 1L,
                        width   = "100%"
                    )
                )
            ),
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["flowSOMPieSize"]],
                    shiny::sliderInput(
                        inputId = "flowSOMPieSize",
                        label   = "Size per pie (px)",
                        min     = 50L, max = 500L, value = 150L, step = 10L,
                        width   = "100%"
                    )
                )
            ),
            shiny::column(
                width = 3L,
                shiny::div(
                    title = .clusterSelectorTooltips[["flowSOMPieMax"]],
                    shiny::numericInput(
                        inputId = "flowSOMPieMax",
                        label   = "Max pies",
                        value   = 5L, min = 1L, max = 50L, step = 1L,
                        width   = "100%"
                    )
                )
            ),
            shiny::column(
                width = 3L,
                htmltools::tags$p(
                    paste(
                        "Each pie = one selected SOM node.",
                        "Slice size = mean marker expression.",
                        "Colours are consistent across nodes."
                    ),
                    style = "color:#666; font-size:12px; margin-top:8px;"
                )
            )
        ),
        shiny::fluidRow(shiny::column(
            width = 12L,
            ## Height is set dynamically in the server so each pie keeps a
            ## consistent size regardless of how many nodes are selected.
            shiny::uiOutput("flowSOMPieUI")
        ))
    )
}

#' Build FlowSOM Stars UI
#'
#' Interactive FlowSOM star plot. Each node shows marker expression as a star,
#' background colour shows the metaclustering. Click or lasso nodes to select
#' clusters.
#'
#' @param fsom Optional \code{FlowSOM} object. If \code{NULL}, returns
#'   \code{NULL}.
#'
#' @return A \code{shinydashboardPlus::box} UI object or \code{NULL}.
#' @keywords internal
.buildFlowSOMStarsBox <- function(fsom) {
    if (is.null(fsom)) return(NULL)

    shinydashboardPlus::box(
        title = "FlowSOM stars", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(shiny::column(
            width = 12,
            htmltools::tags$p(
                "Interactive FlowSOM star plot. Each node shows marker expression",
                "as a star; background colour shows metaclustering. Click or",
                "lasso nodes to select clusters."
            ),
            plotly::plotlyOutput("flowSOMStars") %>% shinyjqui::jqui_resizable()
        ))
    )
}


#' Build Stats UI
#'
#' Statistics panel: compare cell counts, compute percentages, run sample-group
#' t-tests, and display the results as tables and bar plots.
#'
#' @param metaD \code{metadata(sce)} list used to build numeric/factor choices.
#' @param somCodesName Name of the SOM codes metadata slot.
#' @param outputList Named list of cluster groupings.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildStatsBox <- function(metaD, somCodesName, outputList) {
    numCols <- unlist(lapply(metaD$experiment_info, is.numeric), use.names = FALSE)
    expInfo <- metaD$experiment_info[, numCols, drop = FALSE]
    eI <- apply(expInfo, 2, as.numeric)
    if (any(is.na(eI))) {
        stop("NAs produced when converting experiment_info to numeric")
    }
    num_choices <- c("none", colnames(eI))

    factCols <- unlist(lapply(metaD$experiment_info, is.factor), use.names = FALSE)
    fact_choices <- c("none", colnames(metaD$experiment_info)[factCols])

    shinydashboardPlus::box(
        title = "Stats", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(
            shiny::column(
                width = 4,
                shiny::verbatimTextOutput("somClusters")
            ),
            shiny::column(
                width = 8,
                shiny::fluidRow(
                    shiny::column(
                        width = 4,
                        shiny::div(
                            title = .clusterSelectorTooltips[["compareStatsTo"]],
                            shiny::selectInput("compareStatsTo",
                                "select comparison Stats",
                                choices = num_choices,
                                selected = "none",
                                multiple = FALSE, selectize = TRUE
                            )
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::div(
                            title = .clusterSelectorTooltips[["relativeTo"]],
                            shiny::selectInput("relativeTo",
                                "Parent population for percentages",
                                choices = num_choices,
                                selected = "none",
                                multiple = FALSE, selectize = TRUE
                            )
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::div(
                            title = .clusterSelectorTooltips[["singleNode"]],
                            shiny::numericInput("singleNode",
                                "single SOM node counts",
                                min = 1,
                                max = nrow(metaD[[somCodesName]]),
                                value = 1, step = 1
                            )
                        )
                    )
                ),
                shiny::fluidRow(
                    shiny::column(
                        width = 4,
                        shiny::div(
                            title = .clusterSelectorTooltips[["groupsVar"]],
                            shiny::selectInput("groupsVar",
                                "Select groups variable for t-test",
                                choices = fact_choices,
                                selected = "none",
                                multiple = FALSE, selectize = TRUE
                            )
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::div(
                            title = .clusterSelectorTooltips[["group1"]],
                            shiny::selectInput("group1",
                                "group 1",
                                choices = c(),
                                selected = "",
                                multiple = TRUE, selectize = TRUE
                            )
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::div(
                            title = .clusterSelectorTooltips[["group2"]],
                            shiny::selectInput("group2",
                                "group 2",
                                choices = c(),
                                selected = "",
                                multiple = TRUE, selectize = TRUE
                            )
                        )
                    )
                ),
                shiny::fluidRow(
                    shinydashboardPlus::box(
                        title = "Stats", solidHeader = TRUE, width = 12, status = "primary",
                        collapsible = TRUE, collapsed = TRUE,
                        shiny::fluidRow(
                            DT::DTOutput("cellCounts") %>% shinyjqui::jqui_resizable(),
                            shiny::verbatimTextOutput("cellPercentages") %>% shinyjqui::jqui_resizable(),
                        ),
                        shiny::fluidRow(
                            shiny::column(
                                width = 6,
                                shiny::plotOutput("CountBar") %>% shinyjqui::jqui_resizable()
                            ),
                            shiny::column(
                                width = 6,
                                shiny::plotOutput("PercentBar") %>% shinyjqui::jqui_resizable()
                            )
                        )
                    )
                ),
                shiny::fluidRow(
                    shinydashboardPlus::box(
                        title = "t-Test", solidHeader = TRUE, width = 12, status = "primary",
                        collapsible = TRUE, collapsed = TRUE,
                        shiny::fluidRow(
                            shiny::column(
                                width = 12,
                                shiny::verbatimTextOutput("ttestResult")
                            )
                        )
                    )
                )
            )
        )
    )
}


#' Build SOM Raster UI
#'
#' Raster heatmap of SOM node marker expression plus a small interactive plot
#' for selecting nodes from the SOM grid.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildSOMRasterBox <- function() {
    selectionHint <- htmltools::tags$p(
        "Tip: click or lasso nodes in the right-hand plot to add them to the",
        "current selection.",
        style = "color:#888; font-style:italic; font-size:12px; margin:4px 0;"
    )

    shinydashboardPlus::box(
        title = "som plots", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(shiny::column(
            width = 12,
            shiny::plotOutput("somRaster", height = "1200px") %>% shinyjqui::jqui_resizable()
        )),
        shiny::fluidRow(shiny::column(
            width = 12,
            selectionHint,
            plotly::plotlyOutput("somRasterSelect", height = "400px") %>% shinyjqui::jqui_resizable()
        ))
    )
}


#' Build Violin Plots UI
#'
#' Marker expression violin plots for the saved cluster groups. Users choose
#' which markers are shown.
#'
#' @param somCodes SOM codes matrix; column names become marker choices.
#' @param colsUsed Default marker selection.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildViolinBox <- function(somCodes, colsUsed) {
    shinydashboardPlus::box(
        id = "violinBox",
        title = "violin plots", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(
            shiny::column(width = 12),
            shiny::div(
                title = .clusterSelectorTooltips[["violinSelection"]],
                shiny::selectInput("violinSelection",
                    "Select markers to show",
                    choices = colnames(somCodes),
                    selected = colsUsed,
                    multiple = TRUE, selectize = TRUE
                )
            )
        ),
        shiny::fluidRow(
            shiny::column(
                width = 6,
                shiny::plotOutput("VlnPlot") %>% shinyjqui::jqui_resizable()
            ),
            shiny::column(
                width = 6,
                shiny::plotOutput("VlnPlot2") %>% shinyjqui::jqui_resizable()
            )
        )
    )
}


#' Build UpSet Plot UI
#'
#' UpSet-style overlap visualisation of saved cluster groups. Users choose which
#' groups to compare.
#'
#' @param outputList Named list of cluster groupings.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildUpSetBox <- function(outputList) {
    default <- intersect(c("Rest", "selected"), names(outputList))
    if (length(default) == 0) default <- names(outputList)

    shinydashboardPlus::box(
        id = "upsetBox",
        title = "UpSet plot", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(
            shiny::column(width = 12),
            shiny::div(
                title = .clusterSelectorTooltips[["upsetSelection"]],
                shiny::selectInput("upsetSelection",
                    "Select groups to show",
                    choices = names(outputList),
                    selected = default,
                    multiple = TRUE, selectize = TRUE
                )
            )
        ),
        shiny::fluidRow(shiny::column(
            width = 12,
            shiny::plotOutput("UpSet") %>% shinyjqui::jqui_resizable()
        ))
    )
}
