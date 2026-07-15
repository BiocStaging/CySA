# CySA: Interactive Cluster Selector for Cytometry Data.
# Derived from the clusterSelector Shiny module originally developed in CyDa.
# Refactored for Bioconductor with assistance from the opencode AI coding
# assistant. All code is redistributed under the package LICENSE.

# ui-helpers.R ----
# UI construction helpers for clusterSelector(). Splitting the dashboard into
# small builder functions makes the main factory easier to read and test.

#' Build clusterSelector Sidebar UI
#'
#' @param sce Full \code{SingleCellExperiment} for marker choices.
#' @param dList List of marker pairs used to label the plot index selector.
#' @param outputList Named list of cluster groupings.
#' @param nPlots Number of 2D SOM plots.
#'
#' @return A \code{shinydashboard::dashboardSidebar} UI object.
#'
#' @keywords internal
.buildClusterSelectorSidebar <- function(sce, dList, outputList, nPlots) {
    choices <- lapply(dList, FUN = function(x) paste(x, collapse = " - ")) %>%
        unlist() %>%
        unname()

    shinydashboard::dashboardSidebar(
        shinyjs::useShinyjs(),
        shinyjs::extendShinyjs(
            text = "shinyjs.closeWindow = function() { window.close(); }",
            functions = c("closeWindow")
        ),
        shinydashboard::sidebarMenu(
            id = "tabs",
            shinydashboard::menuItem("Parameters to adjust", tabName = "parameters", icon = shiny::icon("line-chart"))
        ),
        shiny::hr(),
        shiny::conditionalPanel(
            "input.tabs == 'parameters'",
            shiny::fluidRow(
                shiny::column(1),
                shiny::column(
                    10,
                    htmltools::tags$p("The node you most recently clicked:"),
                    shiny::verbatimTextOutput("str")
                ),
                shinydashboardPlus::box(
                    title = "dim modifications", solidHeader = TRUE, width = 12,
                    collapsible = TRUE, collapsed = FALSE,
                    shiny::uiOutput("dimUI")
                ),
                shiny::selectizeInput(
                    inputId = "d2Axes",
                    label = "Choose a plot index :",
                    choices = choices,
                    selected = choices[seq_len(min(6, length(choices)))],
                    multiple = TRUE,
                    size = 6
                ),
                shiny::actionButton("applyDimSelection", "apply dim selection"),
                shiny::radioButtons(
                    "selectMode",
                    "select mode",
                    choices = c("view", "remove others", "add", "remove"),
                    selected = "view",
                    inline = FALSE
                ),
                shiny::selectInput("samples2plot",
                    "samples to plot",
                    choices = unique(as.character(sce$sample_id)),
                    selected = unique(as.character(sce$sample_id)),
                    multiple = TRUE, selectize = TRUE
                ),
                shiny::sliderInput("scatterPercentile",
                    "Scatter auto-zoom percentile",
                    min = 0.5, max = 1, value = 0.99, step = 0.01
                ),
                shiny::selectInput("somColorVar",
                    "SOM 2D color by",
                    choices = c("n", "mean", "median", "rdQu", "max"),
                    selected = "n", multiple = FALSE
                ),
                shiny::selectInput("somSizeVar",
                    "SOM 2D size by",
                    choices = c("n", "mean", "median", "rdQu", "max"),
                    selected = "max", multiple = FALSE
                ),
                shiny::textInput("clusterNumbers",
                    "cluster numbers",
                    value = "1",
                    width = NULL,
                    placeholder = NULL
                ),
                shiny::actionButton("applyclusterNumbers", "apply cluster numbers"),
                shiny::textInput("clusterName",
                    "name selection",
                    value = "",
                    width = NULL,
                    placeholder = NULL
                ),
                shiny::actionButton("applyName", "apply name"),
                shiny::selectInput("groupRM",
                    "select to remove",
                    choices = names(outputList),
                    multiple = TRUE, selectize = TRUE
                ),
                shiny::actionButton("rmGroups", "remove groups"),
                shiny::selectInput("clusterNameSelect",
                    "select named",
                    choices = names(outputList),
                    multiple = TRUE, selectize = TRUE
                ),
                shiny::selectInput("clusterNameRM",
                    "select to remove",
                    choices = names(outputList),
                    multiple = FALSE, selectize = TRUE
                ),
                shiny::actionButton("rmGrp", "remove"),
                shiny::downloadButton("downloadPlots", "Download Plots"),
                shiny::actionButton("close", "Close window")
            )
        )
    )
}


#' Build the First Dashboard Body Row
#'
#' @param colTree Optional collapsible tree object.
#' @param nPlots Number of 2D SOM plots (kept for signature consistency).
#'
#' @return A \code{shiny::fluidRow} object.
#'
#' @keywords internal
.buildFirstBodyRow <- function(colTree, nPlots) {
    shiny::fluidRow(
        if (!is.null(colTree)) {
            shiny::column(
                width = 3,
                shinydashboardPlus::box(
                    title = "interactive tree", solidHeader = TRUE, width = 12, status = "primary",
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
                plotly::plotlyOutput("scatter") %>% shinyjqui::jqui_resizable()
            )
        ),
        shiny::column(
            width = 3,
            shinydashboardPlus::box(
                title = "dendrogram", solidHeader = TRUE, width = 12, status = "primary",
                collapsible = TRUE, collapsed = TRUE,
                shiny::plotOutput("dend") %>% shinyjqui::jqui_resizable()
            ),
            shinydashboardPlus::box(
                title = "interactive dendrogram", solidHeader = TRUE, width = 12, status = "primary",
                collapsible = TRUE, collapsed = TRUE,
                plotly::plotlyOutput("dendPlotly") %>% shinyjqui::jqui_resizable()
            )
        )
    )
}


#' Build SOM 2D Plots UI
#'
#' @param nPlots Number of 2D SOM plots to display.
#' @param colsUsed Character vector of SOM column names used for dim red.
#' @param outputList Named list of cluster groupings.
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildSOM2DPlotsBox <- function(nPlots, colsUsed, outputList) {
    cn <- setdiff(colsUsed, c("label", "clusterid"))

    first_row <- shiny::fluidRow(
        shiny::column(
            width = 4,
            plotly::plotlyOutput("somData1") %>% shinyjqui::jqui_resizable()
        ),
        shiny::column(
            width = 4,
            plotly::plotlyOutput("somData2") %>% shinyjqui::jqui_resizable()
        ),
        shiny::column(
            width = 4,
            plotly::plotlyOutput("somData3") %>% shinyjqui::jqui_resizable()
        )
    )

    second_row <- if (nPlots > 3) {
        shiny::fluidRow(
            if (nPlots > 3) {
                shiny::column(
                    width = 4,
                    plotly::plotlyOutput("somData4") %>% shinyjqui::jqui_resizable()
                )
            },
            if (nPlots > 4) {
                shiny::column(
                    width = 4,
                    plotly::plotlyOutput("somData5") %>% shinyjqui::jqui_resizable()
                )
            },
            if (nPlots > 5) {
                shiny::column(
                    width = 4,
                    plotly::plotlyOutput("somData6") %>% shinyjqui::jqui_resizable()
                )
            }
        )
    } else {
        NULL
    }

    controls <- shiny::fluidRow(
        shiny::column(
            width = 8,
            shiny::selectInput("colorbyGroups",
                "Select groups to color by",
                choices = names(outputList),
                multiple = TRUE, selectize = TRUE
            ),
            shiny::selectInput("dimRedSelection",
                "Select markers to use for dim. Red.",
                choices = cn,
                selected = cn,
                multiple = TRUE, selectize = TRUE
            )
        ),
        shiny::column(
            width = 2,
            shiny::numericInput("perplexity",
                "Perplexity",
                value = 30, min = 1, max = 500
            ),
            shiny::checkboxInput("showlegend", "show legend", value = FALSE)
        ),
        shiny::column(
            width = 2,
            shiny::numericInput("n_neighbors",
                "n_neighbors",
                value = 4, min = 2, max = 500
            )
        )
    )

    dim_red_row <- shiny::fluidRow(
        shiny::column(
            width = 4,
            plotly::plotlyOutput("tsne") %>% shinyjqui::jqui_resizable()
        ),
        shiny::column(
            width = 4,
            plotly::plotlyOutput("umap") %>% shinyjqui::jqui_resizable()
        ),
        shiny::column(
            width = 4,
            plotly::plotlyOutput("pca") %>% shinyjqui::jqui_resizable()
        )
    )

    shinydashboardPlus::box(
        title = "som 2D plots", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(shiny::column(
            width = 2,
            shiny::checkboxInput("showGroups", "color groups", value = FALSE)
        )),
        first_row,
        second_row,
        controls,
        dim_red_row
    )
}


#' Build FlowSOM Marker Pies UI
#'
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildFlowSOMPieBox <- function() {
    shinydashboardPlus::box(
        title = "FlowSOM marker pies", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(shiny::column(
            width = 12,
            htmltools::tags$p(
                "Each pie corresponds to one currently selected SOM node. Slices represent markers, coloured consistently across nodes, with slice size showing the node's mean marker expression."
            ),
            shiny::plotOutput("flowSOMPie") %>% shinyjqui::jqui_resizable()
        ))
    )
}


#' Build FlowSOM Stars UI
#'
#' @param fsom Optional \code{FlowSOM} object. If \code{NULL}, returns
#'   \code{NULL}.
#'
#' @return A \code{shinydashboardPlus::box} UI object or \code{NULL}.
#'
#' @keywords internal
.buildFlowSOMStarsBox <- function(fsom) {
    if (is.null(fsom)) return(NULL)

    shinydashboardPlus::box(
        title = "FlowSOM stars", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(shiny::column(
            width = 12,
            htmltools::tags$p(
                "Interactive FlowSOM star plot. Each node shows marker expression as a star, background colour shows the metaclustering. Click or lasso nodes to select clusters."
            ),
            plotly::plotlyOutput("flowSOMStars") %>% shinyjqui::jqui_resizable()
        ))
    )
}


#' Build Stats UI
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
                        shiny::selectInput("compareStatsTo",
                            "select comparison Stats",
                            choices = num_choices,
                            selected = "none",
                            multiple = FALSE, selectize = TRUE
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::selectInput("relativeTo",
                            "Parent population for percentages",
                            choices = num_choices,
                            selected = "none",
                            multiple = FALSE, selectize = TRUE
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::numericInput("singleNode",
                            "single SOM node counts",
                            min = 1,
                            max = nrow(metaD[[somCodesName]]),
                            value = 1, step = 1
                        )
                    )
                ),
                shiny::fluidRow(
                    shiny::column(
                        width = 4,
                        shiny::selectInput("groupsVar",
                            "Select groups variable for t-test",
                            choices = fact_choices,
                            selected = "none",
                            multiple = FALSE, selectize = TRUE
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::selectInput("group1",
                            "group 1",
                            choices = c(),
                            selected = "",
                            multiple = TRUE, selectize = TRUE
                        )
                    ),
                    shiny::column(
                        width = 4,
                        shiny::selectInput("group2",
                            "group 2",
                            choices = c(),
                            selected = "",
                            multiple = TRUE, selectize = TRUE
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
#' @return A \code{shinydashboardPlus::box} UI object.
#'
#' @keywords internal
.buildSOMRasterBox <- function() {
    shinydashboardPlus::box(
        title = "som plots", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(shiny::column(
            width = 12,
            shiny::plotOutput("somRaster", height = "1200px") %>% shinyjqui::jqui_resizable()
        )),
        shiny::fluidRow(shiny::column(
            width = 3,
            plotly::plotlyOutput("somRasterSelect", height = "400px") %>% shinyjqui::jqui_resizable()
        ))
    )
}


#' Build Violin Plots UI
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
            shiny::selectInput("violinSelection",
                "Select markers to show",
                choices = colnames(somCodes),
                selected = colsUsed,
                multiple = TRUE, selectize = TRUE
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
            shiny::selectInput("upsetSelection",
                "Select groups to show",
                choices = names(outputList),
                selected = default,
                multiple = TRUE, selectize = TRUE
            )
        ),
        shiny::fluidRow(shiny::column(
            width = 12,
            shiny::plotOutput("UpSet") %>% shinyjqui::jqui_resizable()
        ))
    )
}
