# Tests for observers_clusterSelector.R
# These test the observer logic in the clusterSelector app.

library(shiny)
library(testthat)

# =============================================================================
# .registerClusterSelectorObservers() - General observers
# =============================================================================
test_that("observers: dimPairSelect updates activePlot", {
  test_app <- make_test_app()
  shiny::testServer(app = test_app$app, expr = {
    suppressWarnings(session$setInputs(dimPairSelect = 2L))
    suppressWarnings(session$flushReact())
    expect_equal(activePlot(), 2L)
  })
})

test_that("observers: clusterNumbers input syncs with rsUsed", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Set cluster numbers via input
      suppressWarnings(session$setInputs(clusterNumbers = "1, 2, 3"))
      suppressWarnings(session$flushReact())

      # rsUsed should be updated
      rs <- rsUsed()
      expect_true(all(c(1L, 2L, 3L) %in% rs))
    }
  )
})


test_that("observers: applyName creates a new named group", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Initial state
      expect_equal(sort(names(rv$outputList)), "Rest")

      # Select clusters and apply a name
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2,3",
        clusterName = "MyCluster",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      # Verify the new group was created
      expect_true("MyCluster" %in% names(rv$outputList))
      expect_equal(sort(rv$outputList$MyCluster), c(1L, 2L, 3L))
    }
  )
})


test_that("observers: applyName does nothing with empty name", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Note: "selected" is auto-added by the sync observer
      initial_names <- sort(names(rv$outputList))

      # Try to apply with empty name
      suppressWarnings(session$setInputs(
        clusterNumbers = "1,2",
        clusterName = "",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      # State should be unchanged (except for "selected" which is auto-added)
      expect_true(all(initial_names %in% sort(names(rv$outputList))))
    }
  )
})


test_that("observers: rmGroups removes multiple groups", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # First create two groups
      suppressWarnings(session$setInputs(
        clusterName = "Group1",
        clusterNumbers = "1,2",
        applyName = 1
      ))
      suppressWarnings(session$flushReact())

      suppressWarnings(session$setInputs(
        clusterName = "Group2",
        clusterNumbers = "3,4",
        applyName = 2
      ))
      suppressWarnings(session$flushReact())

      # Verify groups exist
      expect_true("Group1" %in% names(rv$outputList))
      expect_true("Group2" %in% names(rv$outputList))

      # Select both groups for removal and trigger rmGroups
      suppressWarnings(session$setInputs(
        groupRM = c("Group1", "Group2"),
        rmGroups = 1
      ))
      suppressWarnings(session$flushReact())

      # Verify groups were removed
      expect_false("Group1" %in% names(rv$outputList))
      expect_false("Group2" %in% names(rv$outputList))
    }
  )
})


test_that("observers: rmGroups does nothing when no groups selected", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Note: "selected" is auto-added by the sync observer
      initial_names <- sort(names(rv$outputList))

      # Trigger rmGroups without selecting any groups
      suppressWarnings(session$setInputs(rmGroups = 1))
      suppressWarnings(session$flushReact())

      # State should be unchanged (except for "selected" which is auto-added)
      expect_true(all(initial_names %in% sort(names(rv$outputList))))
    }
  )
})


# test_that("observers: clusterNameSelect updates clusterNumbers input", {
#     test_app <- make_test_app()
#
#     shiny::testServer(
#         app = test_app$app,
#         expr = {
#             # First create a group
#             session$setInputs(
#                 clusterName = "TestGroup",
#                 clusterNumbers = "5,6,7",
#                 applyName = 1
#             )
#             session$flushReact()
#
#             # Select the group via clusterNameSelect
#             session$setInputs(clusterNameSelect = "TestGroup")
#             session$flushReact()
#
#             # The clusterNumbers input should be updated (tested indirectly via observer)
#             # This triggers the observer that updates the text input
#         }
#     )
# })


# test_that("observers: groupsVar updates group1 and group2 choices", {
#     test_app <- make_test_app()
#     sce <- test_app$sce
#     metaD <- S4Vectors::metadata(sce)
#
#     # Find a factor column in experiment_info for testing
#     fact_cols <- unlist(lapply(metaD$experiment_info, is.factor), use.names = FALSE)
#     if (any(fact_cols)) {
#         factor_col <- names(fact_cols)[fact_cols][1]
#
#         shiny::testServer(
#             app = test_app$app,
#             expr = {
#                 # Set the groups variable
#                 session$setInputs(groupsVar = factor_col)
#                 session$flushReact()
#
#                 # The group1 and group2 inputs should be updated with factor levels
#                 # This is tested by the observer running without error
#             }
#         )
#     }
# })


# test_that("observers: group1 and group2 are mutually exclusive", {
#     test_app <- make_test_app()
#     sce <- test_app$sce
#     metaD <- S4Vectors::metadata(sce)
#
#     fact_cols <- unlist(lapply(metaD$experiment_info, is.factor), use.names = FALSE)
#     if (any(fact_cols)) {
#         factor_col <- names(fact_cols)[fact_cols][1]
#         levs <- levels(metaD$experiment_info[, factor_col])
#
#         shiny::testServer(
#             app = test_app$app,
#             expr = {
#                 # Set up groups
#                 session$setInputs(
#                     groupsVar = factor_col,
#                     group1 = levs[1],
#                     group2 = levs[2]
#                 )
#                 session$flushReact()
#
#                 # Select group1 = levs[1], group2 should not include levs[1]
#                 # This is tested by the observer running without error
#             }
#         )
#     }
# })


test_that("observers: sample2PlotDb updates dfPlot", {
  test_app <- make_test_app()
  sce <- test_app$sce

  sample_ids <- unique(as.character(sce$sample_id))

  shiny::testServer(
    app = test_app$app,
    expr = {
      # Change samples to plot
      suppressWarnings(session$setInputs(samples2plot = sample_ids[1]))
      suppressWarnings(session$flushReact())

      # dfPlot should be updated
      expect_true(inherits(dfPlot(), "data.frame"))
    }
  )
})


test_that("observers: selected group triggers selectedUpdate2", {
  test_app <- make_test_app()

  shiny::testServer(
    app = test_app$app,
    expr = {
      initial_count <- selectedUpdate2()

      # Select "selected" in colorbyGroups
      suppressWarnings(session$setInputs(colorbyGroups = "selected"))
      suppressWarnings(session$flushReact())

      # selectedUpdate2 should have incremented
      expect_true(selectedUpdate2() > initial_count)
    }
  )
})


# Note: Testing env$outputList sync requires access to the internal env
# which is not exposed in testServer. This is tested indirectly via rv$outputList.
