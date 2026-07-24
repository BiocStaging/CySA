test_that("server initialises reactive outputList with Rest group", {
    test_app <- make_test_app()

    shiny::testServer(
        app = test_app$app,
        expr = {
            expect_true("Rest" %in% names(rv$outputList))
            expect_length(rv$outputList$Rest, 10)
        }
    )
})

test_that("server updates rv$outputList when group is named", {
    test_app <- make_test_app()

    shiny::testServer(
        app = test_app$app,
        expr = {
            # Initial state: only Rest.
            expect_equal(sort(names(rv$outputList)), "Rest")

            # Simulate typing a name and applying a cluster selection.
            session$setInputs(
                clusterName = "MyGroup",
                clusterNumbers = "1,2",
                applyName = 1
            )
            session$flushReact()

            expect_true("MyGroup" %in% names(rv$outputList))
            expect_true(all(rv$outputList$MyGroup == c(1L, 2L)))
        }
    )
})

test_that("server reactive values are isolated from external env changes", {
    test_app <- make_test_app()

    shiny::testServer(
        app = test_app$app,
        expr = {
            expect_s4_class(test_app$sce, "SingleCellExperiment")
            expect_true("Rest" %in% names(rv$outputList))
            expect_type(dListRV(), "list")
        }
    )
})
