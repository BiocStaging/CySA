test_that("package loads without errors", {
    expect_true(requireNamespace("CySA", quietly = TRUE))
})
