test_that("package loads without errors", {
    expect_true(requireNamespace("CySA", quietly = TRUE))
})

test_that("source code avoids seq(n) in favour of seq_len(n)", {
    r_files <- list.files(system.file("..", "R", package = "CySA"), pattern = "\\.R$", full.names = TRUE)
    for (f in r_files) {
        lines <- readLines(f)
        # Match seq( followed by a bare numeric symbol, not seq_len/along/int/from/length.out
        bad <- grep("\\bseq\\s*\\(\\s*[A-Za-z_][A-Za-z0-9_]*\\s*\\)", lines, value = TRUE)
        bad <- grep("seq_(len|along|int)|\\bseq\\s*\\(\\s*from\\s*=|\\bseq\\s*\\(\\s*length\\.out\\s*=", bad, value = TRUE, invert = TRUE)
        expect(
            length(bad) == 0,
            failure_message = paste("Found seq(n) style call in", f, ":", paste(bad, collapse = "; "))
        )
    }
})
