context("ReactomeFIService")


service <- new("ReactomeFIService")

test_that("queryFIs constructs a network", {
    genes <- c("TP53", "PTEN", "EGFR")
    fis <- data.frame(
        first.protein = c("EGFR", "PTEN"),
        second.protein = c("TP53", "TP53"))
    expect_that(queryFIs(service, genes), equals(fis))
})
