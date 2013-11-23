context("ReactomeFIService")


network <- ReactomeFINetwork("2012")

test_that("network build does nothing if number of genes <= 1", {
    # queryFIs run on this gene return multiplt FIs
    genes <- c("CEP135")
    gene.fis <- fis(build(network, genes))
    expect_that(gene.fis, equals(data.frame()))
})
