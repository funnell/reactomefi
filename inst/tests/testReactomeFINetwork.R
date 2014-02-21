context("ReactomeFINetwork")


network <- ReactomeFINetwork("2012")

test_that("network build does nothing if number of genes <= 1", {
    # queryFIs run on this gene return multiplt FIs
    genes <- c("CEP135")
    gene.fis <- fis(build(network, genes))
    expect_that(gene.fis, equals(data.frame()))
})

test_that("network cluster properly filters modules below size threshold", {
    genes <- c("TP53", "PTEN", "EGFR", "ATM", "CLTCL1", "GRM8", "GRM7")
    test.network <- build(network, genes)
    test.network <- cluster(test.network)

    enriched.modules <- annotateModules(test.network)$module
    enriched.modules <- enriched.modules[!duplicated(enriched.modules)]
    expect_that(length(enriched.modules), equals(2))
    expect_that(enriched.modules, equals(c(0, 1)))

    enriched.modules <- annotateModules(test.network, min.module.size=4)$module
    enriched.modules <- enriched.modules[!duplicated(enriched.modules)]
    expect_that(length(enriched.modules), equals(1))
    expect_that(enriched.modules, equals(c(0)))
})
