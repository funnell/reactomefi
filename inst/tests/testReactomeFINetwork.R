context("ReactomeFINetwork")


network <- ReactomeFINetwork("2013")

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
    expect_equal(length(enriched.modules), 2)
    expect_equal(enriched.modules, c(0, 1))

    enriched.modules <- annotateModules(test.network, min.module.size=4)$module
    enriched.modules <- enriched.modules[!duplicated(enriched.modules)]
    expect_equal(length(enriched.modules), 1)
    expect_equal(enriched.modules, c(0))

    enriched.modules <- annotateModules(test.network, min.module.size=10)
    expect_equal(nrow(enriched.modules), 0)
})

test_that("networks can be plotted", {
    genes <- c("TP53", "PTEN", "EGFR", "ATM", "CLTCL1", "GRM8", "GRM7")
    test.network <- build(network, genes)

    gg <- plot(test.network)
    expect_equal(class(gg), c("gg", "ggplot"))
    
    base <- ggplot_build(gg)
    expect_equal(length(base$data), 3)
    expect_equal(nrow(base$data[[1]]), 7)
    expect_equal(nrow(base$data[[2]]), 7)
    expect_equal(nrow(base$data[[3]]), 10)
    expect_equal(unique(base$data[[1]]$group), 1)

    test.network <- cluster(test.network)
    gg <- plot(test.network)
    expect_equal(class(gg), c("gg", "ggplot"))
    
    base <- ggplot_build(gg)
    expect_equal(length(base$data), 3)
    expect_equal(nrow(base$data[[1]]), 7)
    expect_equal(nrow(base$data[[2]]), 7)
    expect_equal(nrow(base$data[[3]]), 10)
    expect_equal(unique(base$data[[1]]$group), c(1, 2))

    gg <- plot(test.network, colour.modules = FALSE)
    expect_equal(class(gg), c("gg", "ggplot"))
    
    base <- ggplot_build(gg)
    expect_equal(length(base$data), 3)
    expect_equal(nrow(base$data[[1]]), 7)
    expect_equal(nrow(base$data[[2]]), 7)
    expect_equal(nrow(base$data[[3]]), 10)
    expect_equal(unique(base$data[[1]]$group), 1)
})
