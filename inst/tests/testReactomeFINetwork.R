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

    anns <- annotateModules(test.network, min.module.size = 4)
    enriched.modules <- anns$module
    enriched.modules <- enriched.modules[!duplicated(enriched.modules)]
    expect_equal(length(enriched.modules), 1)
    expect_equal(enriched.modules, c(0))

    expect_warning(ann <- annotateModules(test.network, min.module.size = 10))
    expect_equal(dim(ann), c(0, 0))
})

test_that("network can be build with linkers", {
    test.fis <- data.frame(
        first.protein = c("XPO1", "XPO1"),
        second.protein = c("BECN1", "RAE1"),
        stringsAsFactors = FALSE)

    genes <- c("BECN1", "RAE1")
    test.network <- build(network, genes, use.linkers = TRUE)

    expect_equal(fis(test.network), test.fis)
})

test_that("network annotate can ignore linkers", {
    genes <- c("BECN1", "RAE1")
    test.network <- build(network, genes, use.linkers = TRUE)

    annotations <- annotate(test.network, "Pathway", include.linkers = TRUE)
    expect_equal(class(annotations), "data.frame")
    expect_equal(ncol(annotations), 7)
    expect_true(nrow(annotations) > 0)
    expect_true("XPO1" %in% annotations$hits)

    annotations <- annotate(test.network, "Pathway")
    expect_equal(class(annotations), "data.frame")
    expect_equal(ncol(annotations), 7)
    expect_true(nrow(annotations) > 0)
    expect_false("XPO1" %in% annotations$hits)
})

test_that("clustered network annotate can ignore linkers", {
    genes <- c("BECN1", "RAE1")
    test.network <- build(network, genes, use.linkers = TRUE)
    test.network <- cluster(test.network)

    annotations <- annotateModules(test.network, "Pathway",
                                   include.linkers = TRUE)
    expect_equal(class(annotations), "data.frame")
    expect_equal(ncol(annotations), 8)
    expect_true(nrow(annotations) > 0)
    expect_true("XPO1" %in% annotations$hits)

    annotations <- annotateModules(test.network, "Pathway")
    expect_equal(class(annotations), "data.frame")
    expect_equal(ncol(annotations), 8)
    expect_true(nrow(annotations) > 0)
    expect_false("XPO1" %in% annotations$hits)
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

test_that("networks with linkers can be plotted", {
    genes <- c("BECN1", "RAE1")
    test.network <- build(network, genes, use.linkers = TRUE)

    gg <- plot(test.network)
    expect_equal(class(gg), c("gg", "ggplot"))
    base <- ggplot_build(gg)
    expect_true("shape" %in% colnames(base$data[[1]]))

    gg <- plot(test.network, indicate.linkers = FALSE)
    expect_equal(class(gg), c("gg", "ggplot"))
    base <- ggplot_build(gg)
    expect_false("shape" %in% colnames(base$data[[1]]))
})
