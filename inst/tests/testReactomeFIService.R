context("ReactomeFIService")


service <- new("ReactomeFIService")

test_that("queryFIs returns empty data.frame when no FIs found", {
    genes <- c("FOO", "BAR")
    fis <- data.frame(
        first.protein = character(0),
        second.protein = character(0))
    expect_that(queryFIs(service, genes), equals(fis))
})

test_that("queryFIs constructs a network", {
    genes <- c("TP53", "PTEN", "EGFR")
    fis <- data.frame(
        first.protein = c("EGFR", "PTEN"),
        second.protein = c("TP53", "TP53"))
    expect_that(queryFIs(service, genes), equals(fis))
})

test_that("queryBuildNetwork constructs a network with linkers", {
    genes <- c("BECN1", "RAE1")
    fis <- data.frame(
        first.protein = c("XPO1", "XPO1"),
        second.protein = c("BECN1", "RAE1"))
    expect_that(queryBuildNetwork(service, genes), equals(fis))
})

test_that("queryCluster clusters a network", {
    fis <- data.frame(
        first.protein = c("A", "B", "A", "A", "D", "E", "D"),
        second.protein = c("B", "C", "C", "D", "E", "F", "F"))
    modules <- data.frame(
        gene = c("D", "E", "F", "A", "B", "C"),
        module = c(0, 0, 0, 1, 1, 1),
        stringsAsFactors = FALSE)
    expect_that(queryCluster(service, fis), equals(modules))
})

test_that("queryFIsBetween returns FIs between certain genes", {
    gene.pairs <- data.frame(
        first.protein = c("EGFR", "PTEN", "EGFR"),
        second.protein = c("TP53", "TP53", "PTEN"))
    fis <- data.frame(
        first.protein = c("EGFR", "PTEN"),
        second.protein = c("TP53", "TP53"))
    expect_that(queryFIsBetween(service, gene.pairs), equals(fis))
})

test_that("queryEdge returns detailed information for an edge", {
    info.subset <- data.frame(
        accession = c("P00533", "P04637"),
        db.name = c("UniProt", "UniProt"),
        short.name = c("EGFR", "TP53"))
    subset.cols <- c("accession", "db.name", "short.name")
    expect_that(queryEdge(service, "TP53", "EGFR")[subset.cols],
                equals(info.subset))
})

test_that("queryAnnotateGeneSet returns annotations", {
    genes <- c("TP53", "PTEN", "EGFR")
    annotations <- queryAnnotateGeneSet(service, genes, "Pathway")
    expect_that(class(annotations), equals("data.frame"))
    expect_that(dim(annotations), equals(c(125, 7)))
    expect_that(typeof(annotations$topic), equals("character"))
    expect_that(typeof(annotations$hit.num), equals("double"))
    expect_that(typeof(annotations$number.in.topic), equals("double"))
    expect_that(typeof(annotations$ratio.of.topic), equals("double"))
    expect_that(typeof(annotations$p.value), equals("double"))
    expect_that(typeof(annotations$fdr), equals("double"))
    expect_that(typeof(annotations$hits), equals("character"))
})

test_that("queryAnnotateModules returns annotations", {
    modules <- data.frame(
        gene = c("NDEL1", "PLK1", "CENPC1", "NUF2", "CENPE", "CSNK1D", "ODF2",
                 "CEP164", "CDK5RAP2"),
        module = c(0, 0, 0, 0, 0, 1, 1, 1, 1), stringsAsFactors=F)
    annotations <- queryAnnotateModules(service, modules, "Pathway")
    expect_that(class(annotations), equals("data.frame"))
    expect_that(dim(annotations), equals(c(21, 8)))
    expect_that(typeof(annotations$module), equals("double"))
    expect_that(unique(annotations$module), equals(c(0,1)))
    expect_that(typeof(annotations$topic), equals("character"))
    expect_that(typeof(annotations$hit.num), equals("double"))
    expect_that(typeof(annotations$number.in.topic), equals("double"))
    expect_that(typeof(annotations$ratio.of.topic), equals("double"))
    expect_that(typeof(annotations$p.value), equals("double"))
    expect_that(typeof(annotations$fdr), equals("double"))
    expect_that(typeof(annotations$hits), equals("character"))
})

test_that("queryHotNetAnalysis returns modules", {
    mod.genes <- list(c("ING4", "MSH6", "NBN", "MSH2", "AIFM1", "TP53", "MLH1",
                        "CHEK1", "TRIM24", "TRRAP", "ATM", "TOP1", "BAX",
                        "PMS2", "MDM4", "EP400"),
                      c("PIK3C2B", "PDGFRA", "PDGFRB", "PIK3CA", "KIT", "ZEB1",
                        "PTEN", "PIK3R1"),
                      c("KLF6", "BCL11A","KLF4"), c("PRKD2", "NF1", "DST"),
                      c("PAX5", "MYCN"), c("PTCH1", "SHH"),
                      c("NOTCH1", "DTX3"), c("EPHA7", "EPHA2"))
    data(GlioblastomaMAF)
    gene.scores <- maf2genescores(glioblastoma.maf)
    hotnet <- HotNet(gene.scores, "2012")
    res.mod.genes <- lapply(modules(hotnet), function(x) x$genes)
    
    expect_that(res.mod.genes, equals(mod.genes))

    mask <- sapply(modules(hotnet), function(x) x$fdr <= 0.25)
    expect_that(res.mod.genes[mask], equals(mod.genes[1:2]))
})
