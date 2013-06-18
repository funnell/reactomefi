context("ReactomeFIService")


service <- new("ReactomeFIService")

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
