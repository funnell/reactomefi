setMethod("version", signature("ReactomeFIService"), function(object) {
    object@version
})

setMethod("serviceURL", signature("ReactomeFIService"), function(object) {
    version <- version(object)
    base.url = "http://reactomews.oicr.on.ca:8080/"
    if (version == "2012") {
        return(paste(base.url, "caBigR3WebApp2012/FIService/network/", sep=""))
    } 
    return(paste(base.url, "caBigR3WebApp/FIService/network/", sep=""))
})

setMethod("queryFIs",
          signature("ReactomeFIService", genes = "character"),
          function(object, genes) {
    service.url <- paste(serviceURL(object), "queryFIs", sep="")
    genes.str <- paste(genes, collapse = "\t")
    opts <- list(httpheader = c("Content-Type" = "text/plain;charset=UTF-8",
                                Accept = "application/xml"))
    xml <- postForm(service.url,  genes = genes.str, .opts = opts)
    doc <- xmlInternalTreeParse(xml)
    interactions <- xpathApply(doc, "//interaction", function(x) {
        info <- xmlChildren(x)
        first.protein <- xmlValue(xmlChildren(info$firstProtein)$name)
        second.protein <- xmlValue(xmlChildren(info$secondProtein)$name)
        data.frame(first.protein = first.protein,
                   second.protein = second.protein)
    })
    do.call(rbind, interactions)
})

fis2str <- function(fis) {
    fis[, "first.protein"] <- as.character(fis[, "first.protein"])
    fis[, "second.protein"] <- as.character(fis[, "second.protein"])
    fis <- cbind(data.frame(id=1:nrow(fis)), fis)

    fis.list = c()
    for (i in 1:nrow(fis)) {
        first.protein <- fis[i, "first.protein"]
        second.protein <- fis[i, "second.protein"]
        if (first.protein < second.protein) {
            fi.str <- paste(first.protein, second.protein, sep = "\t")
        } else {
            fi.str <- paste(second.protein, first.protein, sep = "\t")
        }
        fi.str <- paste(fis[i, "id"], fi.str, sep="\t")
        fis.list <- c(fis.list, fi.str)
    }
    fis.str <- paste(fis.list, collapse="\n")
    return(fis.str)
}

setMethod("cluster",
          signature("ReactomeFIService", "data.frame"),
          function(object, fis) {
    fis.str <- fis2str(fis)

    service.url <- paste(serviceURL(object), "cluster", sep="")
    opts <- list(httpheader = c("Content-Type" = "text/plain;charset=UTF-8",
                                Accept = "application/xml"))
    xml <- postForm(service.url,  queryFIs = fis.str, .opts = opts,
                    style="post", .contentEncodeFun=function(x) x)
    doc <- xmlInternalTreeParse(xml)
    clusters <- xpathApply(doc, "//geneClusterPairs", function(x) {
        info <- xmlChildren(x)
        cluster <- xmlValue(info$cluster)
        gene <- xmlValue(info$geneId)
        data.frame(cluster = cluster, gene = gene)
    })
    do.call(rbind, clusters)
})
