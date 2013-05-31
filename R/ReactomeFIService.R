setMethod("version", signature("ReactomeFIService"), function(object) {
    return(object@version)
})

setMethod("serviceURL", signature("ReactomeFIService"), function(object) {
    version <- version(object)
    base.url = "http://reactomews.oicr.on.ca:8080/"
    if (version == "2012") {
        return(paste(base.url, "caBigR3WebApp2012/FIService/network/", sep=""))
    } 
    return(paste(base.url, "caBigR3WebApp/FIService/network/", sep=""))
})

getPostXML <- function(url, body) {
    text.gatherer <- basicTextGatherer()
    opts <- list(httpheader = c("Content-Type" = "text/plain;charset=UTF-8",
                                "Accept" = "application/xml"))
    curlPerform(postfields = body, url = url, .opts = opts,
                writefunction = text.gatherer$update)
    xml <- xmlInternalTreeParse(text.gatherer$value())
    return(xml)
}

setMethod("queryFIs",
          signature("ReactomeFIService", genes = "character"),
          function(object, genes) {
    service.url <- paste(serviceURL(object), "queryFIs", sep="")
    genes.str <- paste(genes, collapse = "\t")
    doc <- getPostXML(service.url, genes.str)
    interactions <- xpathApply(doc, "//interaction", function(x) {
        info <- xmlChildren(x)
        first.protein <- xmlValue(xmlChildren(info$firstProtein)$name)
        second.protein <- xmlValue(xmlChildren(info$secondProtein)$name)
        data.frame(first.protein = first.protein,
                   second.protein = second.protein)
    })
    return(do.call(rbind, interactions))
})

setMethod("annotateGeneSet",
          signature("ReactomeFIService", "character", "character"),
          function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    type <- match.arg(type)
    service.url <- paste(serviceURL(object), "annotateGeneSet/", type, sep="")
    genes.str <- paste(genes, collapse = "\n")
    doc <- getPostXML(service.url, genes.str)
    annotations <- xpathApply(doc, "//annotations", function(x) {
        info <- xmlChildren(x)
        data.frame(topic = xmlValue(info$topic),
                   hit.num = xmlValue(info$hitNumber),
                   number.in.topic = xmlValue(info$numberInTopic),
                   ratio.of.topic = xmlValue(info$ratioOfTopic),
                   p.value = xmlValue(info$PValue),
                   fdr = xmlValue(info$fdr),
                   hits = paste(xpathSApply(x, "./hitIds", xmlValue),
                                collapse = ","),
                   stringsAsFactors = FALSE)
    })
    annotations <- do.call(rbind, annotations)
    annotations$hit.num <- as.numeric(annotations$hit.num)
    annotations$number.in.topic <- as.numeric(annotations$number.in.topic)
    annotations$ratio.of.topic <- as.numeric(annotations$ratio.of.topic)
    annotations$p.value <- as.numeric(annotations$p.value)
    annotations$fdr <- gsub("<", "", annotations$fdr)
    annotations$fdr <- as.numeric(annotations$fdr)

    return(annotations[order(annotations$fdr), ])
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
    service.url <- paste(serviceURL(object), "cluster", sep="")
    fis.str <- fis2str(fis)
    doc <- getPostXML(service.url, fis.str)
    clusters <- xpathApply(doc, "//geneClusterPairs", function(x) {
        info <- xmlChildren(x)
        cluster <- xmlValue(info$cluster)
        gene <- xmlValue(info$geneId)
        data.frame(cluster = cluster, gene = gene)
    })
    return(do.call(rbind, clusters))
})
