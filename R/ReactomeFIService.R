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

#' Get POST Query XML Response
#'
#' Make a POST request to a URL with a query string to be provided to
#' curlPerform's postfields argument. curlPerform is used instead of postForm
#' to allow sending unkeyed data.
#'
#' @param url Character string containing URL to post query to.
#' @param body Character string containing data to send in the post request.
#' @return XMLDocument XML document generated from the POST query response.
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

#' FIs to String
#'
#' Convert a FI data frame into a string according to conventions used in the
#' ReactomeFI API.
#'
#' @param fis data frame of FIs (each row contains two gene names)
#' @return character Character string in TSV format where rows are separated
#'  by "\\n" and columns within rows are separated by "\\t"
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

setMethod("queryCluster",
          signature("ReactomeFIService", "data.frame"),
          function(object, fis) {
    service.url <- paste(serviceURL(object), "cluster", sep="")
    fis.str <- fis2str(fis)
    doc <- getPostXML(service.url, fis.str)
    modules <- xpathApply(doc, "//geneClusterPairs", function(x) {
        info <- xmlChildren(x)
        module <- xmlValue(info$cluster)
        gene <- xmlValue(info$geneId)
        data.frame(gene = gene, module = module, stringsAsFactors = F)
    })
    modules <- do.call(rbind, modules)
    modules$module <- as.numeric(modules$module)
    return(modules)
})

extractAnnotations <- function(xml.node) {
    annotations <- xpathApply(xml.node, "./annotations", function(x) {
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
    if (length(annotations) == 0) return(NA)
    annotations <- do.call(rbind, annotations)
    annotations$hit.num <- as.numeric(annotations$hit.num)
    annotations$number.in.topic <- as.numeric(annotations$number.in.topic)
    annotations$ratio.of.topic <- as.numeric(annotations$ratio.of.topic)
    annotations$p.value <- as.numeric(annotations$p.value)
    annotations$fdr <- gsub("<", "", annotations$fdr)
    annotations$fdr <- as.numeric(annotations$fdr)
    return(annotations[order(annotations$fdr), ])
}

setMethod("queryAnnotateGeneSet",
          signature("ReactomeFIService", "character", "character"),
          function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    type <- match.arg(type)
    service.url <- paste(serviceURL(object), "annotateGeneSet/", type, sep="")
    genes.str <- paste(genes, collapse = "\n")
    doc <- getPostXML(service.url, genes.str)
    annot.node <- xmlChildren(doc)$moduleGeneSetAnnotations
    annot.node <- xmlChildren(annot.node)$moduleGeneSetAnnotation
    annotations <- extractAnnotations(annot.node)
    return(annotations)
})

#' Data Frame to TSV
#'
#' Convert a data frame (not including headers) into a TSV string.
#'
#' @param dat Data frame to be converted to a TSV string.
#' @return character Each row is separated by "\\n" and each column within a
#'  row is separated by a "\\t".
df2tsv <- function(dat) {
    tsv <- apply(dat, 1, function(x) paste(x, collapse = "\t"))
    tsv <- paste(tsv, collapse = "\n")
    return(tsv)
}

setMethod("queryAnnotateModules",
          signature("ReactomeFIService", "data.frame", "character"),
          function(object, module.nodes,
                   type = c("Pathway", "BP", "CC", "MF")) {
    type <- match.arg(type)
    service.url <- paste(serviceURL(object), "annotateModules/", type, sep="")
    query <- df2tsv(module.nodes)
    doc <- getPostXML(service.url, query)
    module.annotations <- xpathApply(doc, "//moduleGeneSetAnnotation",
                                     function(x) {
        module <- xmlValue(xmlChildren(x)$module)
        annotations <- extractAnnotations(x)
        if (all(is.na(annotations))) return(annotations)
        cbind(data.frame(module = module), annotations)
    })
    module.annotations <- module.annotations[!is.na(module.annotations)]
    module.annotations <- do.call(rbind, module.annotations)
    module.annotations$module <- as.numeric(module.annotations$module)
    return(module.annotations)
})
