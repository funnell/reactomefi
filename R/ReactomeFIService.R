#' @rdname version-methods
#' @aliases version,ReactomeFIService-method
setMethod("version", signature("ReactomeFIService"), function(object) {
    return(object@version)
})

#' @rdname serviceURL-methods
#' @aliases serviceURL,ReactomeFIService-method
setMethod("serviceURL", signature("ReactomeFIService"), function(object) {
    version <- version(object)
    base.url = "http://reactomews.oicr.on.ca:8080/"
    if (version == "2012") {
        return(paste(base.url, "caBigR3WebApp2012/FIService/network/",
                     sep = ""))
    } 
    return(paste(base.url, "caBigR3WebApp/FIService/network/", sep = ""))
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
    curlPerform(postfields = body, url = url, .opts = opts, .encoding="UTF-8",
                writefunction = text.gatherer$update)
    xml <- xmlInternalTreeParse(text.gatherer$value())
    return(xml)
}

#' Extract FIs From XML Document
#'
#' Extract FI gene pairs from an XML document returned from a server query.
#'
#' @param doc XML document
#' @return data.frame Data frame where each row represents an FI and each of
#'  the two columns in the data frame contains a gene involved in the FI.
extractFIs <- function(doc) {
    interactions <- xpathApply(doc, "//interaction", function(x) {
        info <- xmlChildren(x)
        first.protein <- xmlValue(xmlChildren(info$firstProtein)$name)
        second.protein <- xmlValue(xmlChildren(info$secondProtein)$name)
        data.frame(first.protein = first.protein,
                   second.protein = second.protein, stringsAsFactors = FALSE)
    })
    fis <- do.call(rbind, interactions)
    if (is.null(fis)) return(data.frame())
    return(fis)
}

#' @rdname queryFIs-methods
#' @aliases queryFIs,ReactomeFIService,character-method
setMethod("queryFIs",
          signature("ReactomeFIService", "character"),
          function(object, genes) {
    service.url <- paste0(serviceURL(object), "queryFIs")
    genes.str <- paste(genes, collapse = "\t")
    doc <- getPostXML(service.url, genes.str)
    return(extractFIs(doc))
})

#' @rdname queryBuildNetwork-methods
#' @aliases queryBuildNetwork,ReactomeFIService,character-method
setMethod("queryBuildNetwork",
          signature("ReactomeFIService", "character"),
          function(object, genes) {
    service.url <- paste0(serviceURL(object), "buildNetwork")
    genes.str <- paste(genes, collapse = "\t")
    doc <- getPostXML(service.url, genes.str)
    return(extractFIs(doc))
})

#' @rdname queryFIsBetween-methods
#' @aliases queryFIsBetween,ReactomeFIService,data.frame-method
setMethod("queryFIsBetween",
          signature("ReactomeFIService", "data.frame"),
          function(object, gene.pairs) {
    service.url <- paste0(serviceURL(object), "queryFIsBetween")
    gene.pairs <- as.matrix(gene.pairs)
    first.str <- paste(gene.pairs[, 1], collapse = ",")
    second.str <- paste(gene.pairs[, 2], collapse = ",")
    pairs.str <- paste(first.str, second.str, sep = "\n")
    doc <- getPostXML(service.url, pairs.str)
    return(extractFIs(doc))
})

#' extract Protein Info
#'
#' Extract protein information including accession ID and DB name, protein
#'  name, and sequence.
#'
#' @param protein.node XML node containing protein information
#' @return data.frame Data frame where each row corresponds to a protein and
#'  the columns contain the information mentioned above.
extractProteinInfo <- function(protein.node) {
    accession <- xmlValue(xmlChildren(protein.node)$accession)
    db.name <- xmlValue(xmlChildren(protein.node)$dbName)
    name <- xmlValue(xmlChildren(protein.node)$name)
    short.name <- xmlValue(xmlChildren(protein.node)$shortName)
    prot.seq <- xmlValue(xmlChildren(protein.node)$sequence)
    info <- data.frame(accession = accession,
                       db.name = db.name,
                       short.name = short.name,
                       name = name,
                       sequence = prot.seq,
                       stringsAsFactors = FALSE)
    return(info)
}

#' @rdname queryEdge-methods
#' @aliases queryEdge,ReactomeFIService,character,character-method
setMethod("queryEdge",
          signature("ReactomeFIService", "character", "character"),
          function(object, name1, name2) {
    service.url <- paste0(serviceURL(object), "queryEdge")
    edge.str <- paste(name1, name2, sep = "\t")
    doc <- getPostXML(service.url, edge.str)
    first.prot <- getNodeSet(doc, "//firstProtein", fun = extractProteinInfo)
    second.prot <- getNodeSet(doc, "//secondProtein", fun = extractProteinInfo)
    return(do.call(rbind, c(first.prot, second.prot)))
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
    fis <- cbind(data.frame(id = 1:nrow(fis)), fis)

    fis.list = c()
    for (i in 1:nrow(fis)) {
        first.protein <- fis[i, "first.protein"]
        second.protein <- fis[i, "second.protein"]
        if (first.protein < second.protein) {
            fi.str <- paste(first.protein, second.protein, sep = "\t")
        } else {
            fi.str <- paste(second.protein, first.protein, sep = "\t")
        }
        fi.str <- paste(fis[i, "id"], fi.str, sep = "\t")
        fis.list <- c(fis.list, fi.str)
    }
    fis.str <- paste(fis.list, collapse = "\n")
    return(fis.str)
}

#' @rdname queryCluster-methods
#' @aliases queryCluster,ReactomeFIService,data.frame-method
setMethod("queryCluster",
          signature("ReactomeFIService", "data.frame"),
          function(object, fis) {
    service.url <- paste0(serviceURL(object), "cluster")
    fis.str <- fis2str(fis)
    doc <- getPostXML(service.url, fis.str)
    modules <- xpathApply(doc, "//geneClusterPairs", function(x) {
        info <- xmlChildren(x)
        module <- xmlValue(info$cluster)
        gene <- xmlValue(info$geneId)
        data.frame(gene = gene, module = module, stringsAsFactors = FALSE)
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
    annotations <- do.call(rbind, annotations)

    if (is.null(annotations)) return(data.frame())

    annotations$hit.num <- as.numeric(annotations$hit.num)
    annotations$number.in.topic <- as.numeric(annotations$number.in.topic)
    annotations$ratio.of.topic <- as.numeric(annotations$ratio.of.topic)
    annotations$p.value <- as.numeric(annotations$p.value)
    annotations$fdr <- gsub("<", "", annotations$fdr)
    annotations$fdr <- as.numeric(annotations$fdr)
    return(annotations[order(annotations$fdr), ])
}

#' @rdname queryAnnotateGeneSet-methods
#' @aliases queryAnnotateGeneSet,ReactomeFIService,character,character-method
setMethod("queryAnnotateGeneSet",
          signature("ReactomeFIService", "character", "character"),
          function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    type <- match.arg(type)
    service.url <- paste0(serviceURL(object), "annotateGeneSet/", type)
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
    for (i in 1:ncol(dat)) {
        if (class(dat[, i]) == "numeric") {
            dat[, i] <- format(dat[, i], trim = TRUE)
        }
    }
    tsv <- apply(dat, 1, function(x) paste(x, collapse = "\t"))
    tsv <- paste(tsv, collapse = "\n")
    return(tsv)
}

#' @rdname queryAnnotateModules-methods
#' @aliases queryAnnotateModules,ReactomeFIService,data.frame,character-method
setMethod("queryAnnotateModules",
          signature("ReactomeFIService", "data.frame", "character"),
          function(object, module.nodes,
                   type = c("Pathway", "BP", "CC", "MF")) {
    type <- match.arg(type)
    service.url <- paste0(serviceURL(object), "annotateModules/", type)
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
    module.annotations$module <- module.annotations$module - 1
    return(module.annotations)
})

#' @rdname queryHotNetAnalysis-methods
#' @aliases queryHotNetAnalysis,ReactomeFIService-method
setMethod("queryHotNetAnalysis",
          signature("ReactomeFIService"),
          function(object, gene.scores, delta, fdr, permutations, auto.delta) {
    service.url <- paste0(serviceURL(object), "hotnetAnalysis")
    gene.score.str <- df2tsv(gene.scores)
    query <- paste0(gene.score.str, "\n")
    query <- paste0(query, "fdrCutoff:", fdr, "\n")
    query <- paste0(query, "permutationNumber:", permutations, "\n")
    if (!auto.delta) {
        delta <- format(delta, scientific = F)
        query <- paste0(query, "delta:", delta, "\n")
    }
    doc <- getPostXML(service.url, query)

    result <- xmlChildren(doc)$hotNetResult
    delta <- as.numeric(xmlValue(xmlChildren(result)$delta))
    fdr <- as.numeric(xmlValue(xmlChildren(result)$fdrThreshold))
    permutations <- as.numeric(xmlValue(xmlChildren(result)$permutation))
    auto.delta <- as.logical(xmlValue(xmlChildren(result)$useAutoDelta))
    modules <- xpathApply(doc, "//modules", function(x) {
        fdr <- as.numeric(xmlValue(xmlChildren(x)$fdr))
        p.value <- as.numeric(xmlValue(xmlChildren(x)$pvalue))
        genes <- unlist(xpathApply(x, "genes", xmlValue))
        list(fdr = fdr, p.value = p.value, genes = genes)
    })
    return(list(delta = delta, fdr = fdr, permutations = permutations,
                auto.delta = auto.delta, modules = modules))
})
