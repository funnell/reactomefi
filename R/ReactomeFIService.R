setMethod("version", signature("ReactomeFIService"), function(object) {
    object@version
})

setMethod("url", signature("ReactomeFIService"), function(object) {
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
    url <- paste(url(object), "queryFIs", sep="")
    genes.str <- paste(genes, collapse = "\t")
    opts <- list(httpheader = c("Content-Type" = "text/plain;charset=UTF-8",
                                Accept = "application/xml"))
    xml <- postForm(url,  genes = genes.str, .opts = opts)
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
