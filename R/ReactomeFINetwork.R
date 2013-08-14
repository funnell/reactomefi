setMethod("version", signature("ReactomeFINetwork"), function(object) {
    return(version(object@service))
})

setMethod("service", signature("ReactomeFINetwork"), function(object) {
    return(object@service)
})

setMethod("fis", signature("ReactomeFINetwork"), function(object) {
    object@fis
})

setMethod("fis<-", signature("ReactomeFINetwork", "data.frame"),
          function(object, value) {
    object@fis <- value
    object
})

setMethod("modules", signature("ReactomeFINetwork"), function(object) {
    object@modules
})

setMethod("modules<-", signature("ReactomeFINetwork", "data.frame"),
          function(object, value) {
    object@modules <- value
    object
})

#' Build FI Network
#'
#' Build FI network from a list of genes.
#'
#' @param object ReactomeFINetwork object.
#' @param genes Character vector of gene names
#' @return ReactomeFINetwork ReactomeFINetwork object with fis attribute set
setMethod("build", signature("ReactomeFINetwork"),
          function(object, genes) {
    service <- service(object)
    fis(object) <- queryFIs(service, genes)
    object
})

setMethod("cluster", signature("ReactomeFINetwork"), function(object) {
    if (nrow(fis(object)) == 0) {
        warning("No FI network data found. Please build the network first.")
        return(object)
    }
    service <- service(object)
    modules(object) <- queryCluster(service, fis(object))
    object
})

setMethod("annotate", signature("ReactomeFINetwork", "character"),
          function(object, type = c("Pathway", "BP", "CC", "MF")) {
    if (nrow(fis(object)) == 0) {
        warning("No FI network data found. Please build the network first.")
        return(object)
    }
    type <- match.arg(type)
    service <- service(object)
    fis <- fis(object)
    genes <- union(fis$first.protein, fis$second.protein)
    return(queryAnnotateGeneSet(service, genes, type))
})

setMethod("annotateModules", signature("ReactomeFINetwork"),
          function(object, type = c("Pathway", "BP", "CC", "MF")) {
    if (nrow(modules(object)) == 0) {
        message <- paste("No FI network module data found. Please cluster the",
                         "network first.")
        warning(message)
        return(object)
    }
    type <- match.arg(type)
    service <- service(object)
    return(queryAnnotateModules(service, modules(object), type))
})
