##
# ReactomeFIService generics
##

#' Reactome FI Network Version (2009 or 2012)
#' 
#' Retrieve the Reactome FI network version
#' 
#' @return character
setGeneric("version", function(object) {
    standardGeneric("version")
})

#' RESTful API URL
#'
#' Retrieve the base URL for the ReactomeFI RESTful API
#'
#' @return character
setGeneric("serviceURL", function(object) {
    standardGeneric("serviceURL")
})

#' Query Functional Interactions
#'
#' Query the RESTful API for FIs between genes in a provided vector
#'
#' @param object ReactomeFIService object
#' @param genes character vector of gene names
#' @return data.frame
setGeneric("queryFIs", function(object, genes) {
    standardGeneric("queryFIs")
})

#' Cluster Functional Interaction Network
#'
#' Query the RESTful API to cluster a FI network. The network nodes (genes)
#' and the network modules they belong to are returned.
#'
#' @param object ReactomeFIService object
#' @param fis data frame of functional interactions (see output of queryFIs)
#' @return data.frame
setGeneric("cluster", function(object, fis) {
    standardGeneric("cluster")
})

#' Annotate Gene Set
#'
#' Annotate a gene set with enriched pathways, or GO terms
#'
#' @param object ReactomeFIService object
#' @param genes character vector of gene names
#' @return data.frame
setGeneric("annotateGeneSet",
           function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("annotateGeneSet")
})

#' Annotate FI Network Module Gene Set
#'
#' Annotate a gene set from a FI network module with enriched pathways, or GO
#' terms
#'
#' @param object ReactomeFIService object
#' @param module.nodes data frame with network nodes (genes) and their module
#' @return data.frame
setGeneric("annotateModules",
           function(object, module.nodes,
                    type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("annotateModules")
})
