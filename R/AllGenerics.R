##
# ReactomeFIService generics
##

#' Reactome FI Network Version (2009 or 2012)
#' 
#' Retrieve the Reactome FI network version.
#' 
#' @return character
setGeneric("version", function(object) {
    standardGeneric("version")
})

#' RESTful API URL
#'
#' Retrieve the base URL for the ReactomeFI RESTful API.
#'
#' @return character
setGeneric("serviceURL", function(object) {
    standardGeneric("serviceURL")
})

#' Query Functional Interactions
#'
#' Query the RESTful API for FIs between genes in a provided vector.
#'
#' @param object ReactomeFIService object.
#' @param genes Character vector of gene names.
#' @return data.frame Each row represents a functional interaction and
#'  comprises two columns - one for each gene in the interaction.
setGeneric("queryFIs", function(object, genes) {
    standardGeneric("queryFIs")
})

#' Cluster Functional Interaction Network
#'
#' Query the RESTful API to cluster a FI network. The network nodes (genes)
#' and the network modules they belong to are returned.
#'
#' @param object ReactomeFIService object.
#' @param fis Data frame of functional interactions (see output of queryFIs).
#' @return data.frame Each row contains the gene name and the module id.
setGeneric("cluster", function(object, fis) {
    standardGeneric("cluster")
})

#' Annotate Gene Set
#'
#' Annotate a gene set with enriched pathways, or GO terms
#'
#' @param object ReactomeFIService object.
#' @param genes Character vector of gene names.
#' @return data.frame Each row represents an annotation of the provided type
#'  and includes related information such as the p-value and FDR generated
#'  from enrichment analysis.
setGeneric("annotateGeneSet",
           function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("annotateGeneSet")
})

#' Annotate FI Network Module Gene Set
#'
#' Annotate a gene set from a FI network module with enriched pathways, or GO
#' terms
#'
#' @param object ReactomeFIService object.
#' @param module.nodes Data frame with network nodes (genes) and their module.
#' @return data.frame Each row represents an annotation of the provided type
#'  and includes related information such as the p-value and FDR generated
#'  from enrichment analysis and the module the annotation belongs to.
setGeneric("annotateModules",
           function(object, module.nodes,
                    type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("annotateModules")
})
