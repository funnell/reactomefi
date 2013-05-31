##
# ReactomeFIService generics
##

setGeneric("version", function(object) {
    standardGeneric("version")
})

setGeneric("serviceURL", function(object) {
    standardGeneric("serviceURL")
})

setGeneric("queryFIs", function(object, genes) {
    standardGeneric("queryFIs")
})

setGeneric("annotateGeneSet",
           function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("annotateGeneSet")
})

setGeneric("cluster", function(object, fis) {
    standardGeneric("cluster")
})

setGeneric("annotateModules",
           function(object, module.nodes,
                    type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("annotateModules")
})
