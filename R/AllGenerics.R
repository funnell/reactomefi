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

setGeneric("cluster", function(object, fis) {
    standardGeneric("cluster")
})
