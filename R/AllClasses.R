#' ReactomeFIService Class
#'
#' Represents an interface to the ReactomeFI RESTful API. Available versions
#' include 2009 and 2012.
#'
#' @import RCurl
#' @import XML
#' @rdname ReactomeFIService
setClass("ReactomeFIService",
    representation(version = "character"),
    prototype(version = "2009"),
    validity = function(object) {
        version <- object@version
        if (version != "2009" && version != "2012")
            return("version must be either 2009(default) or 2012")
        TRUE
    }
)

#' ReactomeFIService
#'
#' ReactomeFIService constructor.
#'
#' @param version character version of Reactome FI network (2009 or 2012)
#' @return ReactomeFIService
#' 
#' @rdname ReactomeFIService
ReactomeFIService <- function(version = c("2009", "2012")) {
    version <- match.arg(version)
    return(new("ReactomeFIService", version = version))
}


#' ReactomeFINetwork Class
#'
#' Represents a functional interaction network generated using the Reactome
#'  database and a given list of genes.
#'
#' @rdname ReactomeFINetwork
setClass("ReactomeFINetwork",
    representation(service = "ReactomeFIService", fis = "data.frame",
                   modules = "data.frame"),
    prototype(service = ReactomeFIService(), fis = data.frame(),
              modules = data.frame()),
    validity = function(object) {
        if (class(object@service) != "ReactomeFIService") {
            return("service must be of class ReactomeFIService")
        }
        TRUE
    }
)

#' ReactomeFINetwork
#'
#' ReactomeFINetwork constructor
#'
#' @param genes Character vector of gene names.
#' @param version Version of ReactomeFI network (2009 or 2012).
#' @return ReactomeFINetwork ReactomeFINetwork S4 object with FI network
#'  generated from gene list.
#'
#' @export
#' @rdname ReactomeFINetwork
ReactomeFINetwork <- function(genes, version = c("2009", "2012")) {
    version <- match.arg(version)
    service <- ReactomeFIService(version)
    network <- new("ReactomeFINetwork", service = service)
    network <- build(network, genes)
    return(network)
}
