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
#' @export
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
#' @param version Version of ReactomeFI network (2009 or 2012).
#' @return ReactomeFINetwork ReactomeFINetwork S4 object with FI network
#'  generated from gene list.
#'
#' @export
#' @rdname ReactomeFINetwork
ReactomeFINetwork <- function(version = c("2009", "2012")) {
    version <- match.arg(version)
    service <- ReactomeFIService(version)
    network <- new("ReactomeFINetwork", service = service)
    return(network)
}


#' HotNet Class
#'
#' Represents the results of a HotNet analysis
#'
#' @rdname HotNet
setClass("HotNet",
    representation(service = "ReactomeFIService", gene.scores = "data.frame",
                   delta = "numeric", fdr = "numeric",
                   permutations = "numeric", auto.delta = "logical",
                   modules = "list"),
    prototype(service = ReactomeFIService(), gene.scores = data.frame(),
              delta = 1e-4, fdr = 0.25, permutations = 100, auto.delta = F,
              modules = data.frame()),
    validity = function(object) {
        if (object@delta < 0) {
            return("delta must be positive")
        }
        if (object@fdr < 0 || object@fdr > 1) {
            return("fdr threshold must be >= 0 and <= 1")
        }
        if ((ceiling(object@permutations) != floor(object@permutations)) ||
            (object@permutations < 1) || (object@permutations > 1000)) {
            return("permutations must be a positive integer less than 1000")
        }
        TRUE
    }
)

#' HotNet
#'
#' HotNet constructor
#'
#' @param gene.scores data.frame where the first column contains gene names and
#'  the second column contains scores representing the proportion of samples
#'  in which the gene is mutated.
#' @param version Version of ReactomeFI network (2009 or 2012).
#' @param delta HotNet delta value
#' @param fdr False discovery rate threshold
#' @param permutations Number of permutations. Largest value is 1000.
#' @param auto.delta If true, algorithm will select a delta value. This option
#'  will make the analysis take more time to finish.
#' @return HotNet HotNet S4 object containing the results of
#'  HotNet analysis. These results should be filtered before using them to
#'  generate a network
#'
#' @export
#' @rdname HotNet
HotNet <- function(gene.scores, version = c("2009", "2012"), delta = 1e-4,
                   fdr = 0.25, permutations = 100, auto.delta = F) {
    version <- match.arg(version)
    service <- ReactomeFIService(version)
    res <- queryHotNetAnalysis(service, gene.scores, delta, fdr, permutations,
                               auto.delta)
    hotnet <- new("HotNet", service = service, gene.scores = gene.scores,
                  delta = res$delta, fdr = res$fdr,
                  permutations = res$permutations, auto.delta = res$auto.delta,
                  modules = res$modules)
    return(hotnet)
}
