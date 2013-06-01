#' ReactomeFIService Class
#'
#' Represents an interface to the ReactomeFI RESTful API. Available versions
#' include 2009 and 2012
#'
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
#' ReactomeFIService constructor
#'
#' @param version character version of Reactome FI network (2009 or 2012)
#' @return ReactomeFIService
#' 
#' @rdname ReactomeFIService
ReactomeFIService <- function(version = "2009") {
    return(new("ReactomeFIService", version = version))
}
