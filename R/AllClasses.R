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
