#' @include AllClasses.R
#' @include AllGenerics.R
#' @include ReactomeFIService.R
NULL

#' @rdname version-methods
#' @aliases version,ReactomeFINetwork-method
setMethod("version", signature("ReactomeFINetwork"), function(object) {
    return(version(object@service))
})

#' @rdname service-methods
#' @aliases service,ReactomeFINetwork-method
setMethod("service", signature("ReactomeFINetwork"), function(object) {
    return(object@service)
})

#' @rdname genes-methods
#' @aliases genes,ReactomeFINetwork-method
setMethod("genes", signature("ReactomeFINetwork"), function(object) {
    object@genes
})

#' @rdname genes-methods
#' @aliases genes<-,ReactomeFINetwork,character-method
setMethod("genes<-", signature("ReactomeFINetwork", "character"),
          function(object, value) {
    object@genes <- value
    object
})

#' @rdname fis-methods
#' @aliases fis,ReactomeFINetwork-method
setMethod("fis", signature("ReactomeFINetwork"), function(object) {
    object@fis
})

#' @rdname fis-methods
#' @aliases fis<-,ReactomeFINetwork,data.frame-method
setMethod("fis<-", signature("ReactomeFINetwork", "data.frame"),
          function(object, value) {
    object@fis <- value
    object
})

#' @rdname modules-methods
#' @aliases modules,ReactomeFINetwork-method
setMethod("modules", signature("ReactomeFINetwork"), function(object) {
    object@modules
})

#' @rdname modules-methods
#' @aliases modules<-,ReactomeFINetwork,data.frame-method
setMethod("modules<-", signature("ReactomeFINetwork", "data.frame"),
          function(object, value) {
    object@modules <- value
    object
})

#' @rdname build-methods
#' @aliases build,ReactomeFINetwork,character-method
setMethod("build", signature("ReactomeFINetwork", "character"),
          function(object, genes, use.linkers = FALSE) {
    if (length(genes) > 1) {
        genes(object) <- genes
        service <- service(object)

        if (use.linkers) {
            fis(object) <- queryBuildNetwork(service, genes)
        } else {
            fis(object) <- queryFIs(service, genes)
        }
    }
    object
})

#' @rdname cluster-methods
#' @aliases cluster,ReactomeFINetwork-method
setMethod("cluster", signature("ReactomeFINetwork"), function(object) {
    if (nrow(fis(object)) == 0) {
        warning("No FI network data found. Please build the network first.")
        return(object)
    }
    service <- service(object)
    modules(object) <- queryCluster(service, fis(object))
    object
})

#' @rdname annotate-methods
#' @aliases annotate,ReactomeFINetwork,character-method
setMethod("annotate", signature("ReactomeFINetwork", "character"),
          function(object, type = c("Pathway", "BP", "CC", "MF"),
                   include.linkers = FALSE) {
    if (nrow(fis(object)) == 0) {
        warning("No FI network data found. Please build the network first.")
        return(object)
    }

    fis <- fis(object)
    fi.genes <- union(fis$first.protein, fis$second.protein)
    
    if (!include.linkers) {
        fi.genes <- fi.genes[fi.genes %in% genes(object)]
    }

    type <- match.arg(type)
    service <- service(object)
    return(queryAnnotateGeneSet(service, fi.genes, type))
})

#' @rdname annotateModules-methods
#' @aliases annotateModules,ReactomeFINetwork-method
setMethod("annotateModules", signature("ReactomeFINetwork"),
          function(object, type = c("Pathway", "BP", "CC", "MF"),
                   min.module.size = 1, include.linkers = FALSE) {
    if (nrow(modules(object)) == 0) {
        message <- paste("No FI network module data found. Please cluster the",
                         "network first.")
        warning(message)
        return(data.frame())
    }

    module.size.filter <- function(x) { if (nrow(x) >= min.module.size) x }
    network.modules <- ddply(modules(object), .(module), module.size.filter)
    if (nrow(network.modules) == 0) {
        warning("No modules left to annotate.")
        return(data.frame())
    }

    if (!include.linkers) {
        network.modules <- subset(network.modules, gene %in% genes(object))
    }

    type <- match.arg(type)
    service <- service(object)
    return(queryAnnotateModules(service, network.modules, type))
})

#' Layout ReactomeFINetwork
#
#' Retrieve the coordinates for the network vertices according to the
#'  specified layout algorithm
#
#' @param adj.mat Adjacency matrix
#' @param layout Layout algorithm as defined by sna's gplot.layout
#' @return data.frame DataFrame containing x and y coordinates of network
#'  vertices
layout.net <- function(adj.mat, layout.type) {
    layout.fname <- paste0("gplot.layout.", layout.type)
    if (exists(layout.fname)) {
        layout.f <- get(layout.fname)
        vertices <- layout.f(adj.mat, NULL)
        vertices <- data.frame(vertices)
        colnames(vertices) <- c("x", "y")
    } else {
        warning(paste("invalid network layout type:", layout.type))
        vertices <- data.frame(x = numeric(), y = numeric())
    }
    return(vertices)
}

#' ggplot Network
#
#' Plot network edges and vertices using ggplot
#
#' @param vertex.coords DataFrame containing vertex x,y coordinates, gene
#'  names, and optionally module labels.
#' @param edge.coords DataFrame containing edge line end coordinates
#' @param colour.modules Set to TRUE to colour nodes according to their module
#' @param node.alpha Value between 0 and 1 indicating nodes' transparency
#' @param edge.alpha Value between 0 and 1 indicating the edges' transparency
#' @param indicate.linkers Set to TRUE to visualise linker nodes as a diamond
#' @return ggplot ggplot object
ggplot.net <- function(vertex.coords, edge.coords, colour.modules, node.alpha,
                       edge.alpha, indicate.linkers) {
    if (colour.modules && "module" %in% colnames(vertex.coords)) {
        vertex.coords["module"] <- factor(vertex.coords$module)
        node.colour <- "module"
    } else {
        node.colour <- NULL
    }

    if (indicate.linkers && "linker" %in% colnames(vertex.coords)) {
        node.shape <- "linker"
    } else {
        node.shape <- NULL
    }

    node.aes <- aes_string(colour = node.colour, shape = node.shape)

    gg <- ggplot(vertex.coords, aes(x, y))
    gg <- gg + geom_point(node.aes, size = 10, alpha = node.alpha)
    gg <- gg + geom_text(aes(label = gene))
    gg <- gg + geom_line(aes(x, y, group = edge.id), data = edge.coords,
                         alpha = edge.alpha)
    gg <- gg + theme_minimal()
    gg <- gg + theme(panel.grid = element_blank(),
                     axis.ticks = element_blank(),
                     axis.line  = element_blank(),
                     axis.text  = element_blank(),
                     axis.title = element_blank())
    gg <- gg + scale_x_continuous(expand = c(0.10, 0))
    gg <- gg + scale_shape_manual(values = c(16, 18))
    return(gg)
}

#' Plot Network
#
#' Plot the interaction network.
#
#' @param object ReactomeFINetwork object
#' @param color.modules Set to FALSE to turn off module colouring
#' @param indicate.linkers Set to TRUE to visualise linker nodes as a diamond
#' @return ggplot ggplot object containing a visualization of the given
#'  network
#
#' @importFrom sna gplot.layout.adj
#' @importFrom sna gplot.layout.circle
#' @importFrom sna gplot.layout.circrand
#' @importFrom sna gplot.layout.eigen
#' @importFrom sna gplot.layout.fruchtermanreingold
#' @importFrom sna gplot.layout.geodist
#' @importFrom sna gplot.layout.hall
#' @importFrom sna gplot.layout.kamadakawai
#' @importFrom sna gplot.layout.mds
#' @importFrom sna gplot.layout.princoord
#' @importFrom sna gplot.layout.random
#' @importFrom sna gplot.layout.rmds
#' @importFrom sna gplot.layout.segeo
#' @importFrom sna gplot.layout.seham
#' @importFrom sna gplot.layout.spring
#' @importFrom sna gplot.layout.springrepulse
#' @importFrom sna gplot.layout.target
#' @import ggplot2
#
#' @export
#' @rdname plot-methods
#' @aliases plot,ReactomeFINetwork,missing-method
setMethod("plot", signature(x = "ReactomeFINetwork", y = "missing"),
           function(x, layout = "kamadakawai", colour.modules = TRUE,
                    min.module.size = 1, node.alpha = 0.5, edge.alpha = 0.25,
                    indicate.linkers = TRUE) {
    edgelist <- fis(x)

    if (min.module.size > 1 && nrow(modules(x)) > 0) {
        module.size.filter <- function(m) if (nrow(m) >= min.module.size) m
        plot.modules <- ddply(modules(x), .(module), module.size.filter)

        gene.filter <- function(e) if (all(e %in% plot.modules$gene)) e
        edgelist <- apply(edgelist, 1, gene.filter)
        edgelist <- data.frame(do.call(rbind, edgelist))
    }

    genes <- unlist(edgelist)
    genes <- genes[!duplicated(genes)]

    edgelist["first.protein"] <- factor(edgelist$first.protein, levels=genes)
    edgelist["first.protein"] <- as.numeric(edgelist$first.protein)
    edgelist["second.protein"] <- factor(edgelist$second.protein, levels=genes)
    edgelist["second.protein"] <- as.numeric(edgelist$second.protein)
    edgelist <- as.matrix(edgelist)

    adj.mat <- matrix(0, length(genes), length(genes))
    adj.mat[edgelist] <- 1

    vertex.coords <- layout.net(adj.mat, layout)
    vertex.coords["gene"] <- genes

    edge.coords1 <- vertex.coords[edgelist[, "first.protein"], ]
    edge.coords2 <- vertex.coords[edgelist[, "second.protein"], ]
    edge.coords <- rbind(edge.coords1, edge.coords2)
    edge.coords["edge.id"] <- rep(1:nrow(edgelist), 2)

    if (colour.modules && nrow(modules(x)) > 0) {
        vertex.coords <- merge(vertex.coords, modules(x))
    }

    if (indicate.linkers && any(!genes %in% genes(x))) {
        vertex.coords["linker"] <- !vertex.coords$gene %in% genes(x)
    }

    gg <- ggplot.net(vertex.coords, edge.coords, colour.modules, node.alpha,
                     edge.alpha, indicate.linkers)
    return(gg)
})

#' @rdname buildCytoscapeGraph-methods
#' @aliases buildCytoscapeGraph,ReactomeFINetwork-method
setMethod("buildCytoscapeGraph", signature("ReactomeFINetwork"),
          function(object, layout = "force-directed") {
    if (!require(RCytoscape)) stop("RCytoscape is not installed.")

    edgelist <- fis(object)
    edgelist <- edgelist[!duplicated(edgelist), ]
    genes <- unlist(edgelist)
    genes <- genes[!duplicated(genes)]

    cyto.g <- graphNEL(nodes = as.character(genes))
    for (i in 1:nrow(edgelist)) {
        cyto.g <- addEdge(edgelist[i, 1], edgelist[i, 2], cyto.g)
    }

    cy <- CytoscapeConnection()
    window.title <- deparse(substitute(object))
    if (window.title %in% as.character(getWindowList(cy))) {
        deleteWindow(cy, window.title)
    }
    cw = new.CytoscapeWindow(window.title, cyto.g)
    displayGraph(cw)
    layoutNetwork(cw, layout.name = layout)
    setNodeLabelRule(cw, "label")

    return(cw)
})
