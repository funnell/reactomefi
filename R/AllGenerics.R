##
# ReactomeFIService generics
##

#' Reactome FI Network Version
#' 
#' Retrieve the Reactome FI network version.
#' 
#' @param object ReactomeFIService object.
#' @return character
#'
#' @docType methods
#' @rdname version-methods
setGeneric("version", function(object) standardGeneric("version"))

#' RESTful API URL
#'
#' Retrieve the base URL for the ReactomeFI RESTful API.
#'
#' @param object ReactomeFIService object.
#' @return character
#'
#' @docType methods
#' @rdname serviceURL-methods
setGeneric("serviceURL", function(object) standardGeneric("serviceURL"))

#' Query Functional Interactions
#'
#' Query the RESTful API for FIs between genes in a provided vector.
#'
#' @param object ReactomeFIService object.
#' @param genes Character vector of gene names.
#' @return data.frame Each row represents a functional interaction and
#'  comprises two columns - one for each gene in the interaction.
#'
#' @docType methods
#' @rdname queryFIs-methods
setGeneric("queryFIs", function(object, genes) standardGeneric("queryFIs"))

#' Query Build Network
#'
#' Query the RESTful API for FIs between genes in a provided vector. Uses
#' linker genes.
#'
#' @param object ReactomeFIService object.
#' @param genes Character vector of gene names.
#' @return data.frame Each row represents a functional interaction and
#'  comprises two columns - one for each gene in the interaction.
#'
#' @docType methods
#' @rdname queryBuildNetwork-methods
setGeneric("queryBuildNetwork", function(object, genes) {
    standardGeneric("queryBuildNetwork")
})

#' Query FIs Between Genes
#'
#' Query FIs between a list of pairs of genes.
#'
#' @param object ReactomeFIService object
#' @param gene.pairs Data frame or matrix of gene pairs in which to look for
#'  FIs. Each row contains two columns - one for each gene in the pair.
#' @return data.frame Each row represents a functional interaction and
#'  comprises two columns - one for each gene in the interaction.
#'
#' @docType methods
#' @rdname queryFIsBetween-methods
setGeneric("queryFIsBetween", function(object, gene.pairs) {
    standardGeneric("queryFIsBetween")
})

#' Query Edge
#'
#' Query detailed information for a network edge.
#'
#' @param object ReactomeFIService object
#' @param name1 Name of the first gene in the interaction.
#' @param name2 Name of the second gene in the interaction.
#' @return data.frame Each of the two rows corresponds to a protein specified
#'  in the input parameters. The columns contain information regarding the
#'  proteins including accession ID, database name, protein name and sequence.
#'
#' @docType methods
#' @rdname queryEdge-methods
setGeneric("queryEdge", function(object, name1, name2) {
    standardGeneric("queryEdge")
})

#' Query Cluster Functional Interaction Network
#'
#' Query the RESTful API to cluster a FI network. The network nodes (genes)
#' and the network modules they belong to are returned.
#'
#' @param object ReactomeFIService object.
#' @param fis Data frame of functional interactions. Should be two columns per
#'  row indicating the two nodes in the interaction
#' @return data.frame Each row contains the gene name and the module id.
#'
#' @docType methods
#' @rdname queryCluster-methods
setGeneric("queryCluster", function(object, fis) {
    standardGeneric("queryCluster")
})

#' Query Annotate Gene Set
#'
#' Query the RESTful API to annotate a gene set with enriched pathways, or GO
#'  terms
#'
#' @param object ReactomeFIService object.
#' @param genes Character vector of gene names.
#' @param type Gene annotation enrichment type (Pathway, BP, CC, MF).
#' @return data.frame Each row represents an annotation of the provided type
#'  and includes related information such as the p-value and FDR generated
#'  from enrichment analysis.
#'
#' @docType methods
#' @rdname queryAnnotateGeneSet-methods
setGeneric("queryAnnotateGeneSet",
           function(object, genes, type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("queryAnnotateGeneSet")
})

#' Query Annotate FI Network Module Gene Set
#'
#' Query the RESTful API to annotate a gene set from a FI network module with
#'  enriched pathways, or GO terms
#'
#' @param object ReactomeFIService object.
#' @param module.nodes Data frame with network nodes (genes) and their module.
#' @param type Gene annotation enrichment type (Pathway, BP, CC, MF).
#' @return data.frame Each row represents an annotation of the provided type
#'  and includes related information such as the p-value and FDR generated
#'  from enrichment analysis and the module the annotation belongs to.
#'
#' @docType methods
#' @rdname queryAnnotateModules-methods
setGeneric("queryAnnotateModules",
           function(object, module.nodes,
                    type = c("Pathway", "BP", "CC", "MF")) {
    standardGeneric("queryAnnotateModules")
})

#' Query HotNet Analysis
#'
#' Query the RESTful API to do HotNet analysis. The ReactomeFI API implements
#'  the "HotNet" algorithm for doing cancer mutation analysis developed by
#'  Raphael's group at Brown University.
#'
#' @param object ReactomeFIService object.
#' @param gene.scores Data frame containing gene, score pairs.
#' @param delta Numeric delta value.
#' @param fdr FDR cutoff.
#' @param permutations Number of permutations.
#' @param auto.delta Set to TRUE to automatically select a delta value.
#' @return data.frame
#'
#' @docType methods
#' @rdname queryHotNetAnalysis-methods
setGeneric("queryHotNetAnalysis",
           function(object, gene.scores, delta, fdr, permutations,
                    auto.delta) {
    standardGeneric("queryHotNetAnalysis")
})


##
# Network generics
##

#' Retrieve ReactomeFIService Object
#'
#' Retrieve the underlying ReactomeFIService object.
#'
#' @param object ReactomeFINetwork or HotNet object.
#' @return ReactomeFIService
#'
#' @export
#' @docType methods
#' @rdname service-methods
setGeneric("service", function(object)  standardGeneric("service"))

#' Retrieve Network Genes
#'
#' Retreive network gene list
#'
#' @param object ReactomeFINetwork.
#' @return character
#'
#' @export
#' @docType methods
#' @rdname genes-methods
setGeneric("genes", function(object) standardGeneric("genes"))

#' Set Network Genes
#'
#' Set network gene list.
#'
#' @param value Character gene list
#'
#' @export
#' @docType methods
#' @rdname genes-methods
setGeneric("genes<-", function(object, value) standardGeneric("genes<-"))

#' Retrieve Network Module Data
#'
#' Retreive network module data.
#'
#' @param object ReactomeFINetwork or HotNet object.
#' @return data.frame
#'
#' @export
#' @docType methods
#' @rdname modules-methods
setGeneric("modules", function(object) standardGeneric("modules"))

#' Set Network Module Data
#'
#' Set network module data. ReactomeFINetwork objects accept data.frames and
#'  HotNet objects accept lists
#'
#' @param value Network module information
#'
#' @export
#' @docType methods
#' @rdname modules-methods
setGeneric("modules<-", function(object, value) standardGeneric("modules<-"))


##
# ReactomeFINetwork generics
##

#' Retrieve FI Network Data
#'
#' Retrieve the stored FI network interactions.
#'
#' @param object ReactomeFINetwork object.
#' @return data.frame FI network data. Each row represents an interaction
#'  between two genes and consists of two columns - one for each gene.
#'
#' @export
#' @docType methods
#' @rdname fis-methods
setGeneric("fis", function(object) standardGeneric("fis"))

#' Set FI Network Data
#'
#' Set the FI Network Data.
#'
#' @param value Data frame representing the Reactome FI network data. Each row
#'  represents an interaction between two genes and consists of two columns -
#'  one for each gene.
#' @return ReactomeFINetwork
#'
#' @export
#' @docType methods
#' @rdname fis-methods
setGeneric("fis<-", function(object, value) standardGeneric("fis<-"))

#' Build Network
#'
#' Build FI network from a list of genes.
#'
#' @param object ReactomeFINetwork object.
#' @param genes Character vector of gene names
#' @param use.linkers Set to TRUE to build a network using linker genes
#'  (default: FALSE)
#' @return ReactomeFINetwork ReactomeFINetwork object with fis attribute set
#'
#' @export
#' @docType methods
#' @rdname build-methods
setGeneric("build",
           function(object, genes, use.linkers = FALSE) {
    standardGeneric("build")
})

#' Cluster Network
#'
#' Cluster FI network using its functional interaction data. This method uses
#'  the spectral partition based network clustering algorithm by M.E. Newman.
#'
#' @param object ReactomeFINetwork object.
#' @return ReactomeFINetwork ReactomeFINetwork object with module attribute
#'  set
#'
#' @export
#' @docType methods
#' @rdname cluster-methods
setGeneric("cluster", function(object) standardGeneric("cluster"))

#' Annotate Network
#'
#' Perform gene set enrichment analysis with pathways or GO terms.
#'
#' @param object ReactomeFINetwork object.
#' @param type Character string containing the type of annotation to use.
#'  Accepted values are "Pathway", "BP" for biological process, "CC" for
#'  cellular component, and "MF" for molecular function.
#' @param include.linkers Set to TRUE if linker genes in the network should be
#'  included in network annotation. This may bias results. (default: FALSE)
#' @return data.frame Results of the gene set enrichment analysis.
#'
#' @export
#' @docType methods
#' @rdname annotate-methods
setGeneric("annotate",
           function(object, type = c("Pathway", "BP", "CC", "MF"),
                    include.linkers = FALSE) {
    standardGeneric("annotate")
})

#' Annotate Network Modules
#'
#' Perform gene set enrichment analysis on clustered network modules with
#'  pathways or GO terms.
#'
#' @param object ReactomeFINetwork object.
#' @param type Character string containing the type of annotation to use.
#'  Accepted values are "Pathway", "BP" for biological process, "CC" for
#'  cellular component, and "MF" for molecular function.
#' @param min.module.size Minimum module size to consider for annotation
#'  (default: 1).
#' @param include.linkers Set to TRUE if linker genes in the network should be
#'  included in network annotation. This may bias results. (default: FALSE)
#' @return data.frame Results of the gene set enrichment analysis. The output
#'  will be the same as the \code{\link{annotate}} method plus another column
#'  for the module the annotation corresponds to.
#'
#' @importFrom plyr ddply
#' @importFrom plyr .
#'
#' @export
#' @docType methods
#' @rdname annotateModules-methods
setGeneric("annotateModules",
           function(object, type = c("Pathway", "BP", "CC", "MF"),
                    min.module.size = 1, include.linkers = FALSE) {
    standardGeneric("annotateModules")
})

#' Build Cytoscape Graph
#'
#' Construct and visualize a ReactomeFINetwork object in Cytoscape.
#'
#' @param object ReactomeFINetwork object.
#' @param layout Cytoscape network layout method.
#' @return CytoscapeWindowClass Cytoscape window object.
#'
#' @export
#' @docType methods
#' @rdname buildCytoscapeGraph-methods
setGeneric("buildCytoscapeGraph", function(object, layout = "force-directed") {
    standardGeneric("buildCytoscapeGraph")
})


##
# HotNetAnalysis generics
##

#' Create HotNet ReactomeFI Network
#'
#' Create a ReactomeFI Network from a subset of the HotNet modules.
#'
#' @param object HotNet object.
#' @param fdr False discovery rate threshold for HotNet modules
#' @return ReactomeFINetwork
#'
#' @export
#' @docType methods
#' @rdname subnet-methods
setGeneric("subnet", function(object, fdr = 0.05) standardGeneric("subnet"))
