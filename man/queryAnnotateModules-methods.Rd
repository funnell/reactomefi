\docType{methods}
\name{queryAnnotateModules}
\alias{queryAnnotateModules}
\alias{queryAnnotateModules,ReactomeFIService,data.frame,character-method}
\title{Query Annotate FI Network Module Gene Set}
\usage{
queryAnnotateModules(object, module.nodes, type = c("Pathway", "BP", "CC",
  "MF"))

\S4method{queryAnnotateModules}{ReactomeFIService,data.frame,character}(object,
  module.nodes, type = c("Pathway", "BP", "CC", "MF"))
}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{module.nodes}{Data frame with network nodes (genes)
  and their module.}

  \item{type}{Gene annotation enrichment type (Pathway, BP,
  CC, MF).}
}
\value{
data.frame Each row represents an annotation of the
provided type and includes related information such as the
p-value and FDR generated from enrichment analysis and the
module the annotation belongs to.
}
\description{
Query the RESTful API to annotate a gene set from a FI
network module with enriched pathways, or GO terms
}

