\name{annotateModules}
\alias{annotateModules}
\title{Annotate FI Network Module Gene Set}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{module.nodes}{Data frame with network nodes (genes)
  and their module.}
}
\value{
  data.frame Each row represents an annotation of the
  provided type and includes related information such as
  the p-value and FDR generated from enrichment analysis
  and the module the annotation belongs to.
}
\description{
  Annotate a gene set from a FI network module with
  enriched pathways, or GO terms
}

