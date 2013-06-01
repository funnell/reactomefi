\name{annotateGeneSet}
\alias{annotateGeneSet}
\title{Annotate Gene Set}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{genes}{Character vector of gene names.}
}
\value{
  data.frame Each row represents an annotation of the
  provided type and includes related information such as
  the p-value and FDR generated from enrichment analysis.
}
\description{
  Annotate a gene set with enriched pathways, or GO terms
}

