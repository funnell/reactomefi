\docType{methods}
\name{queryAnnotateGeneSet}
\alias{queryAnnotateGeneSet}
\alias{queryAnnotateGeneSet,ReactomeFIService,character,character-method}
\title{Query Annotate Gene Set}
\usage{
queryAnnotateGeneSet(object, genes, type = c("Pathway", "BP", "CC", "MF"))

\S4method{queryAnnotateGeneSet}{ReactomeFIService,character,character}(object,
  genes, type = c("Pathway", "BP", "CC", "MF"))
}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{genes}{Character vector of gene names.}
}
\value{
data.frame Each row represents an annotation of the
provided type and includes related information such as the
p-value and FDR generated from enrichment analysis.
}
\description{
Query the RESTful API to annotate a gene set with enriched
pathways, or GO terms
}

