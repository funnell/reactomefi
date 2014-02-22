\docType{methods}
\name{annotate}
\alias{annotate}
\alias{annotate,ReactomeFINetwork,character-method}
\title{Annotate Network}
\arguments{
  \item{object}{ReactomeFINetwork object.}

  \item{type}{Character string containing the type of
  annotation to use.  Accepted values are "Pathway", "BP"
  for biological process, "CC" for cellular component, and
  "MF" for molecular function.}
}
\value{
  data.frame Results of the gene set enrichment analysis.
}
\description{
  Perform gene set enrichment analysis with pathways or GO
  terms.
}

