\docType{methods}
\name{annotateModules}
\alias{annotateModules}
\alias{annotateModules,ReactomeFINetwork-method}
\title{Annotate Network Modules}
\arguments{
  \item{object}{ReactomeFINetwork object.}

  \item{type}{Character string containing the type of
  annotation to use.  Accepted values are "Pathway", "BP"
  for biological process, "CC" for cellular component, and
  "MF" for molecular function.}
}
\value{
  data.frame Results of the gene set enrichment analysis.
  The output will be the same as the \code{\link{annotate}}
  method plus another column for the module the annotation
  corresponds to.
}
\description{
  Perform gene set enrichment analysis on clustered network
  modules with pathways or GO terms.
}

