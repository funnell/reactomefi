\docType{methods}
\name{annotateModules}
\alias{annotateModules}
\alias{annotateModules,ReactomeFINetwork-method}
\title{Annotate Network Modules}
\usage{
annotateModules(object, type = c("Pathway", "BP", "CC", "MF"),
  min.module.size = 1, include.linkers = FALSE)

\S4method{annotateModules}{ReactomeFINetwork}(object, type = c("Pathway",
  "BP", "CC", "MF"), min.module.size = 1, include.linkers = FALSE)
}
\arguments{
  \item{object}{ReactomeFINetwork object.}

  \item{type}{Character string containing the type of
  annotation to use.  Accepted values are "Pathway", "BP"
  for biological process, "CC" for cellular component, and
  "MF" for molecular function.}

  \item{include.linkers}{Set to TRUE if linker genes in the
  network should be included in network annotation. This
  may bias results. (default: FALSE)}
}
\value{
data.frame Results of the gene set enrichment analysis. The
output will be the same as the \code{\link{annotate}}
method plus another column for the module the annotation
corresponds to.
}
\description{
Perform gene set enrichment analysis on clustered network
modules with pathways or GO terms.
}

