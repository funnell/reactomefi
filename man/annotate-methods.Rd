\docType{methods}
\name{annotate}
\alias{annotate}
\alias{annotate,ReactomeFINetwork,character-method}
\title{Annotate Network}
\usage{
annotate(object, type = c("Pathway", "BP", "CC", "MF"),
  include.linkers = FALSE)

\S4method{annotate}{ReactomeFINetwork,character}(object, type = c("Pathway",
  "BP", "CC", "MF"), include.linkers = FALSE)
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
data.frame Results of the gene set enrichment analysis.
}
\description{
Perform gene set enrichment analysis with pathways or GO
terms.
}

