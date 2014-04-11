\docType{methods}
\name{build}
\alias{build}
\alias{build,ReactomeFINetwork,character-method}
\title{Build Network}
\usage{
build(object, genes, use.linkers = FALSE)

\S4method{build}{ReactomeFINetwork,character}(object, genes,
  use.linkers = FALSE)
}
\arguments{
  \item{object}{ReactomeFINetwork object.}

  \item{genes}{Character vector of gene names}

  \item{use.linkers}{Set to TRUE to build a network using
  linker genes (default: FALSE)}
}
\value{
ReactomeFINetwork ReactomeFINetwork object with fis
attribute set
}
\description{
Build FI network from a list of genes.
}

