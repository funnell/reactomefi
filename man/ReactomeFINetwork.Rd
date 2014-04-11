\docType{class}
\name{ReactomeFINetwork-class}
\alias{ReactomeFINetwork}
\alias{ReactomeFINetwork-class}
\title{ReactomeFINetwork Class}
\usage{
ReactomeFINetwork(version = c("2009", "2012", "2013"), genes = character(),
  use.linkers = FALSE, cluster = FALSE)
}
\arguments{
  \item{version}{Version of ReactomeFI network (2009, 2012,
  2013).}

  \item{genes}{Character vector of gene names.}

  \item{use.linkers}{Set to TRUE to build a network using
  linker genes (default: FALSE)}

  \item{cluster}{Set to TRUE to cluster the network if gene
  names have been provided (default FALSE).}
}
\value{
ReactomeFINetwork ReactomeFINetwork S4 object
}
\description{
Represents a functional interaction network generated using
the Reactome database and a given list of genes.

ReactomeFINetwork constructor
}

