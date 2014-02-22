\name{ReactomeFINetwork}
\alias{ReactomeFINetwork}
\title{ReactomeFINetwork Class}
\usage{
  ReactomeFINetwork(version = c("2009", "2012"),
    genes = character(), cluster = FALSE)
}
\arguments{
  \item{version}{Version of ReactomeFI network (2009 or
  2012).}

  \item{genes}{Character vector of gene names.}

  \item{cluster}{Set to TRUE to cluster the network if gene
  names have been provided (default FALSE).}
}
\value{
  ReactomeFINetwork ReactomeFINetwork S4 object
}
\description{
  Represents a functional interaction network generated
  using the Reactome database and a given list of genes.

  ReactomeFINetwork constructor
}

