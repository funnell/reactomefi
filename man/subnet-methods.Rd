\docType{methods}
\name{subnet}
\alias{subnet}
\alias{subnet,HotNet-method}
\title{Create HotNet ReactomeFI Network}
\usage{
subnet(object, fdr = 0.05)

\S4method{subnet}{HotNet}(object, fdr = 0.05)
}
\arguments{
  \item{object}{HotNet object.}

  \item{fdr}{False discovery rate threshold for HotNet
  modules}
}
\value{
ReactomeFINetwork
}
\description{
Create a ReactomeFI Network from a subset of the HotNet
modules.
}

