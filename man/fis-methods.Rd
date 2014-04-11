\docType{methods}
\name{fis}
\alias{fis}
\alias{fis,ReactomeFINetwork-method}
\alias{fis<-}
\alias{fis<-,ReactomeFINetwork,data.frame-method}
\title{Retrieve FI Network Data}
\usage{
fis(object)

fis(object) <- value

\S4method{fis}{ReactomeFINetwork}(object)

\S4method{fis}{ReactomeFINetwork,data.frame}(object) <- value
}
\arguments{
  \item{object}{ReactomeFINetwork object.}

  \item{object}{ReactomeFINetwork object.}

  \item{value}{Data frame representing the Reactome FI
  network data. Each row represents an interaction between
  two genes and consists of two columns - one for each
  gene.}
}
\value{
data.frame FI network data. Each row represents an
interaction between two genes and consists of two columns -
one for each gene.

ReactomeFINetwork
}
\description{
Retrieve the stored FI network interactions.

Set the FI Network Data.
}

