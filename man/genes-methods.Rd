\docType{methods}
\name{genes}
\alias{genes}
\alias{genes,ReactomeFINetwork-method}
\alias{genes<-}
\alias{genes<-,ReactomeFINetwork,character-method}
\title{Retrieve Network Genes}
\usage{
genes(object)

genes(object) <- value

\S4method{genes}{ReactomeFINetwork}(object)

\S4method{genes}{ReactomeFINetwork,character}(object) <- value
}
\arguments{
  \item{object}{ReactomeFINetwork.}

  \item{value}{Character gene list}
}
\value{
character
}
\description{
Retreive network gene list

Set network gene list.
}

