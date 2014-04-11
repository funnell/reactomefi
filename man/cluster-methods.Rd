\docType{methods}
\name{cluster}
\alias{cluster}
\alias{cluster,ReactomeFINetwork-method}
\title{Cluster Network}
\usage{
cluster(object)

\S4method{cluster}{ReactomeFINetwork}(object)
}
\arguments{
  \item{object}{ReactomeFINetwork object.}
}
\value{
ReactomeFINetwork ReactomeFINetwork object with module
attribute set
}
\description{
Cluster FI network using its functional interaction data.
This method uses the spectral partition based network
clustering algorithm by M.E. Newman.
}

