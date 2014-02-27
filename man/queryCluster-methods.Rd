\docType{methods}
\name{queryCluster}
\alias{queryCluster}
\alias{queryCluster,ReactomeFIService,data.frame-method}
\title{Query Cluster Functional Interaction Network}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{fis}{Data frame of functional interactions. Should
  be two columns per row indicating the two nodes in the
  interaction}
}
\value{
  data.frame Each row contains the gene name and the module
  id.
}
\description{
  Query the RESTful API to cluster a FI network. The
  network nodes (genes) and the network modules they belong
  to are returned.
}
