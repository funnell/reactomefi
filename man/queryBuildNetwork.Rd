\name{queryBuildNetwork}
\alias{queryBuildNetwork}
\title{Query Build Network}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{genes}{Character vector of gene names.}
}
\value{
  data.frame Each row represents a functional interaction
  and comprises two columns - one for each gene in the
  interaction.
}
\description{
  Query the RESTful API for FIs between genes in a provided
  vector. Uses linker genes.
}

