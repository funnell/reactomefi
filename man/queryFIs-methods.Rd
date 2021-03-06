\docType{methods}
\name{queryFIs}
\alias{queryFIs}
\alias{queryFIs,ReactomeFIService,character-method}
\title{Query Functional Interactions}
\usage{
queryFIs(object, genes)

\S4method{queryFIs}{ReactomeFIService,character}(object, genes)
}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{genes}{Character vector of gene names.}
}
\value{
data.frame Each row represents a functional interaction and
comprises two columns - one for each gene in the
interaction.
}
\description{
Query the RESTful API for FIs between genes in a provided
vector.
}

