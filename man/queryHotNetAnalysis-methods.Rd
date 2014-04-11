\docType{methods}
\name{queryHotNetAnalysis}
\alias{queryHotNetAnalysis}
\alias{queryHotNetAnalysis,ReactomeFIService-method}
\title{Query HotNet Analysis}
\usage{
queryHotNetAnalysis(object, gene.scores, delta, fdr, permutations, auto.delta)

\S4method{queryHotNetAnalysis}{ReactomeFIService}(object, gene.scores, delta,
  fdr, permutations, auto.delta)
}
\arguments{
  \item{object}{ReactomeFIService object.}

  \item{gene.scores}{Data frame containing gene, score
  pairs.}

  \item{delta}{Numeric delta value.}

  \item{fdr}{FDR cutoff.}

  \item{permutations}{Number of permutations.}
}
\value{
data.frame
}
\description{
Query the RESTful API to do HotNet analysis. The ReactomeFI
API implements the "HotNet" algorithm for doing cancer
mutation analysis developed by Raphael's group at Brown
University.
}

