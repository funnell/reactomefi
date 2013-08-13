\name{queryFIsBetween}
\alias{queryFIsBetween}
\title{Query FIs Between Genes}
\arguments{
  \item{object}{ReactomeFIService object}

  \item{gene.pairs}{Data frame or matrix of gene pairs in
  which to look for FIs. Each row contains two columns -
  one for each gene in the pair.}
}
\value{
  data.frame Each row represents a functional interaction
  and comprises two columns - one for each gene in the
  interaction.
}
\description{
  Query FIs between a list of pairs of genes.
}

