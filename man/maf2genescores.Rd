\name{maf2genescores}
\alias{maf2genescores}
\title{NCI MAF to genescore}
\usage{
  maf2genescores(maf)
}
\arguments{
  \item{maf}{NCI Mutation Annotation File.}
}
\value{
  data.frame Two column data.frame containing gene names
  and scores.
}
\description{
  Converts a NCI Mutation Annotation File to a two column
  data.frame where the first column contains gene names and
  the second column contains scores representing the
  proportion of samples in which a gene is mutated.
}

