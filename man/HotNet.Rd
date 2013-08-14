\name{HotNet}
\alias{HotNet}
\title{HotNet Class}
\usage{
  HotNet(version = c("2009", "2012"))
}
\arguments{
  \item{version}{Version of ReactomeFI network (2009 or
  2012).}

  \item{genescores}{data.frame where the first column
  contains gene names and the second column contains scores
  representing the proportion of samples in which the gene
  is mutated.}

  \item{delta}{HotNet delta value}

  \item{fdr.threshold}{}

  \item{permutations}{Number of permutations. Largest value
  is 1000.}

  \item{auto.delta}{If true, algorithm will select a delta
  value. This option will make the analysis take more time
  to finish.}
}
\value{
  HotNet HotNet S4 object containing the results of HotNet
  analysis. These results should be filtered before using
  them to generate a network
}
\description{
  Represents the results of a HotNet analysis

  HotNet constructor
}

