\name{df2tsv}
\alias{df2tsv}
\title{Data Frame to TSV}
\usage{
  df2tsv(dat)
}
\arguments{
  \item{dat}{Data frame to be converted to a TSV string.}
}
\value{
  character Each row is separated by "\\n" and each column
  within a row is separated by a "\\t".
}
\description{
  Convert a data frame (not including headers) into a TSV
  string.
}

