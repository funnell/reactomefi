\name{fis2str}
\alias{fis2str}
\title{FIs to String}
\usage{
  fis2str(fis)
}
\arguments{
  \item{fis}{data frame of FIs (each row contains two gene
  names)}
}
\value{
  character Character string in TSV format where rows are
  separated by "\\n" and columns within rows are separated
  by "\\t"
}
\description{
  Convert a FI data frame into a string according to
  conventions used in the ReactomeFI API.
}

