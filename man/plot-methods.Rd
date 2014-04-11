\docType{methods}
\name{plot,ReactomeFINetwork,missing-method}
\alias{plot,ReactomeFINetwork,missing-method}
\title{Plot Network
Plot the interaction network.}
\usage{
\S4method{plot}{ReactomeFINetwork,missing}(x, layout = "kamadakawai",
  colour.modules = TRUE, min.module.size = 1, node.alpha = 0.5,
  edge.alpha = 0.25, indicate.linkers = TRUE)
}
\arguments{
  \item{object}{ReactomeFINetwork object}

  \item{color.modules}{Set to FALSE to turn off module
  colouring}

  \item{indicate.linkers}{Set to TRUE to visualise linker
  nodes as a diamond}

  \item{x}{the coordinates of points in the plot.
  Alternatively, a single plotting structure, function or
  \emph{any \R object with a \code{plot} method} can be
  provided.}
}
\value{
ggplot ggplot object containing a visualization of the
given network
}
\description{
Plot Network Plot the interaction network.
}

