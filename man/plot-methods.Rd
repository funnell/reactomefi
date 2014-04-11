\docType{methods}
\name{plot,ReactomeFINetwork,missing-method}
\alias{plot,ReactomeFINetwork,missing-method}
\title{Plot Network
Plot the interaction network. Layout methods are currently provided by the
 gplot.layout methods in the sna package.}
\usage{
\S4method{plot}{ReactomeFINetwork,missing}(x, layout = "kamadakawai",
  colour.modules = TRUE, min.module.size = 1, node.alpha = 0.5,
  edge.alpha = 0.25, indicate.linkers = TRUE)
}
\arguments{
  \item{x}{ReactomeFINetwork object}

  \item{layout}{Name of a layout method.}

  \item{colour.modules}{Set to FALSE to turn off module
  colouring.}

  \item{min.module.size}{Minimum size of a module to be
  plotted (default: 1).}

  \item{node.alpha}{Alpha value for nodes (default: 0.25).}

  \item{edge.alpha}{Alpha value for edges (default: 0.25).}

  \item{indicate.linkers}{Set to TRUE to visualise linker
  nodes as a diamond.}
}
\value{
ggplot ggplot object containing a visualization of the
given network.
}
\description{
Plot Network Plot the interaction network. Layout methods
are currently provided by the gplot.layout methods in the
sna package.
}

