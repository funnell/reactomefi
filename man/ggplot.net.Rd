\name{ggplot.net}
\alias{ggplot.net}
\title{ggplot Network}
\usage{
  ggplot.net(vertex.coords, edge.coords, colour.modules,
    node.alpha, edge.alpha)
}
\arguments{
  \item{vertex.coords}{DataFrame containing vertex x,y
  coordinates, gene names, and optionally module labels.}

  \item{edge.coords}{DataFrame containing edge line end
  coordinates}

  \item{colour.modules}{Set to TRUE to colour nodes
  according to their module}

  \item{node.alpha}{Value between 0 and 1 indicating nodes'
  transparency}

  \item{edge.alpha}{Value between 0 and 1 indicating the
  edges' transparency}
}
\value{
  ggplot ggplot object
}
\description{
  Plot network edges and vertices using ggplot
}

