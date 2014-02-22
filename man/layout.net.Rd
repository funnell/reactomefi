\name{layout.net}
\alias{layout.net}
\title{Layout ReactomeFINetwork}
\usage{
  layout.net(adj.mat, layout)
}
\arguments{
  \item{adj.mat}{Adjacency matrix}

  \item{layout}{Layout algorithm as defined by sna's
  gplot.layout}
}
\value{
  data.frame DataFrame containing x and y coordinates of
  network vertices
}
\description{
  Retrieve the coordinates for the network vertices
  according to the specified layout algorithm
}

