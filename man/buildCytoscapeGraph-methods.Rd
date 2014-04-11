\docType{methods}
\name{buildCytoscapeGraph}
\alias{buildCytoscapeGraph}
\alias{buildCytoscapeGraph,ReactomeFINetwork-method}
\title{Build Cytoscape Graph}
\usage{
buildCytoscapeGraph(object, layout = "force-directed")

\S4method{buildCytoscapeGraph}{ReactomeFINetwork}(object,
  layout = "force-directed")
}
\arguments{
  \item{object}{ReactomeFINetwork object.}

  \item{layout}{Cytoscape network layout method.}
}
\value{
CytoscapeWindowClass Cytoscape window object.
}
\description{
Construct and visualize a ReactomeFINetwork object in
Cytoscape.
}

