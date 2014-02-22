\docType{methods}
\name{queryEdge}
\alias{queryEdge}
\alias{queryEdge,ReactomeFIService,character,character-method}
\title{Query Edge}
\arguments{
  \item{object}{ReactomeFIService object}

  \item{name1}{Name of the first gene in the interaction.}

  \item{name2}{Name of the second gene in the interaction.}
}
\value{
  data.frame Each of the two rows corresponds to a protein
  specified in the input parameters. The columns contain
  information regarding the proteins including accession
  ID, database name, protein name and sequence.
}
\description{
  Query detailed information for a network edge.
}

