\docType{methods}
\name{modules}
\alias{modules}
\alias{modules,HotNet-method}
\alias{modules,ReactomeFINetwork-method}
\alias{modules<-}
\alias{modules<-,ReactomeFINetwork,data.frame-method}
\title{Retrieve Network Module Data}
\usage{
modules(object)

modules(object) <- value

\S4method{modules}{HotNet}(object)

\S4method{modules}{ReactomeFINetwork}(object)

\S4method{modules}{ReactomeFINetwork,data.frame}(object) <- value
}
\arguments{
  \item{object}{ReactomeFINetwork or HotNet object.}

  \item{object}{ReactomeFINetwork or HotNet object.}

  \item{value}{Network module information}
}
\value{
data.frame
}
\description{
Retreive network module data.

Set network module data. ReactomeFINetwork objects accept
data.frames and HotNet objects accept lists
}

