\name{getPostXML}
\alias{getPostXML}
\title{Get POST Query XML Response
Make a POST request to a URL with a query string to be provided to
curlPerform's postfields argument. curlPerform is used instead of postForm
to allow sending unkeyed data.}
\usage{
getPostXML(url, body)
}
\arguments{
  \item{url}{Character string containing URL to post query
  to.}

  \item{body}{Character string containing data to send in
  the post request.}
}
\value{
XMLDocument XML document generated from the POST query
response.
}
\description{
Get POST Query XML Response Make a POST request to a URL
with a query string to be provided to curlPerform's
postfields argument. curlPerform is used instead of
postForm to allow sending unkeyed data.
}

