\name{gplast.get}
\alias{gplast.get}

\title{
Get information from individual slots in an OGP object.
}

\description{
Get information from individual slots in an OGP object and any available results from previous analysis.
}

\usage{
gplast.get(object, what="status")
}

\arguments{
\item{object}{
an object of class 'OGP' \code{\link[geneplast:OGP-class]{OGP-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. Options: "cogids", "sspids", "orthodist", "results", and "status".
}
}

\value{
slot content from an object of class 'OGP' \code{\link[geneplast:OGP-class]{OGP-class}}.
}

\author{
Rodrigo Dalmolin, Mauro Castro
}

\examples{

#load datasets used for demonstration
data(gpdata.gs)

#create and object of class 'OGP'
ogp <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogids)

## run the gplast function
## this example uses the especies/COGs listed in the gpdata object
ogp <- gplast(ogp)
res <- gplast.get(ogp,what="results")

}

\keyword{methods}
