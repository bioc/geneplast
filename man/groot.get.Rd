\name{groot.get}
\alias{groot.get}

\title{
Get information from individual slots in an OGR object.
}

\description{
Get information from individual slots in an OGR object and any available results from a previous analysis.
}

\usage{
groot.get(object, what="status")
}

\arguments{
\item{object}{
an object of class 'OGR' \code{\link[geneplast:OGR-class]{OGR-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. Options: "cogids", "spbranches", "orthoroot" ,"results" and "status".
}
}

\value{
slot content from an object of class 'OGR' \code{\link[geneplast:OGR-class]{OGR-class}}.
}

\author{
Rodrigo Dalmolin, Mauro Castro
}

\examples{

#load datasets used for demonstration
data(gpdata.gs)

#create and object of class 'OGR' for H. sapiens
ogr <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogids)

## run the groot function
## this example uses the orthologous groups listed in the gpdata object
ogr <- groot(ogr, nPermutations=100)
res <- groot.get(ogr, what="results")

}

\keyword{methods}
