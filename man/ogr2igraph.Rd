\name{ogr2igraph}
\alias{ogr2igraph}

\title{Integrates evolutionary rooting information and graphs.}

\description{
This function adds evolutionary rooting information to an 'igraph' class object.
}

\usage{
ogr2igraph(ogr, cogdata, g, idkey = "ENTREZ")
}

\arguments{
\item{ogr}{
a processed object of class 'OGR' evaluated by the \code{\link[geneplast:groot]{groot}} method.
}
\item{cogdata}{
a data.frame with COG to protein mapping, with at least one gene annotation type listed in the 'rtni' object (e.g. ENTREZ gene ID).
}
\item{g}{
a 'igraph' object.
}
\item{idkey}{
 a single character value specifying a vertex attribute name used to map 'cogdata' annotation to 'g'.
}
}

\value{
Return an updated 'igraph' object with evolutionary root information.
}

\author{
Rodrigo Dalmolin, Mauro Castro, Sheyla Trefflich
}


\references{
Dalmolin RJ and Castro, MA. Geneplast: Evolutionary rooting using orthologous groups distribution. Journal Paper (in preparation), 2016.
}

\seealso{
  \code{\link[geneplast:OGR-class]{OGR-class}}
}


\examples{

\dontrun{
#This example requires the geneplast.data.string.v91 package! (currently available under request)
library(geneplast.data.string.v91)
data(gpdata_string_v91)
}

}

\keyword{methods}

