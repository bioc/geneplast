\name{ogr2tni}
\alias{ogr2tni}

\title{Integrates evolutionary rooting information and regulatory networks.}

\description{
This function adds evolutionary rooting information to regulons of a TNI class object.
}

\usage{
ogr2tni(ogr, cogdata, tni)
}

\arguments{
\item{ogr}{
a processed object of class 'OGR' evaluated by the \code{\link[geneplast:groot]{groot}} method.
}
\item{cogdata}{
a data.frame with COG to protein mapping, with at least one gene annotation type listed in the 'rtni' object (e.g. ENTREZ gene ID).
}
\item{tni}{
a processed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the methods 
\code{\link[RTN:tni.permutation]{tni.permutation}}, \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} 
and \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}}.
}
}

\value{
Return an updated TNI object with evolutionary root information.
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

# Evolutionary rooting inference for all human genes

library(geneplast.data.string.v91)
data(gpdata_string_v91)

cogids <- cogdata$cog_id[cogdata$ssp_id=="9606"]
ogr <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogids)

ogr <- groot(ogr, nPermutations=100, verbose=TRUE)

# Integrates evolutionary rooting information and regulons

library(RTN)
library(Fletcher2013b)
data("rtni1st")
rtni1st <- ogr2tni(ogr, cogdata, rtni1st)


}

}

\keyword{methods}

