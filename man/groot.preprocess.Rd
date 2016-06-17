\name{groot.preprocess}
\alias{groot.preprocess}

\title{Evolutionary rooting inference.}

\description{
Constructor for the 'OGR-class' object.
}

\usage{
groot.preprocess(cogdata, phyloTree, spid, cogids=NULL, verbose=TRUE)
}

\arguments{
\item{cogdata}{
a data frame with COG's data.
}
\item{phyloTree}{
an object of class "phylo".
}
\item{spid}{
a single character or integer value specifying the reference species to be used in the rooting algorithm. This species should be listed in the 'phyloTree'.
}
\item{cogids}{
an optional data frame with COG's annotation. Alternatively, it can be a character vector with COG IDs.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\details{
This function creates an \code{\link[geneplast:OGR-class]{OGR-class}} object and checks the consistency of the input data for the evolutionary root pipeline. Internally, the function access the presence and absence of orthologs for each species in a given OG and computes the \code{orthoct} data.frame.
}

\value{
A preprocessed object of class 'OGR'.
}

\author{
Rodrigo Dalmolin, Mauro Castro
}


\references{
Dalmolin RJ and Castro, MA. Geneplast: Evolutionary rooting using orthologous groups distribution. Journal Paper (in preparation), 2016.
}

\seealso{
  \code{\link[geneplast:OGR-class]{OGR-class}}
}


\examples{

#load datasets used for demonstration
data(gpdata.gs)

#create and object of class 'OGR' for H. sapiens
ogr <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogids)

}

\keyword{methods}

