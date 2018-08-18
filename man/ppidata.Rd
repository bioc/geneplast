\name{ppi}
\alias{ppi.gs}
\docType{data}

\title{A pre-processed igraph object for the geneplast package.}

\description{
An igraph object used to demonstrate geneplast main functions.
}

\usage{data(ppi.gs)}

\format{
  \code{ppi.gs}{an igraph object.}
}

\details{
Protein-protein interactions (PPI) mapped to H. sapiens apoptosis and genome-stability genes. The PPI information is derived from the STRING database, release 9.1, using combined score >=900 (Franceschini et al., 2013).
}

\value{
an igraph object.
}

\references{
Franceschini A, Szklarczyk D, Frankild S, Kuhn M, Simonovic M, Roth A, Lin J, Minguez P, Bork P, von Mering C, Jensen LJ. STRING v9.1: protein-protein interaction networks, with increased coverage and integration. Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29.
}

\examples{
  data(ppi.gs)
}

\keyword{dataset}

