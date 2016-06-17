\name{gpdata}
\alias{gpdata.gs}
\docType{data}

\title{A pre-processed dataset for the geneplast package.}

\description{
A dataset used to demonstrate geneplast main functions.
}

\usage{data(gpdata.gs)}

\format{
  \code{gpdata.gs}{ data frame containing information about clusters of orthologous groups.}
}

\details{

The dataset consists of 4 R objects used in the geneplast vignettes.

\describe{
  \item{cogdata}{
  A data frame with three columns containing the Clusters of Orthologous Groups information retrieved
  from the STRING database (http://string-db.org/), release 9.1. 
  Column 1 = Ensembl protein ID; column 2 = NCBI species ID;column 3 = OG ID.
  }
  \item{sspids}{
  A data frame containing the species listed in STRING database (http://string-db.org/), release 9.1.
  Column 1 = NCBI species ID; column 2 = NCBI species name;column 3 = species domain (eukaryotes).
  }
  \item{cogids}{
  A one-column data.frame with orthologous groups IDs available in the 'cogdata' object.
  }
  \item{phyloTree}{
  An object of class "phylo" listing the eukaryotes in the STRING database, release 9.1.
  }
}
}

\value{
a dataset.
}

\references{
Franceschini A, Szklarczyk D, Frankild S, Kuhn M, Simonovic M, Roth A, Lin J, Minguez P, Bork P, von Mering C, Jensen LJ. STRING v9.1: protein-protein interaction networks, with increased coverage and integration. Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29.
}

\examples{
  data(gpdata.gs)
}

\keyword{dataset}

