\name{gpdata}
\alias{gpdata.gs}
\docType{data}

\title{A pre-processed dataset for the geneplast package.}

\description{
A dataset used to demonstrate geneplast main functions.
}

\usage{data(gpdata.gs)}

\format{
  \code{gpdata.gs}{A data frame containing orthology annotation.}
}

\details{

The dataset consists of 4 R objects to be used for demonstration purposes only in the geneplast vignette.

\describe{
  \item{cogdata}{
  A data frame with three columns listing orthology 
  annotation retrieved from the STRING database (http://string-db.org/), 
  release 9.1. Column 1 = Ensembl protein ID; column 2 = NCBI species ID; 
  column 3 = OG ID. Note: This dataset is to be used for demonstration purposes 
  only as it represents a subset of the STRING database; in order to reduce the 
  dataset size, orthology annotation was mapped to genome stability genes (Castro et al.).
  }
  \item{sspids}{
  A data frame containing the species listed in STRING database (http://string-db.org/), release 9.1.
  Column 1 = NCBI species ID; column 2 = NCBI species name;column 3 = species domain (eukaryotes).
  }
  \item{cogids}{
  A one-column data.frame listing orthologous groups (OGs) available in the 'cogdata' object.
  }
  \item{phyloTree}{
  An object of class "phylo" for the eukaryotes listed in the STRING database, release 9.1.
  }
}
}

\value{
a dataset.
}

\references{
Franceschini et al. STRING v9.1: protein-protein interaction networks, with increased coverage and integration. Nucleic Acids Research 41(Database issue):D808-15, 2013. doi:10.1093/nar/gks1094. Epub 2012 Nov 29.

Castro et al. Evolutionary Origins of Human Apoptosis and Genome-Stability Gene Networks. Nucleic Acids Research 36(19): 6269-83, 2008. doi:10.1093/nar/gkn636.
}

\examples{
  data(gpdata.gs)
}

\keyword{dataset}

