\name{gplast.preprocess}
\alias{gplast.preprocess}

\title{Evolutionary plasticity inference.}

\description{
Constructor for the 'OGP-class' object.
}

\usage{
gplast.preprocess(cogdata, sspids=NULL, cogids=NULL, verbose=TRUE)
}

\arguments{
\item{cogdata}{
a data frame or matrix with COG's data.
}
\item{sspids}{
an optional data frame with species annotation. Alternatively, it can be a character vector with species IDs.
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
This function creates an \code{\link[geneplast:OGP-class]{OGP-class}} object and checks the consistency of the input data for the evolutionary plasticity pipeline. Internally, the function counts the number of orthologs for each species in a given OG and computes the \code{orthodist} matrix.
}

\value{
A preprocessed object of class 'OGP'.
}

\author{
Rodrigo Dalmolin, Mauro Castro
}


\references{
Dalmolin, RJ, Castro, MA, Rybarczyk-Filho JL, Souza LH, de Almeida RM, Moreira JC. Evolutionary plasticity determination by orthologous groups distribution. Biol Direct. 2011 May 17;6:22. DOI: 10.1186/1745-6150-6-22.

Castro MA, Dalmolin RJ, Moreira JC, Mombach JC, de Almeida RM. Evolutionary origins of human apoptosis and genome-stability gene networks. Nucleic Acids Res. 2008 Nov;36(19):6269-83. DOI: 10.1093/nar/gkn636.
}

\seealso{
  \code{\link[geneplast:OGP-class]{OGP-class}}
}


\examples{

#load datasets used for demonstration
data(gpdata.gs)

#create and object of class 'OGP'
ogp <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogids)

}

\keyword{methods}

