\name{gplast}
\alias{gplast}

\title{Evolutionary plasticity inference.}

\description{
Function to calculate abundance, diversity, and evolutionary plasticity of an orthologous group (OG).
}

\usage{
gplast(object, verbose=TRUE)
}

\arguments{
\item{object}{
this argument is an object of class 'OGP' (\code{\link[geneplast:OGP-class]{OGP-class}}).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\details{
This method computes the abundance and diversity of an OG, and derives the evolutionary plasticity as described in Castro et al. (2008) and Dalmolin et al. (2011).
The OG diversity corresponds to the normalized Shannon's diversity index and estimates the distribution of orthologous proteins across the species listed in the input dataset. The  OG abundance is simply the ratio between the number of orthologs of a given OG and the number of organisms listed in the group. Evolutionary Plasticity Index is calculated as described in Dalmolin et al (2011).
}

\value{
A processed object of class 'OGP', including COG's abundance, diversity, and plasticity results.
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

## run the gplast function
## this example uses the especies/COGs listed in the gpdata object
ogp <- gplast(ogp)
res <- gplast.get(ogp,what="results")

}

\keyword{methods}

