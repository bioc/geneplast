\name{groot.plot}
\alias{groot.plot}

\title{Plot the inferred evolutionary root of a given OG or the map of LCAs of a given species.}

\description{
Plot the inferred evolutionary root of a given OG onto the species tree or the map of LCAs of a given species.
}

\usage{
groot.plot(ogr, whichOG, fname="gproot",  width=4.5, height=6.5, cex.lab=0.3, 
  cex.nodes=0.6, adj.tips=c(1, 0.5), lab.offset=1.5, col.tips=c("green2","grey"), 
  col.edges=c("black","grey"), col.root="red", plot.sspnames=TRUE, 
  plot.subtree=FALSE, plot.lcas=FALSE)
}

\arguments{
\item{ogr}{
this argument is an object of class 'OGR' evaluated by the groot \code{\link[geneplast:groot]{groot}} method.
}
\item{whichOG}{
a single character value indicating the OG to be plotted.
}
\item{fname}{
a character string naming a file.
}
\item{width}{
a single numeric value specifying the width of the graphics region in inches.
}
\item{height}{
a single numeric value specifying the height of the graphics region in inches.
}
\item{cex.lab}{
numeric character expansion factor for tip labels.
}
\item{cex.nodes}{
numeric expansion factor for node symbols.
}
\item{adj.tips}{
two numeric values specifying the adjustment of the labels.
}
\item{lab.offset}{
a single numeric value specifying the offset of the labels.
}
\item{col.tips}{
a character vector of length=2 specifying the colors of the tips.
}
\item{col.edges}{
a character vector of length=2 specifying the colors of the edges.
}
\item{col.root}{
a character value specifying the color of the inferred root.
}
\item{plot.sspnames}{
a single logical value specifying whether ssp names should be used to generate the plot.
}
\item{plot.subtree}{
a single logical value specifying whether a sub-species tree should be used to generate the plot.
}
\item{plot.lcas}{
a single logical value specifying whether a species tree should be generated mapping the positions
of all possible roots.
}
}

\value{
a pdf file.
}

\author{
Rodrigo Dalmolin, Mauro Castro
}


\references{
Dalmolin RJ and Castro, MA. Geneplast: Evolutionary rooting using orthologous groups distribution. Journal Paper (in preparation), 2016.
}

\seealso{
  \code{\link[geneplast:groot]{groot}}
}


\examples{

#load datasets used for demonstration
data(gpdata.gs)

#create and object of class 'OGR' for H. sapiens
ogr <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogids)

## run the groot function
ogr <- groot(ogr, nPermutations=100)

## this example plots NOG40170 in the phyloTree
groot.plot(ogr,whichOG="NOG40170")

}

\keyword{methods}

