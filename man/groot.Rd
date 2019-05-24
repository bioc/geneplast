\name{groot}
\alias{groot}

\title{Evolutionary rooting inference.}

\description{
Function to determine the evolutionary root of a gene based on the distribution of its orthologs.
}

\usage{
groot(object, method="BR", penalty=2, cutoff=0.3, nPermutations=1000, 
pAdjustMethod="bonferroni", verbose=TRUE)
}

\arguments{
\item{object}{
this argument is an object of class 'OGR' (\code{\link[geneplast:OGR-class]{OGR-class}}).
}
\item{method}{
a single character value specifying the rooting algorithm. Options: "BR" and "KS" (see details).
}
\item{penalty}{
a single numeric value specifying the penalty used in the rooting algorithm (see details).
}
\item{cutoff}{
a single numeric value in [0,1] specifying the cutoff used in the BR statistics (see details).
}
\item{nPermutations}{
a single integer value specifying the number of permutations used to compute a null distribution for the inferred roots in the species tree.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be used (see 'p.adjust' for details).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\details{
The 'groot' function addresses the problem of finding the evolutionary root of a feature in an phylogenetic tree. The method infers the probability that such feature was present in the Last Common Ancestor (LCA) of a given lineage. Events like horizontal gene transfer, gene deletion, and de novo gene formation add noise to vertical heritage patterns. The 'groot' function assesses the presence and absence of the orthologs in the extant species of the phylogenetic tree in order to build a probability distribution, which is used to identify vertical heritage patterns.

The 'penalty' argument weighs gene gain and loss; penalty=1 indicates equal probability; penalty > 1 indicates higher probability of gene loss while penalty < 1 indicates higher probability of gene gain (penalty value should be greater than zero; default penalty=2). 

After the probability distribution is built for a given lineage, then a rooting algorithm is used to search the LCA that provides the best vertical heritage pattern (default method='BR'). The rooting algorithms finds the optimum point that splits the probability distribution into two components: one enriched with the queried feature (supporting vertical heritage in the lineage) and another with low evidence in favour of the feature's presence. The cutoff sets the tolerance for the discrimination between the two components (default cutoff=0.3). Based on the optimization settings, then the 'groot' function computes the inconsistency score 'Dscore', which assesses the significance of the inferred root against a null distribution derived by permutation analysis.

}

\value{
An processed object of class 'OGR', including results from the rooting algorithm.
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

## run the groot function
## this example uses the orthologous groups listed in the gpdata object
ogr <- groot(ogr, nPermutations=100)
res <- groot.get(ogr, what="results")

\dontrun{
# Option: parallel version with SNOW package!
library(snow)
options(cluster=makeCluster(2, "SOCK"))
ogr <- groot(ogr, nPermutations=100)
stopCluster(getOption("cluster"))
}

}

\keyword{methods}

