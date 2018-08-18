\name{groot}
\alias{groot}

\title{Evolutionary rooting inference.}

\description{
Function to determine the evolutionary root of a gene based on its orthologs distribution.
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
This function addresses the problem of finding the evolutionary root of a feature in an evolutionary tree. The method infers the probability that such feature had been present in each last common ancestor (LCA) of a given existent species in the tree. Regarding genes, this is a delicate task, especially due to mechanisms like horizontal gene transfer, gene deletion, de novo gene formation, among others, which add noise to the vertical heritage pattern. The groot function estimates this information by accessing the presence and absence of orthologs in the tree branches and summarizes the available evidences supporting the ancestor ortholog in the LCA. A probability distribution is derived from this process, which sets the searching space used in the rooting algorithm described next.
The penalty argument allows weighing gene gain and loss during speciation; penalty=1 indicates same probability; penalty > 1 indicates higher probability of gene loss while penalty < 1 indicates higher probability of gene gain. Penalty value should be greater than zero (default penalty=2).
The rooting algorithm is set by the argument method. The current version introduces the Bridge (BR) algorithm, which search the previously computed probability distribution to find the most reliable evolutionary root. The algorithm search for the optimum point that splits the probability distribution in two components: one enriched with the queried feature and another with low evidence in favor of the feature's presence. The cutoff sets the tolerance for the discrimination (default cutoff=0.3). The groot function also computes the inconsistence score Dscore in order to accesses the significance of the observed root, which is tested against a null distribution derived by permutation analysis.
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

