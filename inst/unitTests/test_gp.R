# Unit tests fot OGP-class methods
test_ogp <- function(){
  data(gpdata.gs)
  ogp <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogids)
  ogp <- gplast(ogp, verbose=FALSE)
  res <- gplast.get(ogp, what="results")
  checkTrue(is.data.frame(res) && ncol(res)==3)
}
# Unit tests fot OGR-class methods
test_ogr <- function(){
  data(gpdata.gs)
  ogr <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogids)
  ogr <- groot(ogr, nPermutations=50, verbose=FALSE)
  res <- groot.get(ogr,what="results")
  checkTrue(is.data.frame(res) && ncol(res)==4)
}
