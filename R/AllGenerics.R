##generic functions
setGeneric("gplast", function(object, verbose=TRUE) standardGeneric("gplast"), 
           package="geneplast")
setGeneric("gplast.get", function(object, what="status") standardGeneric("gplast.get"), 
           package="geneplast")
setGeneric("groot", function(object, method="BR", penalty=2, cutoff=0.3, 
                             nPermutations=1000, pAdjustMethod="bonferroni", 
                             verbose=TRUE) standardGeneric("groot"), package="geneplast")
setGeneric("groot.get", function(object, what="status") standardGeneric("groot.get"), 
           package="geneplast")
