######################################
###  All methods
######################################
##------------------------------------------------------------------------------
##initialization
setMethod("initialize","OGP",
          function(.Object, cogids, sspids, orthodist){
            ##-----initialization
            .Object@cogids=cogids
            .Object@sspids=sspids
            .Object@orthodist=orthodist
            .Object@status=c(Preprocess="[x]", Plasticity="[ ]")
            .Object
          }
)
##initialization
setMethod("initialize","OGR",
          function(.Object, cogids, tree, spbranches, orthoct){
            ##-----initialization
            .Object@cogids=cogids
            .Object@tree=tree
            .Object@spbranches=spbranches
            .Object@orthoct=orthoct
            .Object@status=c(Preprocess="[x]", Rooting="[ ]")
            .Object
          }
)

##------------------------------------------------------------------------------
##gplast
setMethod("gplast","OGP",
          function(object, verbose=TRUE){
            
            if(verbose)cat("-Performing plasticity analysis...\n")
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: input 'object' needs preprocessing!")
            
            #---run EPI
            abundance<-cogAbundance(object@orthodist)
            diversity<-cogDiversity(object@orthodist)
            plasticity<-epidx(div=diversity, abd=abundance)
            
            #---return OGP
            object@abundance=abundance
            object@diversity=diversity
            object@plasticity=plasticity
            object@status["Plasticity"]="[x]"
            return(object)
            
          }
)
##------------------------------------------------------------------------------
##groot
setMethod("groot","OGR",
          function(object, method="BR", penalty=2, cutoff=0.3, nPermutations=1000, 
                   pAdjustMethod="bonferroni", verbose=TRUE){
            
            if(object@status["Preprocess"]!="[x]")
              stop("NOTE: input 'object' needs preprocessing!")
            
            if(isParallel()){
              if(verbose)
                cat("-Performing rooting analysis (parallel version - ProgressBar not available)...\n")
            } else {
              if(verbose)cat("-Performing rooting analysis...\n")
            }
            
            #---main checks
            geneplast.checks(name="method",para=method)
            geneplast.checks(name="penalty",para=penalty)
            geneplast.checks(name="cutoff",para=cutoff)
            geneplast.checks(name="nPermutations",para=nPermutations)
            geneplast.checks(name="pAdjustMethod",para=pAdjustMethod)
            
            if(verbose)cat("--For", nrow(object@cogids), "orthologous groups...\n")
            
            #run KS/BR
            if(method=="KS")penalty=1
            rootprobs<-aggregate(object@orthoct[,-1],by=list(object@orthoct[,1]),
                                 FUN=getprobs, penalty=penalty)
            rownames(rootprobs)<-rootprobs[,1]
            rootprobs<-rootprobs[,-1]
            
            #---testing chunk!
            #rootprobs[,]<-sample(c(0,1),prod(dim(rootprobs)),replace=TRUE)
            # #nr<-nrow(rootprobs)
            # for(i in 1:ncol(rootprobs)){
            #   rootprobs[,i]<-as.numeric( rnorm(nr,sd=(1/(1:nr))) > (1/nr) )
            # }
            #plot(1:nr,rowSums(rootprobs))
            #rootprobs[1,]<-1

            #---
            if(method=="KS"){
              orthoroot<-ksmethod(rootprobs,object@orthoct)
            } else {
              orthoroot<-brmethod(rootprobs,cutoff)
            }
            sspCount<-colSums(object@orthoct[,-1])
            
            #---stats
            stats<-runPermutation(orthoroot,rootprobs,nPermutations,verbose)
            stats$adjpvals<-p.adjust(stats$pvals,method=pAdjustMethod)
            
            #---update orthoroot
            orthoroot<-cbind(orthoroot,Zscore=stats$zscore,Pvalue=stats$pvals,
                             AdjPvalue=stats$adjpvals,sspCount=sspCount)
            orthoroot<-orthoroot[sort.list(orthoroot$Root,decreasing=FALSE),]
            orthorootsort<-NULL
            #---
            for(i in unique(orthoroot$Root)){
              tp<-orthoroot[orthoroot$Root==i,]
              tp<-tp[sort.list(tp$Pvalue),]
              orthorootsort<-rbind(if(!is.null(orthorootsort))orthorootsort,tp)
            }
            object@orthoroot=orthorootsort
            object@status["Rooting"]="[x]"
            return(object)
          }
)
##------------------------------------------------------------------------------
##get slots from gplast 
setMethod(
  "gplast.get","OGP",
  function(object, what="status"){
    ##-----check input arguments
    geneplast.checks(name="gplast.get",para=what)
    ##-----get query
    query<-NULL  
    if(what=="cogids"){
      query<-object@cogids
    } else if(what=="sspids"){
      query<-object@sspids
    } else if(what=="orthodist"){
      query<-object@orthodist
    } else if(what=="abundance"){
      query<-object@abundance
    } else if(what=="diversity"){
      query<-object@diversity
    } else if(what=="plasticity"){
      query<-object@plasticity
    } else if(what=="status"){
      query<-object@status
    } else if(what=="results"){
      if(object@status["Plasticity"]!="[x]"){
        warning("NOTE: input 'object' needs 'gplast' evaluation!")
        query<-data.frame()
      } else {
        nms<-names(object@abundance)
        query<-data.frame(
          abundance=round(object@abundance[nms],4),
          diversity=round(object@diversity[nms],4),
          plasticity=round(object@plasticity[nms],4),
          stringsAsFactors=FALSE
        )
      }
    }
    return(query)
  }
)
##------------------------------------------------------------------------------
##get slots from groot 
setMethod(
  "groot.get","OGR",
  function(object, what="status"){
    ##-----check input arguments
    geneplast.checks(name="groot.get",para=what)
    ##-----get query
    query<-NULL  
    if(what=="cogids"){
      query<-object@cogids
    } 
    else if(what=="spbranches"){
      query<-object@spbranches
    }
    else if(what=="tree"){
      query<-object@tree
    } 
    else if(what=="orthoroot"){
      query<-object@orthoroot
    } 
    else if(what=="status"){
      query<-object@status
    } 
    else if(what=="results"){
      if(object@status["Rooting"]!="[x]"){
        warning("NOTE: input 'object' needs 'groot' evaluation!")
        query<-data.frame()
      } else {
        query<-data.frame(
          Root=object@orthoroot$Root,
          Dscore=round(object@orthoroot$Dscore,2),
          Pvalue=signif(object@orthoroot$Pvalue,3),
          AdjPvalue=signif(object@orthoroot$AdjPvalue,3),
          stringsAsFactors=FALSE
        )
        rownames(query)<-rownames(object@orthoroot)
      }
    }
    return(query)
  }
)
##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show","OGP",
  function(object) {
    cat("An OGP (Orthologous Gene Plasticity) object:\n")
    message("--status:")
    print(gplast.get(object, what="status"), quote=FALSE)
  }
)
setMethod(
  "show","OGR",
  function(object) {
    cat("An OGR (Orthologous Gene Root) object:\n")
    message("--status:")
    print(groot.get(object, what="status"), quote=FALSE)
  }
)

