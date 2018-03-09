######################################
###  All constructors
######################################
##------------------------------------------------------------------------------
##OGP constructor
gplast.preprocess<-function(cogdata, sspids=NULL, cogids=NULL, verbose=TRUE){
  
  if(verbose)cat("-Preprocessing for input data...\n")
  
  #---main checks
  cogdata=geneplast.checks(name="cogdata",para=cogdata)
  sspids=geneplast.checks(name="sspids",para=sspids)
  cogids=geneplast.checks(name="cogids",para=cogids)
  #get sspids
  if(is.null(sspids)){
    sspids<-unique(as.character(cogdata[,"ssp_id"]))
    sspids<-sspids[!is.na(sspids)]
    sspids<-sspids[sspids!='']
    sspids<-sort(sspids)
    sspids<-data.frame(ssp_id=sspids,stringsAsFactors=FALSE)
    rownames(sspids)<-sspids[,"ssp_id"]
  } else {
    if(any(!sspids[,"ssp_id"]%in%cogdata[,"ssp_id"])){
      stop("NOTE: 'sspids' not listed in 'cogdata'!")
    }
  }
  #get cogids
  if(is.null(cogids)){
    cogids<-unique(as.character(cogdata[,3]))
    cogids<-cogids[!is.na(cogids)]
    cogids<-cogids[cogids!='']
    cogids<-sort(cogids)
    cogids<-data.frame(cogids=cogids,stringsAsFactors=FALSE)
    rownames(cogids)<-cogids[,1]
  } else {
    if(any(!cogids[,1]%in%cogdata$cog_id)){
      stop("NOTE: 'cogids' not listed in 'cogdata'!")
    }
  }
  #remove non-usefull data and compute orthodist
  cogdata<-cogdata[cogdata[,3]%in%cogids[,1],]
  orthodist<-orthoCount(cogdata=cogdata,cogvec=cogids[,1],
                        sspvec=sspids[,1],verbose)  
  object <- new("OGP",cogids=cogids,sspids=sspids,orthodist=orthodist)
  return(object)
  
}

##------------------------------------------------------------------------------
##OGR constructor
groot.preprocess<-function(cogdata, phyloTree, spid, cogids=NULL, verbose=TRUE){
  
  if(verbose)cat("-Preprocessing for input data...\n")
  
  #---main checks
  cogdata<-geneplast.checks(name="cogdata",para=cogdata)
  spid<-geneplast.checks(name="spid",para=spid)
  cogids<-geneplast.checks(name="cogids",para=cogids)
  geneplast.checks(name="phyloTree",para=phyloTree)
  
  #check phyloTree ids in cogdata
  if(any(!phyloTree$tip.label%in%cogdata$ssp_id)){
    stop("NOTE: all id(s) in 'phyloTree' should be listed in 'cogdata'!")
  }
  
  #---rotate phyloTree to set spid at the top
  tip<-which(phyloTree$tip.label==spid)
  if(phyloTree$edge[Nedge(phyloTree),2]!=tip){
    phyloTree<-rotatePhyloTree(phyloTree,spid)
    if(phyloTree$edge[Nedge(phyloTree),2]!=tip){
      warning("NOTE: spid seems not placed at the top of the phyloTree!")
    }
  }
  
  #---compute spbranches from phyloTree
  spbranches<-spBranches(phyloTree=phyloTree,spid=spid)
  spbranches<-geneplast.checks(name="spbranches",para=spbranches)
  
  #get cogids
  if(is.null(cogids)){
    cogids <- cogdata$cog_id[cogdata$ssp_id==spid]
    cogids<-unique(as.character(cogids))
    cogids<-cogids[!is.na(cogids)]
    cogids<-cogids[cogids!='']
    cogids<-data.frame(cog_id=cogids,stringsAsFactors=FALSE, row.names = cogids)
  } else {
    if(any(!cogids$cog_id%in%cogdata$cog_id)){
      stop("NOTE: one or more 'cogids' not listed in 'cogdata'!")
    }
  }
  
  #remove non-usefull data
  cogdata<-cogdata[cogdata$cog_id%in%cogids$cog_id,]
  checkcogid<-unique(cogdata$cog_id[which(cogdata$ssp_id==spid)]) 
  idx<-cogids$cog_id%in%checkcogid
  if(any(!idx)){
    checkcogid<-cogids$cog_id[!idx]
    cogids<-cogids[idx,,drop=FALSE]
    tp<-paste(checkcogid,collapse =",")
    warning("'spid' not listed in one or more 'cogids':\n", tp,call. = FALSE)
  }
  #compute orthoct
  orthoct<-rootCount(cogdata=cogdata,cogvec=cogids$cog_id,sspvec=spbranches$ssp_id, verbose=verbose)
  orthoct<-data.frame(ssp_id=row.names(orthoct),orthoct,stringsAsFactors=FALSE)
  orthoct<-merge(spbranches,orthoct, by="ssp_id", all=TRUE)
  rownames(orthoct)<-orthoct[,1]
  orthoct<-orthoct[,-c(1,2)]
  orthoct<-orthoct[sort.list(orthoct[,spid]),]
  
  #---return OGR
  object <- new("OGR", cogids=cogids, tree=phyloTree, spbranches=spbranches, orthoct=orthoct)
  
  return(object)
}

