
######################################################################
### ogr2igraph
######################################################################
ogr2igraph <- function(ogr, cogdata, g, idkey="ENTREZ"){
  
  #--- check input objects
  if(class(ogr)!="OGR")stop("NOTE: input 'ogr' should be a 'OGR' class object! ")
  if(class(g)!="igraph")stop("NOTE: input 'g' should be an 'igraph' class object! ")
  cogdata=geneplast.checks(name="cogdata",para=cogdata)
  idkey=geneplast.checks(name="idkey",para=idkey)
  
  #--- check status
  if(ogr@status["Rooting"]!="[x]")
    stop("NOTE: input 'ogr' needs 'groot' evaluation!")
  
  #--- check idkey
  if(!idkey%in%vertex_attr_names(g)){
    stop("NOTE: 'idkey' should be listed as vertex attribute in 'g'!")
  }
  vids <- vertex_attr(g,idkey)
  
  #--- get ssp_id
  spbranch <- groot.get(ogr, what="spbranches")
  ssp_id <- names(spbranch)[3]
  
  #--- get cogSP
  cogSP <- cogdata[cogdata$ssp_id==ssp_id,]
  results <- groot.get(ogr, what="results")
  idx <- match(cogSP$cog_id, rownames(results))
  cogSP$Root <- results$Root[idx]
  
  #--- find matching col
  mcol <-sapply(1:ncol(cogSP), function(i){
    sum(vids%in%cogSP[[i]])
  })
  mcol <- which.max(mcol)
  
  #--- count matches
  nm <- round( (sum(vids%in%cogSP[[mcol]])/length(vids) ) * 100, digits = 0)
  if(nm < 50){
    stop("NOTE: 'idkey' maps only ",nm,"% of the IDs between 'cogdata' and 'g' objects!")
  } else if(nm < 90){
    warning("NOTE: 'idkey' maps ",nm,"% of the IDs between 'cogdata' and 'g' objects!")
  }
  
  #--- transfer rooting information
  idx <- match(vids,cogSP[[mcol]])
  g <- set_vertex_attr(g, name="COGID", value=cogSP$cog_id[idx])
  g <- set_vertex_attr(g, name="Root", value=cogSP$Root[idx])
  return(g)
  
}

######################################################################
### ogr2tni
######################################################################
ogr2tni <- function(ogr, cogdata, tni){
  
  #---
  if(class(ogr)!="OGR")stop("NOTE: input 'ogr' should be a 'OGR' class object! ")
  if(class(tni)!="TNI")stop("NOTE: input 'tni' should be a 'TNI' class object! ")
    
  #---
  if(ogr@status["Rooting"]!="[x]")
    stop("NOTE: input 'ogr' needs 'groot' evaluation!")
  if(tni@status["DPI.filter"]!="[x]")
    stop("NOTE: input 'tni' needs 'tni.dpi.filter' evaluation!")
  
  #--- check agreement between ogr, cogdata, tni!
  
  
  
  
  
  
  #--- get ssp_id
  spbranch <- groot.get(ogr, what="spbranches")
  ssp_id <- names(spbranch)[3]
  
  #---- get cogSP
  cogSP <- cogdata[cogdata$ssp_id==ssp_id,]
  results <- groot.get(ogr, what="results")
  idx <- match(cogSP$cog_id, rownames(results))
  cogSP$Root <- results$Root[idx]
  
  #---- match geneannot and rooting information
  geneannot <- tni@rowAnnotation
  rownames(geneannot) <- geneannot$ENTREZ
  idx <- match(geneannot$ENTREZ,cogSP$gene_id)
  geneannot$COGID <- cogSP$cog_id[idx]
  geneannot$Root <- cogSP$Root[idx]
  
  #--- update tni with rooting information
  tni@rowAnnotation <- cbind(tni@rowAnnotation,geneannot[,c("COGID","Root")])
  validTar <- tni@rowAnnotation$PROBEID[!is.na(tni@rowAnnotation$Root)]
  validReg <- tni@regulatoryElements[tni@regulatoryElements %in% validTar]
  tni@regulatoryElements <- validReg
  tni@rowAnnotation <- tni@rowAnnotation[validTar,]
  tni@gexp <- tni@gexp[validTar,]
  tni@results$tn.ref <- tni@results$tn.ref[validTar,validReg]
  tni@results$tn.dpi <- tni@results$tn.dpi[validTar,validReg]
  
  return(tni)
  
}


######################################################################
### EPI supplements
######################################################################
#count protein number for each cog
orthoCount<-function(cogdata,cogvec,sspvec,verbose){
  DT <- as.data.table(cogdata) 
  len<-length(cogvec)
  cog_id <-NULL
  if(verbose) pb <- txtProgressBar(style=3)
  orthotable<-sapply(1:len,function(i){
    if(verbose) setTxtProgressBar(pb, i/len)
    dt<-DT[cog_id==cogvec[i]]
    sapply(sspvec,function(sp){
      sum(dt$ssp_id==sp)
    })
  })
  if(verbose) close(pb)
  rownames(orthotable)<-sspvec
  colnames(orthotable)<-cogvec
  return(orthotable)
}
#get abundance for each cog
cogAbundance<-function(orthodist){
  res<-rep(0,ncol(orthodist))
  for(i in 1:ncol(orthodist)){
    tp<-orthodist[,i]
    if(sum(tp)>0){
      tp<-tp[tp>0]
      res[i]<-mean(tp)
    }
  }
  names(res)<-colnames(orthodist)
  return(res)
}
#get diversity of each cog
cogDiversity<-function(orthodist){
  #get probs
  sspProb<-function(orthodist){
    res<-sapply(1:ncol(orthodist),function(i){
      tp<-orthodist[,i]
      sm<-sum(tp)
      tp/ifelse(sm==0,1,sm)
    })
    is.nan(res)
    rownames(res)<-rownames(orthodist)
    colnames(res)<-colnames(orthodist)
    return(res)
  }
  probs<-sspProb(orthodist)
  plnp<-probs*log(probs)
  plnp[is.nan(plnp)]<-0
  div<--colSums(plnp)/log(nrow(plnp))
  return(div)
}
#compute epi of each cog
epidx<-function(div, abd){
  res<-1-(div/sqrt(abd))
  res[is.nan(res)]=0
  return(res)
}

######################################################################
###  KS statistic
######################################################################
ksmethod<-function(rootprobs,orthoct){
  #compute roots
  Root<-NULL
  Dscore<-NULL
  D<-NULL
  nroot<-nrow(rootprobs)
  for(i in 1:ncol(rootprobs)){
    proot<-c(rootprobs[,i],0)
    res<-sapply(1:nroot,function(j){
      x<-proot[1:j]
      y<-proot[(j+1):(nroot+1)]
      dr<-0
      if(x[j]>y[1]){
        dr<-droot(x,y,FALSE)
      }
      dr
    })
    rt<-which(res==max(res))
    if(length(rt)>1){
      dscore<-sapply(1:length(rt),function(l){
        mean(proot[1:rt[l]])-mean(proot[(rt[l]+1):(nroot+1)])
      })
      idx<-which(dscore==max(dscore))[1]
      rt<-rt[idx]
      dscore<-dscore[idx]
    } else {
      dscore<-mean(proot[1:rt])-mean(proot[(rt+1):(nroot+1)])
    }
    Root<-c(Root,rt)
    D<-c(D,res[rt])
    Dscore<-c(Dscore,dscore)
  }
  orthoroot<-data.frame(Root=Root,D=D,Dscore=Dscore,stringsAsFactors=FALSE)
  rownames(orthoroot)<-colnames(rootprobs)
  return(orthoroot)
}
# KS statistic
droot<-function (x, y, doplot=FALSE) {
  x=1-x;y=1-y
  x <- x[!is.na(x)]
  n <- length(x)
  y <- y[!is.na(y)]
  n.x <- as.double(n)
  n.y <- length(y)
  n <- n.x * n.y/(n.x + n.y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  if (length(unique(w)) < (n.x + n.y)) {
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
  }
  if(doplot){
    fx <- ecdf(x)
    fy <- ecdf(y)
    xx <- seq(min(x, y), max(x, y), length.out=length(x))
    x0 <- mean(xx[which( fx(xx)-fy(xx) == max(fx(xx)-fy(xx)) )])
    y0 <- fx(x0)
    y1 <- fy(x0)
    plot(fx, xlim=range(c(x, y)), lty=1, do.points=FALSE, verticals=TRUE,main="",
         ylab="Cumulative probability",col="blue")
    plot(fy, add=TRUE, do.points=FALSE, col="black", verticals=TRUE)
    points(c(x0, x0), c(y0, y1), pch=16, col="red")
    segments(x0, y0, x0, y1, col="red", lty=1)
  }
  max(z)
}

######################################################################
###  Bridge statistic
######################################################################
brmethod<-function(rootprobs,cutoff){
  #compute roots
  Root<-NULL
  Dscore<-NULL
  D<-NULL
  nroot<-nrow(rootprobs)+1
  for(i in 1:ncol(rootprobs)){
    proot<-c(rootprobs[,i],0)
    rt<-NA
    dr<-1
    while(is.na(rt)){
      k<-which(proot<cutoff)[1]
      res<-sapply( (k+1):nroot, function(j){
        mean(proot[k:j])
      })
      bridge<-which(res>cutoff)
      if(length(bridge)>0){
        dr<-res[bridge[1]]
        proot[1:(k+bridge[1])]<-1
      } else {
        rt<-k-1
      }
    }
    proot<-c(rootprobs[,i],0)
    dscore<-mean(proot[1:rt])-mean(proot[(rt+1):(nroot)])
    Root<-c(Root,rt)
    D<-c(D,dr)
    Dscore<-c(Dscore,dscore)
  }
  orthoroot<-data.frame(Root=Root, D=D, Dscore=Dscore, stringsAsFactors=FALSE)
  rownames(orthoroot)<-colnames(rootprobs)
  return(orthoroot)
}

######################################################################
###  Count protein number for each cog
######################################################################
rootCount<-function(cogdata,cogvec,sspvec, verbose){
  if(verbose) pb <- txtProgressBar(style=3)
  len <- length(cogvec)
  orthotable<-sapply(1:len,function(i){
    if(verbose) setTxtProgressBar(pb, i/len)
    dt<-cogdata[which(cogdata$cog_id==cogvec[i]),]
    sapply(sspvec,function(sp){sum(dt$ssp_id==sp)})
  })
  if(verbose) close(pb)
  rownames(orthotable)<-sspvec
  colnames(orthotable)<-cogvec
  orthotable[orthotable>0]=1
  return(orthotable)
}

######################################################################
###  Permutation analysis
######################################################################
runPermutation<-function(orthoroot, rootprobs, nPermutations, verbose){
  #function to compute Dscore from random rootprobs
  computeNullRootprobs<-function(rootprobs,orthoroot){
    sapply(1:ncol(rootprobs),function(i){
      gedAdjDscore(proot=rootprobs[,i],rt=orthoroot[i,"Root"],TRUE)
    })
  }
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("computeNullRootprobs","rootprobs",
                         "orthoroot","gedAdjDscore"),envir=environment())
    nulldist<-parSapply(cl, 1:nPermutations, function(i) {
      computeNullRootprobs(rootprobs,orthoroot)
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    nulldist<-sapply(1:nPermutations,function(i){
      if(verbose) setTxtProgressBar(pb, i/nPermutations)
      computeNullRootprobs(rootprobs,orthoroot)
    })
    if(verbose) close(pb)
  }
  adjDscore<-sapply(1:ncol(rootprobs),function(i){
    gedAdjDscore(rootprobs[,i],rt=orthoroot[i,"Root"])
  })
  #z-transformation
  xmd<-apply(nulldist,1,mean)
  xsd<-apply(nulldist,1,sd)
  zscore<-(adjDscore-xmd)/xsd
  pvals<-pnorm(zscore,lower.tail=FALSE)
  return(list(zscore=zscore,pvals=pvals))
}
#---
gedAdjDscore<-function(proot,rt,perm=FALSE){
  nroot<-length(proot)
  #pspace<-c(rep(1, nroot),rep(0, nroot))
  pspace<- c(1:(nroot*2))%%2
  ct<-1+nroot-rt
  mask<-ct:(ct+nroot-1)
  if(perm){
    pspace[mask]<-sample(pspace)[mask]
  } else {
    pspace[mask]<-proot
  }
  rt<-nroot
  nroot<-length(pspace)
  mean(pspace[1:rt])-mean(pspace[(rt+1):(nroot)])
}
# proot<-c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0);rt=8
# proot<-c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1);rt=8
# gedAdjDscore(proot,rt)
# null<-sapply(1:1000,function(i){
#   gedAdjDscore(proot=proot,rt=rt,TRUE)
# })
# max(null)
# plot(density(null))
#---
# testdnoise<-function(nroot=25,nPermutations=1000,noise=0){
#   sapply(1:(nroot),function(rt){
#     proot<-rep(0,nroot);proot[1:rt]=1
#     if(noise>0)proot[sample(nroot,noise)]<-sample(c(0,1),noise,replace=TRUE)
#     gedAdjDscore(proot,rt)-gedAdjDscore(proot,rt,TRUE)
#   })
# }
# res<-testdnoise(nroot=25,nPermutations=1000,noise=25);plot(density(res))
#---
#check nulldist across rootprobs for the best distributions
# testddist<-function(nroot=25,nPermutations=1000){
#   testsd<-function(proot,rt,nPermutations){
#     xd<-sapply(1:nPermutations,function(i){gedAdjDscore(proot,rt,TRUE)})
#     xd<-quantile(xd,probs=1-1/nPermutations,names=FALSE)
#     gedAdjDscore(proot,rt)/xd
#   }
#   sapply(1:(nroot),function(rt){
#     proot<-rep(0,nroot);proot[1:rt]=1
#     testsd(proot,rt,nPermutations)
#   })
# }
# testddist(nroot=25,nPermutations=1000)
#---

######################################################################
###  pars-based probs
######################################################################
getprobs<-function(x,penalty){
  loss<-length(x)-sum(x)
  gain<-sum(x)*penalty
  gain/(gain+loss)
}

######################################################################
###  check if snow cluster is loaded
######################################################################
isParallel<-function(){
  b1<-"package:snow" %in% search()
  b2<-tryCatch({
    cl<-getOption("cluster")
    cl.check<-FALSE
    if(is(cl, "cluster")){
      cl.check <- all( sapply(1:length(cl),function(i){
        isOpen(cl[[i]]$con) }) == TRUE )
    }
    cl.check
  }, error=function(e){ FALSE 
  })
  all(c(b1,b2))
}


######################################################################
### Return LCAs and branches from a 'phylo' obj
######################################################################
spBranches<-function(phyloTree,spid){
  if(all(phyloTree$tip.label!=spid)){
    stop("NOTE: 'spid' should be listed in 'phyloTree'!")
  }
  lcas<-getLCAs(phyloTree)
  spbranches<-getBranches(phyloTree,lcas)
  #---set spid as the 1st branch
  spbranches$branch<-spbranches$branch+1
  spbranches[spid,"branch"]<-1
  #---set obj to 'sspbranches' format
  spbranches<-spbranches[,c("ssp_id","ssp_name","branch")]
  spbranches<-spbranches[sort.list(spbranches$branch),]
  colnames(spbranches)[3]<-spid
  rownames(spbranches)<-spbranches$ssp_id
  return(spbranches)
}
getLCAs<-function(phyloTree){
  ntips<-length(phyloTree$tip.label)
  edgetree<-phyloTree$edge
  tip<-edgetree[nrow(edgetree),2]
  LCAs<-edgetree[which(edgetree[,2]==tip),1]
  while( LCAs[length(LCAs)]>(ntips+1) ){
    idx<-which(edgetree[,2] == LCAs[length(LCAs)] )
    res<-edgetree[idx,1]
    LCAs<-c(LCAs,res)
  }
  return(LCAs)
}
getBranches<-function(phyloTree, LCAs){
  ntips<-length(phyloTree$tip.label)
  edgetree<-phyloTree$edge
  edges<-1:nrow(edgetree)
  lcas<-rev(which(edgetree[,2]%in%LCAs))
  edgetree<-cbind(edgetree,NA)
  for(loc in LCAs){
    idx1<-which(edgetree[,2]==loc)
    idx2<-which(edges>idx1)
    edgetree[edges[idx2],3]<-loc
    edges<-edges[-idx2]
  }
  edgetree[is.na(edgetree[,3]),3]<-edgetree[1,1]
  branches<-edgetree[edgetree[,2]<=ntips,2:3]
  #---
  ids<-phyloTree$tip.label[branches[,1]]
  if(!is.null(phyloTree$tip.alias)){
    alias<-phyloTree$tip.alias[branches[,1]]
  } else {
    alias<-ids
  }
  branches<-data.frame(ids,alias,branches,stringsAsFactors=FALSE)
  colnames(branches)<-c("ssp_id","ssp_name","tip","lca")
  gps<-as.integer(as.factor(rank(-branches$lca)))
  branches$branch<-gps
  rownames(branches)<-branches$ssp_id
  return(branches)
}
tipOrder<-function(phyloTree){
  tporder<-phyloTree$edge[,2]
  tporder<-tporder[tporder<=Ntip(phyloTree)]
  tporder<-as.character(phyloTree$tip.label[tporder])
  return(tporder)
}
rotatePhyloTree<-function(phyloTree,spid){
  tip <- which(phyloTree$tip.label==spid)
  lcas <- mrca(phyloTree)[,spid]
  phyloTree$edge.length <- rep(1, nrow(phyloTree$edge))
  tgroup<-dist.nodes(phyloTree)[,tip]
  tgroup<-tgroup[lcas]
  names(tgroup)<-names(lcas)
  #---
  ct<-1;tp<-tgroup
  for(i in sort(unique(tgroup))){
    tgroup[tp==i]<-ct;ct<-ct+1
  }
  #---
  tord<-rev(rank(rev(tgroup), ties.method = "first"))
  #---
  phyloTree<-rotateConstr(phyloTree,names(sort(tord,decreasing=TRUE)))
  tord<-tord[tipOrder(phyloTree)]
  #---
  tgroup<-tgroup[names(tord)]
  phyloTree$tip.group<-tgroup
  #---
  lcas<-lcas[names(tord)]
  #atualiza lca do spid, troca pelo nodo mais proximo
  tp<-phyloTree$edge[,2]
  lcas[spid]<-phyloTree$edge[which(tp==tip),1]
  phyloTree$tip.lcas<-lcas
  #---
  return(phyloTree)
}
# getLCAs<-function(phyloTree){
#   edgetree<-phyloTree$edge
#   LCAs<-edgetree[1,1]
#   while( !is.na(LCAs[length(LCAs)]) ){
#     idx<-which(edgetree[,1] == LCAs[length(LCAs)] )[2]
#     res<-edgetree[idx,2]
#     LCAs<-c(LCAs,res)
#   }
#   LCAs<-LCAs[1:(length(LCAs)-1)]
#   return(LCAs)
# }

