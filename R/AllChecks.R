##------------------------------------------------------------------------
##This function is used for argument checking
geneplast.checks <- function(name, para) {
  if(name=="gplast.get"){
    opts<-c("cogids","sspids","orthodist","abundance",
            "diversity","plasticity","status","results")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", 
                 paste(opts,collapse = ", ") ) ,call.=FALSE)
  } 
  else if(name=="groot.get"){
    opts<-c("cogids","spbranches","orthoroot","status","results","tree")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", 
                 paste(opts,collapse = ", ") ) ,call.=FALSE)
  } 
  else if(name=="groot"){
    if(!class(para)=="groot")
      stop("'object' should be an object of class 'groot'!\n",call.=FALSE)
  } 
  else if(name=="cogdata"){
    if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<3 ){
      stop("'cogdata' object should be a data frame with length >= 3 !\n",
           call.=FALSE)
    }
    clpars<-c("protein_id","ssp_id","cog_id")
    clname<-tolower(colnames(para))
    if(!all(clpars%in%clname)){
      stop("'cogdata' colnames should include: ",paste(clpars,collapse=", "),
           call.=FALSE)
    }
    colnames(para)<-clname
    para<-para[,c(clpars,clname[which(!clname%in%clpars)])]
    para[,1]<-as.character(para[,1])
    para[,2]<-as.character(para[,2])
    para[,3]<-as.character(para[,3])
    for(i in 1:ncol(para)){
      if(!is.numeric(para[,i]) || !is.integer(para[,i])){
        para[,i]<-as.character(para[,i])
      }
    }
    if(is.matrix(para))para<-as.data.frame(para,stringsAsFactors=FALSE,
                                           check.names=FALSE)
    for(i in 1:3){
      idx<-!is.na(para[,i]) & !para[,i]==""
      para<-para[which(idx),]
    }
    if(nrow(para)<=3){
      stop("'cogdata' has no useful data!\n",call.=FALSE)
    }
    return(para)
  } 
  else if(name=="cogids"){
    if(is.null(para)){
      return(para)
    } else {
      b1<-is.character(para)
      b2<-is.matrix(para) || is.data.frame(para)
      if(!b1 && !b2){
        stop("'cogids' should be a vector of data frame! \n",
             call.=FALSE)
      }
      if(b1){
        para<-unique(as.character(para))
        para<-para[!is.na(para)]
        para<-para[para!='']
        para<-sort(para)
        para<-data.frame(cog_id=para,stringsAsFactors=FALSE,check.names=FALSE)
      } else {
        clpars<-c("cog_id")
        clname<-tolower(colnames(para))
        if(!all(clpars%in%clname)){
          stop("'cogids' colnames should include: '",clpars,"'!", call.=FALSE)
        }
        colnames(para)<-clname
        para<-para[,c(clpars,clname[which(!clname%in%clpars)]),drop=FALSE]
        para[,1]<-as.character(para[,1])
        for(i in 1:ncol(para)){
          if(!is.numeric(para[,i]) || !is.integer(para[,i])){
            para[,i]<-as.character(para[,i])
          }
        }
        para<-data.frame(para,stringsAsFactors=FALSE,check.names=FALSE)
        uni<-unique(para[,1])
        uni<-uni[!is.na(uni)]
        uni<-uni[uni!='']
        uni<-sort(uni)
        para<-para[match(uni,para[,1]), ,drop=FALSE]
      }
      rownames(para)<-para[,1]
      return(para)
    }
  } 
  else if(name=="sspids"){
    if(is.null(para)){
      return(para)
    } else {
      b1<-is.character(para) || is.integer(para) || is.numeric(para)
      b2<-is.matrix(para) || is.data.frame(para)
      if(!b1 && !b2){
        stop("'sspids' should be a vector of characters or dataframe! \n",
             call.=FALSE)
      }
      if(b1){
        para<-unique(as.character(para))
        para<-para[!is.na(para)]
        para<-para[para!='']
        para<-sort(para)
        para<-data.frame(ssp_id=para,stringsAsFactors=FALSE)
      } else {
        clpars<-c("ssp_id")
        clname<-tolower(colnames(para))
        if(!all(clpars%in%clname)){
          stop("'cogdata' colnames should include: ",paste(clpars,collapse=", "),
               call.=FALSE)
        }
        colnames(para)<-clname
        para<-para[,c(clpars,clname[which(!clname%in%clpars)]),drop=FALSE]
        para[,1]<-as.character(para[,1])
        para[,2]<-as.character(para[,2])
        for(i in 1:ncol(para)){
          if(!is.numeric(para[,i]) || !is.integer(para[,i])){
            para[,i]<-as.character(para[,i])
          }
        }
        para<-data.frame(para,stringsAsFactors=FALSE)
        uni<-unique(para[,1])
        uni<-uni[!is.na(uni)]
        uni<-uni[uni!='']
        uni<-sort(uni)
        para<-para[match(uni,para[,1]), ,drop=FALSE]
      }
      rownames(para)<-para[,1]
      return(para)
    }
  } 
  else if(name=="spbranches"){
    b1<-is.matrix(para) || is.data.frame(para)
    if(!b1){
      stop("'spbranches' should be a data.frame with spp branches! \n",
           call.=FALSE)
    }
    clpars<-c("ssp_id","ssp_name")
    clname<-tolower(colnames(para))
    clNOTpars<-clname[which(!clname%in%clpars)]
    if(!all(clpars%in%clname)){
      stop("'spbranches' colnames should include: ",paste(clpars,collapse=" AND "),
           call.=FALSE)
    }
    colnames(para)<-clname
    para<-para[,c(clpars,clNOTpars),drop=FALSE]
    para[,1]<-as.character(para[,1])
    para[,2]<-as.character(para[,2])
    for(i in 1:ncol(para)){
      if(!is.numeric(para[,i]) && !is.integer(para[,i])){
        para[,i]<-as.character(para[,i])
      }
    }
    para<-data.frame(para,stringsAsFactors=FALSE,check.names=FALSE)
    #---remove dupplicated ids!
    #uni<-unique(para[,1])
    #uni<-uni[!is.na(uni)]
    #uni<-uni[uni!='']
    #uni<-sort(uni)
    #para<-para[match(uni,para[,1]), ,drop=FALSE]
    #---
    if(any(duplicated(para[,1])))
      stop("NOTE: 'spbranches' should have unique spp ids!")
    if(any(duplicated(para[,2])))
      stop("NOTE: 'spbranches' should have unique spp names!")
    rownames(para)<-para[,1]
    return(para)
  } 
  else if(name=="phyloTree"){
    if(!"phylo"%in%class(para))stop("'phyloTree' should be an object of class 'phylo'!")
  } 
  else if(name=="spid"){
    if(is.integer(para) || is.numeric(para) ) para<-as.character(para)
    if( !is.character(para) || length(para)!=1 )
      stop("'spid' should be a character or integer value!\n",call.=FALSE)
    para<-as.character(para)
    return(para)
  } 
  else if(name=="idkey"){
    if( !is.character(para) || length(para)!=1 )
      stop("'idkey' should be a character value!\n",call.=FALSE)
    para<-as.character(para)
    return(para)
  } 
  else if(name=="cutoff"){
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || 
       para>1 || para<0)
      stop("'cutoff' should be numeric value >=0 and <=1 !\n",call.=FALSE)
  } 
  else if(name=="penalty"){
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<=0)
      stop("'penalty' should be an integer or numeric value > 0 !\n",call.=FALSE)
  } 
  else if(name=="nPermutations") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || 
       para<1 || round(para,0)!=para)
      stop("'nPermutations' should be an integer >=1 !",call.=FALSE)
  } 
  else if(name=="method"){
    opts<-c("KS","BR")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'method' should be one of \n", paste(opts,collapse = ", ") ) )
  } 
  else if(name=="pAdjustMethod") {
    if(!is.character(para) || length(para)!=1 || 
         !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", 
                       "BY", "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm','hochberg',
           'hommel','bonferroni','BH','BY','fdr' and 'none'!",call.=FALSE)
  } 
  else if(name=="lab.offset"){
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1)
      stop("'lab.offset' should be a single numeric value!\n",call.=FALSE)
  } 
  else if(name=="adj.tips"){
    if(!(is.integer(para) || is.numeric(para)) || length(para)>2 || length(para)<1)
      stop("'adj.tips' should be a numeric vector of length 2 !\n",call.=FALSE)
  }
  else if(name=="col.edges"){
    if( !areColors(para) || length(para)!=2)
      stop("'col.edges' should be a color vector of length 2 !\n",call.=FALSE)
  }
  else if(name=="col.tips"){
    if( !areColors(para) || length(para)!=2)
      stop("'col.tips' should be a color vector of length 2 !\n",call.=FALSE)
  }  
  else if(name=="col.root"){
    if( !areColors(para) || length(para)!=1)
      stop("'col.root' should be a single color value!\n",call.=FALSE)
  }    
  else if(name=="plot.subtree") {
    if(!is.logical(para) || length(para)!=1)
      stop("'plot.subtree' should be a single logical value!",call.=FALSE)
  }
  else if(name=="plot.lcas") {
    if(!is.logical(para) || length(para)!=1)
      stop("'plot.lcas' should be a single logical value!",call.=FALSE)
  }
  else if(name=="plot.sspnames") {
    if(!is.logical(para) || length(para)!=1)
      stop("'plot.sspnames' should be a single logical value!",call.=FALSE)
  } 
  else if(name=="fname"){
    if( !is.character(para) || length(para)!=1 )
      stop("'fname' should be a single character value!\n",call.=FALSE)
  }
  else {
    stop("...<",name,"> arg is missing!!",call. = FALSE)
  }
}
areColors <- function(x) {
  all(sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
  }))
}
