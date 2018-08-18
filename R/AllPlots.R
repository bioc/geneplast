
##################################################
### plot groot
##################################################
groot.plot<-function(ogr, whichOG, fname="gproot",  width=4.5, height=6.5, cex.lab=0.3, 
                     cex.nodes=0.6, adj.tips=c(1, 0.5), lab.offset=1.5, 
                     col.tips=c("green2","grey"),col.edges=c("black","grey"), 
                     col.root="red", plot.sspnames=TRUE, plot.subtree=FALSE, 
                     plot.lcas=FALSE){
  if(!plot.lcas){
    if(ogr@status["Rooting"]!="[x]")
      stop("NOTE: 'ogr' object should be evaluated by the 'groot' function!")
    if(!is.character(whichOG) || length(whichOG)>1)
      stop("'whichOG' should be a single character value!")
  }
  if(class(phyloTree)!="phylo")stop("'phyloTree' should be an object of class 'phylo'!")

  #check args
  if(length(adj.tips)==1)adj.tips<-c(adj.tips,adj.tips) 
  geneplast.checks(name="fname",para=fname)
  geneplast.checks(name="plot.subtree",para=plot.subtree)
  geneplast.checks(name="plot.lcas",para=plot.lcas)
  geneplast.checks(name="plot.sspnames",para=plot.sspnames)
  geneplast.checks(name="lab.offset",para=lab.offset)
  geneplast.checks(name="adj.tips",para=adj.tips)
  geneplast.checks(name="col.tips",para=col.tips)
  geneplast.checks(name="col.edges",para=col.edges)
  geneplast.checks(name="col.root",para=col.root)
  pargs<-list(cex.lab=cex.lab,cex.nodes=cex.nodes,adj.tips=adj.tips, lab.offset=lab.offset,  
              col.tips=col.tips,col.edges=col.edges,col.root=col.root,
              plot.sspnames=plot.sspnames,plot.subtree=plot.subtree)

  #-------------
  # PLOT TREE
  #-------------
  
  if(plot.lcas){
    
    # PLOT LCAS
    
    phyloTree<-ogr@tree
    spbranches<-ogr@spbranches
    refsp<-spbranches$ssp_id[1]
    refspname<-spbranches$ssp_name[1]
    if(plot.sspnames){
      phyloTree$tip.label<-phyloTree$tip.alias
    }
    
    #---get OG root position
    lcas<-getLCAs(phyloTree)
    
    #---plot
    fname<-paste(fname,"_",refsp,"LCAs.pdf",sep="")
    pdf(file=fname, width=width, height=height)
    .plot.lcas(phyloTree, pargs=pargs, lcas=lcas, refsp=refsp, refspname=refspname)
    invisible(dev.off())
    cat("PDF file ",fname, " has been generated!", sep = "'")
    
  } else {
    
    #---get, tree, root of whichOG and other relevant results
    
    phyloTree<-ogr@tree
    orthoroot<-ogr@orthoroot
    if(!whichOG%in%rownames(orthoroot)){
      stop(paste(whichOG,"should be listed in the 'ogr' object!"))
    }
    root<-orthoroot[whichOG,"Root"]
    orthoct<-ogr@orthoct[,whichOG,drop=FALSE]
    lgres<-as.numeric(orthoroot[whichOG,c("Dscore","Pvalue","AdjPvalue")])
    lgres<-c(format(lgres[1], digits=3),format(lgres[2:3], digits=3, scientific=TRUE))
    lgres<-paste(c("Dscore","Pvalue","AdjPvalue"),lgres,sep=" = ")
    spbranches<-ogr@spbranches
    refsp<-spbranches$ssp_id[1]
    
    #---aline orthoct and spbranches with phyloTree
    spbranches<-spbranches[as.character(phyloTree$tip.label),,drop=FALSE]
    orthoct<-orthoct[as.character(phyloTree$tip.label),,drop=FALSE]
    
    if(plot.sspnames){
      rownames(spbranches)<-spbranches$ssp_name
      rownames(orthoct)<-spbranches$ssp_name
      phyloTree$tip.label<-phyloTree$tip.alias
    }
    
    fname<-paste(fname,"_",whichOG,"_",refsp,"LCAs.pdf",sep="")
    
    if(plot.subtree){
      
      # CUT TREE ON THE OG ROOT AND PLOT
      
      warning("please, observe that the OG root is placed in a subtree. 
              In this case, only part of the OG distribution is shown!")
      
      #---filtra tudo acima da raiz!
      subbranch<-spbranches[spbranches[,refsp]<=root,]
      idx<-!phyloTree$tip.label%in%rownames(subbranch)
      subPhyloTree<-drop.tip(phyloTree,tip=phyloTree$tip.label[idx])
      subct<-orthoct[subPhyloTree$tip.label,,drop=FALSE]
      
      #---get OG root position
      rtnode<-subPhyloTree$edge[1,1]
      
      #---plot
      pdf(file=fname, width=width, height=height)
      .plot.phylo(subPhyloTree, whichOG=whichOG, orthoct=subct, rtnode=rtnode, 
                  lgres=lgres, pargs=pargs)
      dev.off()
      cat("PDF file ",fname, " has been generated!", sep = "'")
      
    } else {
      
      # PLOT ALL TREE
      
      #---get OG root position
      lcas<-getLCAs(phyloTree)
      rtnode<-lcas[root-1]
      
      #---plot
      pdf(file=fname, width=width, height=height)
      .plot.phylo(phyloTree, whichOG=whichOG,orthoct=orthoct, rtnode=rtnode, 
                  lgres=lgres, pargs=pargs, root=root, lcas=lcas)
      invisible(dev.off())
      cat("PDF file ",fname, " has been generated!", sep = "'")
    }
    
  }
  
}

##################################################
### Custom (internal) plot.phylo fuction
##################################################
.plot.phylo<-function(x, whichOG, orthoct, rtnode, lgres, pargs, root, lcas){
  
  #---number of segments below the inferred root
  nsegs<- which(x$edge[,2]==rtnode)
  if(length(nsegs)==0)nsegs=0
  
  #---set x,y positions
  # yy <- node.height(x, clado.style = TRUE)
  yy <- .node.height(x)
  xx <- node.depth(x, method = 1)-1
  xx <- max(xx) - xx
  Ntip <- length(x$tip.label)
  Nedge <- dim(x$edge)[1]
  
  #---set plot limits
  y.lim<-c(1, Ntip)
  x.lim <- c(0, NA)
  pin1 <- par("pin")[1]
  strWi <- strwidth(x$tip.label, "inches", cex = pargs$cex.lab)
  xx.tips <- xx[1:Ntip] * 1.04
  alp <- try(uniroot(function(a) max(a * xx.tips + strWi) - pin1, c(0, 1e+06))$root, 
             silent = TRUE)
  if (is.character(alp)){
    tmp <- max(xx.tips) * 1.5
  } else {
    tmp <- max(xx.tips + strWi/alp)
  }
  tmp <- tmp + pargs$lab.offset
  x.lim[2] <- tmp
  
  #---get cols
  col.edges<-pargs$col.edges
  col.root<-pargs$col.root
  col.tips<-rev(pargs$col.tips)
  col.root[2]<-"white"
  col.tips[3]<-colorRampPalette(c(col.tips[1],"white"))(10)[8]
  col.tips[4]<-colorRampPalette(c(col.tips[2],"white"))(10)[8]
  
  #---plot segments
  par(mai = rep(0, 4))
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", 
               axes = FALSE, asp = NA)
  segments(xx[x$edge[, 1]], yy[x$edge[, 1]], xx[x$edge[, 2]], yy[x$edge[, 2]], 
           col=c(rep( col.edges[2],nsegs),rep( col.edges[1],Nedge-nsegs)), 
           lwd=1, lty=1)
  #---plot root
  points(xx[rtnode], yy[rtnode], pch=23, col=col.root[1], 
         bg=col.root[2], cex=pargs$cex.nodes*1.3)
  text(xx[rtnode], yy[rtnode], labels = paste("root",root), cex = 0.5, 
       pos = 2, offset = 0.5)
  
  #---plot tips
  idx<-orthoct[,1]+1
  points(xx[1:Ntip] + pargs$adj.tips[1] - 0.5, yy[1:Ntip] + pargs$adj.tips[2] - 0.5, 
         pch=21, cex=pargs$cex.nodes, lwd=0.5, col=col.tips[1:2][idx],
         bg=col.tips[3:4][idx])
  
  #---plot labels
  text(xx[1:Ntip] + pargs$lab.offset, yy[1:Ntip], x$tip.label, 
       adj = 0, font = 3, srt = 0, cex = pargs$cex.lab, col = "black")
  
  #---plot legend
  legend("topleft",legend=c(whichOG), cex=1.0, bty="n")  
  legend("bottomleft",legend=c("OG presence","OG absence","Inferred root",NA,lgres), 
         col=c(col.tips[2],col.tips[1],col.root[1],NA,NA,NA,NA), 
         pt.bg=c(col.tips[4],col.tips[3],col.root[2],NA,NA,NA,NA),
         pch=c(21,21,23,NA,NA,NA,NA), cex=0.6, pt.cex=0.8, inset = 0.05, bty="n")  
}
# fix a bug (or change) in "node.height" from ape
.node.height <- function(phy){
  n <- length(phy$tip.label)
  m <- phy$Nnode
  N <- dim(phy$edge)[1]
  phy <- reorder(phy)
  yy <- numeric(n + m)
  e2 <- phy$edge[, 2]
  yy[e2[e2 <= n]] <- 1:n
  phy <- reorder(phy, order = "postorder")
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  .C(node_height_clado, as.integer(n), as.integer(e1), 
     as.integer(e2), as.integer(N), double(n + m), as.double(yy))[[6]]
}

##################################################
### Custom (internal) plot.lcas fuction
##################################################
.plot.lcas<-function(x, pargs, lcas, refsp, refspname){

  #---set x,y positions
  # yy <- node.height(x, clado.style = TRUE)
  yy <- .node.height(x)
  xx <- node.depth(x, method = 1)-1
  xx <- max(xx) - xx
  Ntip <- length(x$tip.label)
  Nedge <- dim(x$edge)[1]
  
  #---set plot limits
  y.lim<-c(1, Ntip)
  x.lim <- c(0, NA)
  pin1 <- par("pin")[1]
  strWi <- strwidth(x$tip.label, "inches", cex = pargs$cex.lab)
  xx.tips <- xx[1:Ntip] * 1.04
  alp <- try(uniroot(function(a) max(a * xx.tips + strWi) - pin1, c(0, 1e+06))$root, 
             silent = TRUE)
  if (is.character(alp)){
    tmp <- max(xx.tips) * 1.5
  } else {
    tmp <- max(xx.tips + strWi/alp)
  }
  tmp <- tmp + pargs$lab.offset
  x.lim[2] <- tmp
  
  #---get cols
  col.edges<-pargs$col.edges
  col.root<-pargs$col.root
  col.tips<-rev(pargs$col.tips)
  col.root[2]<-"white"
  col.tips[3]<-colorRampPalette(c(col.tips[1],"white"))(10)[8]
  col.tips[4]<-colorRampPalette(c(col.tips[2],"white"))(10)[8]
  
  #---plot segments
  par(mai = rep(0, 4))
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", ylab = "", 
               axes = FALSE, asp = NA)
  segments(xx[x$edge[, 1]], yy[x$edge[, 1]], xx[x$edge[, 2]], yy[x$edge[, 2]], 
           col=col.edges[2], lwd=1, lty=1)
  
  #---plot root
  text(xx[lcas], yy[lcas], labels = paste("root",2:(length(lcas)+1)), cex = 0.25, 
       pos = 2, offset = 0.1, srt=-45)
  
  #---plot tips
  points(xx[1:Ntip] + pargs$adj.tips[1] - 0.5, yy[1:Ntip] + pargs$adj.tips[2] - 0.5, 
         pch=21, cex=pargs$cex.nodes, lwd=0.5,
         col=col.tips[1],bg=col.tips[3])
  
  #---plot labels
  text(xx[1:Ntip] + pargs$lab.offset, yy[1:Ntip], x$tip.label, 
       adj = 0, font = 3, srt = 0, cex = pargs$cex.lab, col = "black")
  
  #---plot legend
  legend("topleft",legend=paste("REF: ",refspname," (",refsp,")",sep=""), 
         cex=0.8, bty="n")  
  
}
