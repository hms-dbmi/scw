require(scde)
require(Cairo)

t.postprocess <- function(res,name,port,env=go.env,colcols=NULL,n.clusters=10,distance.threshold=0.9, z.score=3,perplexity=20,seed=0,include.aspects=TRUE,top.aspects=50,return.details=F,browse=F) {
  cat("processing ", name,": ")
  pwpca <- res$pwpca; clpca <- res$clpca; cd <- res$cd; prior <- res$prior; knn <- res$knn; varinfo <- res$varinfo;
  
  tam <- pagoda.top.aspects(pwpca[!grepl("custom",names(pwpca))],clpca,use.oe.scale=F,z.score=z.score,adjust.scores=T)
  grcol <- colorRampPalette(c("white","black"),space="Lab")(1024);
  libsize <- colSums(cd)
  scol <- grcol[round((libsize-min(libsize))/diff(range(libsize))*(length(grcol)-1))+1]
  names(scol) <- colnames(cd)
  scol <- scol[rownames(knn)];
  
  if(!is.null(varinfo$batch)) {
    bcol<-rainbow(length(levels(varinfo$batch)),v=0.5)[as.integer(varinfo$batch)] 
    names(bcol) <- names(varinfo$batch)
  } else {
    bcol<-NULL;
  }
  
  
  ## x <- pagoda.cluster.cells(tam,varinfo,include.aspects=include.aspects,verbose=TRUE,return.details=T)
  ## hc <- x$clustering;
  ## cl <- cutree(hc,n.clusters)
  ## clcol <- rainbow(length(unique(cl)))
  ## clcol <- clcol[cl]
  ## require(Cairo)
  ## CairoPNG(file=paste(name,"tam.png",sep="."),width=600,height=600)
  ## tamr <- pagoda.reduce.loading.redundancy(tam,pwpca,clpca,plot=T,labRow=NA,labCol=NA,box=T,margins=c(0.5,0.5),distance.threshold=0.001,cell.clustering=hc,n.cores=1)
  ## dev.off();

  ## CairoPNG(file=paste(name,"tamr.png",sep="."),width=600,height=600)
  ## tamr2 <- pagoda.reduce.redundancy(tamr,distance.threshold=distance.threshold,plot=T,cell.clustering=hc,labRow=NA,labCol=NA,box=T,margins=c(0.5,0.5),top=top.aspects,trim=0,weighted.correlation=T,col.cols=rbind(scol))
  ## dev.off();

  require(Cairo)
  tamr <- pagoda.reduce.loading.redundancy(tam,pwpca,clpca,plot=F,labRow=NA,labCol=NA,box=T,margins=c(0.5,0.5),distance.threshold=0.001,n.cores=1)
  tamr2 <- pagoda.reduce.redundancy(tamr,distance.threshold=distance.threshold,plot=F,labRow=NA,labCol=NA,box=T,margins=c(0.5,0.5),top=top.aspects,trim=0,weighted.correlation=T,col.cols=rbind(scol))
  x <- pagoda.cluster.cells(tamr,varinfo,include.aspects=TRUE,verbose=TRUE,return.details=T)
  hc <- x$clustering;
  cl <- cutree(hc,n.clusters)
  clcol <- rainbow(length(unique(cl)))
  clcol <- clcol[cl]
  
  if(!is.null(colcols)) { # side colors were supplied
    # enforce column order
    if(is.list(colcols)) {
      # new format
      colcols <- lapply(colcols,function(d) {
        so <- match(colnames(varinfo$mat),d$data)
        d$data <- d$data[so]
        if(!is.null(d$text)) { d$text <- d$text[so] }
        d
      })
    } else {
      # old format
      cn <- rownames(colcols);  if(is.null(cn)) { cn <- paste('m',c(1:length(cn)),sep='') }
      colcols <- colcols[,colnames(varinfo$mat),drop=F]
      colcols <- lapply(1:nrows(colcols),function(i) {
        list(data=colcols[i,],legacy=T)
      })
    }
  } else {
    colcols <- list()
  }
  acols <- list(
    "depth"=list(data=colSums(h.dl[[1]]$cd),
                 colors=colorRampPalette(c("white","black"),space="Lab")(1024),
                 quantile.range=0.95),
    "clusters"=list(data=as.factor(cl),
                    colors=rainbow(length(unique(cl)),v=0.6),
                    text=paste('cluster',cl))
  );

  if(!is.null(varinfo$batch)) { acols <- c(acols,list(batch=list(data=as.factor(varinfo$batch),colors=rainbow(length(levels(as.factor(varinfo$batch))),v=0.5)))) }
  colcols <- c(acols,colcols)
  
    
  CairoPNG(file=paste(name,"tamr2.png",sep="."),width=900,height=500)
  pagoda.view.aspects(tamr2,cell.clustering=hc,box=T,labCol=NA,margins=c(0.5,10),row.clustering=NA,col.cols=rbind(clusters=colcols$clusters$colors[as.integer(colcols$clusters$data)]))
  dev.off();
  
  #pagoda.view.aspects(tamr2,cell.clustering=hclust(as.dist(1-cor(tamr$xv)),method='ward'),col.cols=rbind(scol),box=T,labCol=NA,margins=c(0.5,25),ColSideColors.unit.vsize=0.05)
  
  # update the pathway table to include the custom sets
  if(any(grepl("^custom",names(pwpca)))) {
    tamr2$df <- pagoda.top.aspects(pwpca,clpca,return.table=T,z.score=-Inf);
  }
  
  if(!include.aspects) {
    x <- pagoda.cluster.cells(tam,varinfo,include.aspects=F,verbose=T,return.details=T)
  }
  library(Rtsne);
  set.seed(seed);
  #tSNE.pagoda <- tsne(x$distance,k=2,perplexity=20,initial_dims=100)
  tSNE.pagoda <- Rtsne(x$distance,is_distance=T,initial_dims=100,perplexity=perplexity)
  CairoPNG(file=paste(name,"rtSNE.pagoda.png",sep="."),width=350,height=350)
  par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
  plot(tSNE.pagoda$Y,col=clcol,cex=1,pch=19,xlab="PC1",ylab="PC2")
  dev.off()
  emb <- tSNE.pagoda$Y; rownames(emb)<-labels(x$distance);
  
  app <- make.pagoda.app(tamr2,tam,varinfo,env,pwpca,clpca,col.cols=colcols,cell.clustering=hc,title=name,embedding=emb)
  if(!(!missing(port) && is.null(port))) {  # show app unless port was specified to be NULL; missing value is fine
    if(missing(port)) {
      show.app(app,name,browse=browse,port=port)
    } else {
     show.app(app,name,browse=browse)
    }
  }
  saveRDS(app,file=paste(name,"app.rds",sep="."))
  #app <- readRDS(paste(name,"app.rds",sep=".")); 
  #show.app(app,name,port=1461,browse=F)
  
  if(return.details) {
    return(invisible(list(app=app,tam=tam,hc=hc,cl=cl,clcol=clcol,colcols=colcols,tamr=tamr,emb=emb,tamr2=tamr2)))
  } else {
    return(invisible(app))
  }
}


t.process.dataset <- function(dat,name,env=go.env,batch=NULL,k=min(20,ncol(cd)/n.groups),max.model.plots=50,cd=NULL,varinfo=NULL,knn=NULL,max.adj.var=5,skip.pca=FALSE,min.size.entries=1e3,control.for.depth.variation=TRUE,max.quantile=1,seed=0) {
  cat("processing ",name,": ");
  if(is.null(cd)) {
    cd <- dat;
    
    CairoPNG(file=paste(name,"reads.per.cell.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
    hist(colSums(dat)/1e6,col="wheat",xlab="reads per cell (M)",main="Read counts across cells")
    abline(v=1e5/1e6,lty=2,col=2)
    dev.off()
    table(colSums(dat)>=min.cell.reads)
    
    CairoPNG(file=paste(name,"reads.per.gene.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
    hist(log10(rowSums(dat)+1),col="wheat",xlab="reads per gene (log10)",main="Read counts across genes")
    abline(v=log10(min.gene.reads+1),lty=2,col=2)
    dev.off()

    CairoPNG(file=paste(name,"genes.per.cell.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
    hist(log10(colSums(dat>0)+1),col="wheat",xlab="genes per cell (log10)",main="Gene counts across cells")
    abline(v=log10(min.cell.genes+1),lty=2,col=2)
    dev.off()
    
    # filter out low-gene cells
    vi <- colSums(cd)>min.cell.reads; table(vi)
    cd <- cd[,vi]; 
    
    # remove genes that don't have many reads
    vi <- rowSums(cd)>min.gene.reads; table(vi)
    cd <- cd[vi,];
    
    # remove genes that are not seen in a sufficient number of cells
    vi <- rowSums(cd>0)>min.gene.cells; table(vi)
    cd <- cd[vi,];
  }  
  cat("proceeding with ",nrow(cd)," genes across ",ncol(cd)," cells ");
  
  set.seed(seed);
  
  if(is.null(knn) || is.null(varinfo)) {
    knn <- knn.error.models(cd,groups=as.factor(rep(name,ncol(cd))),k=k,n.cores=n.cores,min.count.threshold=1,min.nonfailed=min.nonfailed,verbose=0,max.model.plots=max.model.plots,min.size.entries=min.size.entries)
    cat("models ")
    prior <- scde.expression.prior(models=knn,counts=cd,length.out=400,show.plot=F,max.quantile=max.quantile)
    pdf(file=paste(name,"varnorm.pdf",sep="."),height=4,width=8)
    varinfo <- pagoda.varnorm(knn,counts=cd,trim=trim/ncol(cd),plot=T,verbose=1,prior=prior,max.adj.var=max.adj.var,weight.df.power=1,batch=batch)
    dev.off();
    cat("varinfo ")
    if(control.for.depth.variation) {
      varinfo <- pagoda.subtract.aspect(varinfo,colSums(cd[,rownames(knn)]>0))
    }
  }

  if(!skip.pca) {
    pwpca <- pagoda.pathway.wPCA(varinfo,env,n.components=1,n.cores=n.cores,n.internal.shuffles=0,verbose=1,n.randomizations=5)
    cat("pathways ")
    pdf(file=paste(name,"clvar.pdf",sep="."),height=4,width=8)
    clpca <- pagoda.gene.clusters(varinfo,trim=(trim+5)/ncol(varinfo$mat),n.clusters=150,n.cores=n.cores,verbose=1,plot=T)
    dev.off();
    cat("clusters\n")
  } else {
    pwpca <- clpca <- NULL;
  }
  return(list(cd=cd,knn=knn,prior=prior,varinfo=varinfo,pwpca=pwpca,clpca=clpca,batch=batch))
}
