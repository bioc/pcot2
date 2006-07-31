"pcot2" <-
function(emat, class=NULL, imat, permu="ByColumn", iter=1000, alpha=0.05, adjP.method="BY", var.equal=TRUE, ncomp=2, dist.method="euclidean"){
   if (!is.null(class)) cla <- class else cla <- colnames(emat)
   if (length(cla)!=ncol(emat)) stop("The length of the 'class' argument must equal the number of columns (samples) in the 'emat' (expression matrix) argument")
   if (length(levels(as.factor(cla)))!=2) stop("The 'class' argument may only have two levels")
   lab <- levels(as.factor(cla)); lab1 <- lab[1]; lab2 <- lab[2]; lab.print <- paste("Comparison: ", lab1, "-", lab2, sep="")
   ind1 <- which(cla==lab1); ind2 <- which(cla==lab2)
    
   npath <- ncol(imat)
   num <- apply(imat, 2, sum)
   tstat <- p.nor <- double(npath)

   for (i in 1:npath){
    dat <- emat[as.logical(imat[,i]),]
    newdat <- pco(dat, ncomp,dist.method)
    if (var.equal){
      output <- t2(x=newdat,ind1,ind2)
      tstat[i] <- output[1]
      p.nor[i] <- output[2]} else tstat[i] <- t2.unequ(newdat,ind1,ind2)
  }
   if(var.equal) p.adj <- p.adjust(p.nor, adjP.method)
   p.per <- t2.permu(emat,imat,tstat,ind1,ind2,permu,iter,ncomp,dist.method)
   p.per.adj <- p.adjust(p.per, adjP.method)

   if(var.equal) res <- data.frame(Num=num, T2=tstat, P.nor=p.nor, P.adj=p.adj, P.permu=p.per, P.permu.adj=p.per.adj) else res <- data.frame(Num=num, T2=tstat, P.permu=p.per, P.permu.adj=p.per.adj)
   res.all <- res[order(p.per),]                
   res.sig <- res.all[res.all$P.permu<=alpha,]
   cat(lab.print,"\n")
   list(res.all=res.all, res.sig=res.sig)
 }

