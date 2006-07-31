"t2.permu" <-
function(x,imat,tstat,ind1,ind2,permu,iter,ncomp,dist.method){
  npath <- ncol(imat)
  tmat <- matrix(0, nrow=iter, ncol=npath)
  
  for (k in 1:iter){
    if (trunc(k/50) == k/50) cat("iter=", k, "\n")
    if (permu=="ByColumn"){
      nn <- ncol(x)
      newD <- x[,sample(1:nn,nn)]}
    if (permu=="ByRow"){
      nn <- nrow(x)
      newD <- x[sample(1:nn,nn),]}
    
    dimnames(newD) <- dimnames(x)
    
    for (i in 1:npath){
      dat <- newD[as.logical(imat[,i]),]
      newdat <- pco(dat,ncomp,dist.method)
      tmat[k,i] <- t2(newdat,ind1,ind2)[1]
    }
  }
  
  p.per <- double(npath)
  for (i in 1:npath)
    p.per[i] <- (sum(abs(tmat[-1,i])>=abs(tstat[i]))+1)/iter
  return(p.per)
}

