"t2.permu" <-
function(x,imat,tstat,ind1,ind2,permu,iter,ncomp,dist.method){
  npath <- ncol(imat)
  tmat <- matrix(0, nrow = iter, ncol = npath)
  tmat[1,] <- tstat
  for (k in 2:iter) {
    if (trunc(k/50) == k/50)
      cat("iter=", k, "\n")
    if (permu == "ByColumn") {
      nn <- ncol(x)
      newD <- x[, sample(1:nn, nn)]
    }
    if (permu == "ByRow") {
      nn <- nrow(x)
      newD <- x[sample(1:nn, nn), ]
    }
    dimnames(newD) <- dimnames(x)
    for (i in 1:npath) {
      dat <- newD[as.logical(imat[, i]), ]
      newdat <- pco(dat, ncomp, dist.method)
      tmat[k, i] <- t2(newdat, ind1, ind2)[1]
    }
  }
  p.per <- c()
  for (i in 1:npath) p.per[i] <- sum(abs(tmat[,i]) >= abs(tmat[1,i]))/iter
  return(p.per)
}

