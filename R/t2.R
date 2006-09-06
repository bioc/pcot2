"t2" <-
function(x, ind1, ind2){
  dat1 <- x[ind1,]
  dat2 <- x[ind2,]
  nn <- length(ind1)
  nd <- length(ind2) 
  SN <- var(dat1)
  SD <- var(dat2)
  S <- ((nd-1)*SD + (nn-1)*SN)/(nd+nn-2)

  xbar1 <- colMeans(dat1)
  xbar2 <- colMeans(dat2)
  xdiff <- xbar1 - xbar2

  t2 <- ((nd*nn)/(nd+nn))*(xdiff %*% solve(S) %*% xdiff)
  p <- length(xdiff)
  f <- (nd+nn-p-1)*t2/((nd+nn-2)*p)
  pvalue <- 1-pf(f,p, nd+nn-p-1)
  return(cbind(t2, pvalue))  
}

