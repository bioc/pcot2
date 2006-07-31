"t2.unequ" <-
function(x, ind1, ind2){
  dat1 <- x[ind1,]
  dat2 <- x[ind2,]
  nn <- length(ind1)
  nd <- length(ind2) 
  SN <- var(dat1)
  SD <- var(dat2)
  S <- SD/nd + SN/nn

  xbar1 <- apply(dat1, 2, mean)
  xbar2 <- apply(dat2, 2, mean)
  xdiff <- xbar1 - xbar2
  t2 <- xdiff %*% solve(S) %*% xdiff
  return(t2)
}

