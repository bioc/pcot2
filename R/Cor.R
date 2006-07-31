"Cor" <-
function(x, ind1, ind2){
  dat.n <- x[,ind1]
  dat.d <- x[,ind2]
  nn <- length(ind1)
  nd <- length(ind2) 
  SN <- var(t(dat.n))
  SD <- var(t(dat.d))
  S <- ((nd-1)*SD + (nn-1)*SN)/(nd+nn-2)
  return(cov2cor(S))
}

