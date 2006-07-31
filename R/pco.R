"pco" <-
function(x,ncomp,dist.method){
  d <- Dist(t(x),method=dist.method)
  pco <- cmdscale(d, k=ncomp, eig=TRUE)
  newx <- pco$points 
  return(newx)
}

