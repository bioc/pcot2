"getImat" <-
function(x, pathlist, ms=10){
  gname1 <- rownames(x); gname2 <- names(pathlist)
  gname <- intersect(gname1, gname2)
  if (is.na(gname[1])) stop("The row (gene) names of the data matrix do not match the names in the 'pathlist' argument")
  gset <- unique(unlist(pathlist))
  pname <- gset[!is.na(gset)]
  imat <- matrix(0, nrow=length(gname), ncol=length(pname))

  for (i in 1:length(gname)){
    ind <- match(gname[i],gname2)
    if(!is.na(pathlist[[ind]][1])) imat[i,match(pathlist[[ind]], pname)] <- 1
  }
  dimnames(imat) <- list(gname, pname)
  imat <- imat[,apply(imat,2,sum)>=ms]
  return(imat)
}

