"aveProbe" <-
function(x, imat=NULL, ids){
  if (!setequal(rownames(x), names(ids))) cat("Warnings: The names of identifiers should be matched up with the names of probes in the data", "\n")
  if (is.list(ids)) cat("Warnings: The variable ids should be converted to a vector", "\n")
  uu <- unique(ids);  uu <- uu[!is.na(uu)]   
  if (is.null(imat)){
    newx <- matrix(0, nrow=length(uu), ncol=ncol(x))
    rownames(newx) <- uu; colnames(newx) <- colnames(x)
    for (i in 1:length(uu)){
      ind <- which(uu[i]==ids)
      if (length(ind)==1) newx[i,] <- x[ind,] else newx[i,] <- apply(x[ind,],2,median)}
    list(newx=newx)}
  else{
    newx <- matrix(0, nrow=length(uu), ncol=ncol(x))
    newimat <- matrix(0, nrow=length(uu), ncol=ncol(imat))
    rownames(newx) <- rownames(newimat) <- uu
    colnames(newx) <- colnames(x); colnames(newimat) <- colnames(imat)    
    for (i in 1:length(uu)){
      ind <- which(uu[i]==ids)
      if (length(ind)==1) newx[i,] <- x[ind,] else newx[i,] <- apply(x[ind,],2,median)
      if (length(ind)==1) newimat[i,] <- imat[ind,] else newimat[i,] <- apply(imat[ind,],2,median)}
    list(newx=newx, newimat=newimat)}
}

