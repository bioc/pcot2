"aveExprs" <-
function(x, ind.cla, subgene){
  subx <- x[match(subgene, rownames(x)),ind.cla]
  return(rowMeans(subx))
}

