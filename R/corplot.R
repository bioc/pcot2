"corplot" <-
function(x, sel, cla=NULL, inputP=NULL, main, gene.locator=FALSE, add.name=TRUE, font.size=1, dist.method="euclidean"){
  if (!is.null(cla)) clab <- cla else clab <- colnames(x)
  if(length(levels(as.factor(clab)))!=2) stop("The 'class' argument may only have two levels")
  lab <- levels(as.factor(clab)); lab1 <- lab[1]; lab2 <- lab[2]
  ind1 <- which(clab==lab1); ind2 <- which(clab==lab2)

  mat <- x[match(sel,rownames(x)),]
  rr <- Cor(mat, ind1, ind2)
  d <- Dist(rr, method=dist.method)
  hc<-hclust(d,"ave")
  ord<-hc$order
  lname <- hc$label[ord]
  av1 <- aveExprs(mat, ind1, lname)
  av2 <- aveExprs(mat, ind2, lname)
  av3 <- av1-av2

  range.x <- range(x)  
  dif <- rowMeans(x[,ind1])-rowMeans(x[,ind2])
  range.dif <- c(-max(abs(dif)), max(abs(dif)))

  if (!gene.locator) plotCor(rr[ord,ord], inputP, av1, av2, av3, range.x, range.dif, labels=lname, main=main, add.name=add.name, font.size=font.size) else {
    out <- geneLocator(x=rr[ord,ord],main=main)
    cat("You have chosen", length(out), "genes", "\n")
    return(out)}
}

