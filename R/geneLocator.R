"geneLocator" <-
function(x,nrgcols=15,main=main){
  n <- ncol(x)
  ids <- colnames(x)
  plot.new() 
  par(mar=c(4,3,3,3),mgp=c(1,0,0))
  image(1:n, 1:n, x[,n:1], col = rb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=c(-max(x), max(x)))
  title(main=main)
  axis(1, at = 1:n, labels = ids, las = 2, cex.axis = 0.5, col.axis = 1, cex=.1) 
  aa <- round(locator(2, type="n")$x)
  return(ids[aa[1]:aa[2]])
}

