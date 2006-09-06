"corplot2" <-
function(x, sel, cla=NULL, inputP=NULL, main, gene.locator=FALSE, add.name=TRUE, font.size=1, dist.method="euclidean"){

    lab <- levels(as.factor(cla)); lab1 <- lab[1]; lab2 <- lab[2]
    ind1 <- which(cla==lab1); ind2 <- which(cla==lab2)
    library(amap)
    mm.cor<- max(c(cor(t(x[match(sel, rownames(x)),ind1])),cor(t(x[match(sel, rownames(x)),ind2]))))
    
  if (!gene.locator){
    if (!is.null(inputP)) layout(matrix(c(14,1,3,11,14,5,8,12,14,7,10,12,14,6,9,12,14,4,2,13),4,5), heights=c(0.8,5,5,1), widths=c(4,0.3,0.6,0.3,4)) else if (add.name) layout(matrix(c(14,1,3,11,14,5,8,12,14,7,10,12,14,6,9,12,14,4,2,13),4,5), heights=c(0.8,5,5,1), widths=c(4,0.3,0.5,0.3,4)) else layout(matrix(c(14,1,3,11,14,5,8,12,14,7,10,12,14,6,9,12,14,4,2,13),4,5), heights=c(0.8,5,5,1), widths=c(4,0.3,0.3,0.3,4))
    
    nrgcols <- 15
    par(tcl=0.4)
    
    #1:ClassI
    par(mar=c(1,3,1,1),mgp=c(1,0,0))
    rr <- cor(t(x[match(sel, rownames(x)),ind1]))
    hc <- hclust(Dist(rr, method=dist.method),"ave")
    ord1<-hc$order
    rr.ord <- rr[ord1, ord1]; n <- ncol(rr.ord)
    image(1:n, 1:n, rr.ord[,n:1], col = rb(nrgcols), axes = FALSE, xlab = "", ylab = paste("Ordered by genes in Class",lab1), zlim=c(-mm.cor, mm.cor), main=paste("Class", lab1))

    #2:ClassIV
    par(mar=c(1,1,1,3),mgp=c(0,0,0))
    rr <- cor(t(x[match(sel, rownames(x)),ind2]))
    hc <- hclust(Dist(rr, method=dist.method),"ave")
    ord2<-hc$order
    rr.ord <- rr[ord2, ord2]; n <- ncol(rr.ord)
    image(1:n, 1:n, rr.ord[,n:1], col = rb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=c(-mm.cor, mm.cor), main=paste("Class", lab2))

    #3:ClassI, ordered by genes in ClassIV
    par(mar=c(1,3,1,1),mgp=c(1,0,0))
    rr <- cor(t(x[match(sel, rownames(x)),ind1]))
    hc <- hclust(Dist(rr, method=dist.method),"ave")
    rr.ord <- rr[ord2, ord2]; n <- ncol(rr.ord)
    image(1:n, 1:n, rr.ord[,n:1], col = rb(nrgcols), axes = FALSE, xlab = "", ylab = paste("Ordered by genes in Class",lab2), zlim=c(-mm.cor, mm.cor), main=paste("Class", lab1))

    #4:ClassIV, ordered by genes in ClassI
    par(mar=c(1,1,1,3),mgp=c(0,0,0))
    rr <- cor(t(x[match(sel, rownames(x)),ind2]))
    hc <- hclust(Dist(rr, method=dist.method),"ave")
    rr.ord <- rr[ord1, ord1]; n <- ncol(rr.ord)
    image(1:n, 1:n, rr.ord[,n:1], col = rb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=c(-mm.cor, mm.cor), main=paste("Class", lab2))

    #5:average expression profiles in ClassI, ordered by ClassI
    par(mar=c(1,0,1,1.5),mgp=c(0,0,0))
    range.x <- range(x)
    av1 <- t(rowMeans(x[match(hc$label[ord1], rownames(x)),ind1])[n:1])
    image(1:1, 1:n, as.matrix(av1), col = wb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.x)
    title(paste("",lab1, sep=""))
    box(col="black")

    #6:average expression profiles in ClassIV, ordered by ClassI
    par(mar=c(1,1.5,1,0),mgp=c(0,0,0))
    av2 <- t(rowMeans(x[match(hc$label[ord1], rownames(x)),ind2])[n:1])
    image(1:1, 1:n, as.matrix(av2), col = wb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.x)
    title(paste("",lab2, sep=""))
    box(col="black")

    #7:Difference, ordered by classI
    av3 <- av1-av2
    dif <- rowMeans(x[,ind1])-rowMeans(x[,ind2])  #difference for all the genes
    range.dif <- c(-max(abs(dif)), max(abs(dif)))  
    if(!is.null(inputP)){
       if (add.name) {
         par(mar=c(1,4,1,1),mgp=c(1,0.5,0))
         image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)
         axis(2, at = n:1, labels = hc$label[ord1], las = 2, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)} else {
           par(mar=c(1,0,1,1),mgp=c(1,0.5,0))
           image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)}       
       axis(4, at = n:1, labels = inputP[match(hc$label[ord1],names(inputP))], las = 2, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
       ii <- which(inputP[match(hc$label[ord1],names(inputP))]<=0.05)
       if(length(ii)!=0) axis(4, at = length(hc$label)+1-ii, labels = inputP[match(hc$label[ord1],names(inputP))][ii], las = 2, cex.axis = 0.5, col.axis = "red", cex=.1, col="red", cex.axis=font.size)} else if (add.name){
         par(mar=c(1,4,1,0),mgp=c(1,0.5,0))
         image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)
         axis(2, at = n:1, labels = hc$label[ord1], las = 2, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
       } else {
         par(mar=c(1,0,1,0),mgp=c(1,0.5,0))
         image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)
       }   
    title(paste(lab1, "-", lab2, sep=""))
    box(col="black")

    #8:average expression profiles in ClassI, ordered by ClassIV
    par(mar=c(1,0,1,1.5),mgp=c(0,0,0))
    av1 <- t(rowMeans(x[match(hc$label[ord2], rownames(x)),ind1])[n:1])
    image(1:1, 1:n, as.matrix(av1), col = wb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.x)
    box(col="black")

    #9:average expression profiles in ClassIV, ordered by ClassIV
    par(mar=c(1,1.5,1,0),mgp=c(0,0,0))
    av2 <- t(rowMeans(x[match(hc$label[ord2], rownames(x)),ind2])[n:1])
    image(1:1, 1:n, as.matrix(av2), col = wb(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.x)
    box(col="black")

    #10:Difference, ordered by classIV
    av3 <- av1-av2  
    if(!is.null(inputP)){
       if (add.name){
         par(mar=c(1,4,1,1),mgp=c(1,0.5,0))
         image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)
         axis(2, at = n:1, labels = hc$label[ord2], las = 2, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)} else {
           par(mar=c(1,0,1,1),mgp=c(1,0.5,0))
           image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)}
       axis(4, at = n:1, labels = inputP[match(hc$label[ord2],names(inputP))], las = 2, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
       ii <- which(inputP[match(hc$label[ord2],names(inputP))]<=0.05)
       if(length(ii)!=0) axis(4, at = length(hc$label)+1-ii, labels = inputP[match(hc$label[ord2],names(inputP))][ii], las = 2, cex.axis = 0.5, col.axis = "red", cex=.1, col="red", cex.axis=font.size)
    } else if (add.name){
         par(mar=c(1,4,1,0),mgp=c(1,0.5,0))
         image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)
         axis(2, at = n:1, labels = hc$label[ord2], las = 2, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
       } else {
         par(mar=c(1,0,1,0),mgp=c(1,0.5,0))
         image(1:1, 1:n, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "", zlim=range.dif)
       }
    box(col="black")

    #11: color key (grey)
    par(mar=c(3,3,1,1),mgp=c(1,0.5,0))
    image(1:nrgcols, 1:1,  as.matrix(round(seq(min(range.x), max(range.x), len=nrgcols),1)), col = wb(nrgcols),axes = FALSE, xlab="", ylab = "")
    axis(1, at = 1:nrgcols, labels = round(seq(min(range.x), max(range.x), len=nrgcols),1), las = 1, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
    box(col="black")

    #12: color key (red-green)
    par(mar=c(3,0,1,0),mgp=c(1,0.5,0))
    image(1:nrgcols, 1:1, as.matrix(round(seq(min(range.dif),max(range.dif), len=nrgcols),1)), col = rg(nrgcols), axes = FALSE, xlab="", ylab = "")
    axis(1, at = 1:nrgcols, labels = round(seq(min(range.dif),max(range.dif), len=nrgcols),1), las = 1, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
    box(col="black")

    #13: color key (red-blue)
    par(mar=c(3,1,1,2),mgp=c(1,0.5,0))
    image(1:nrgcols,1:1, as.matrix(round(seq(-mm.cor,mm.cor, len=nrgcols),1)), col = rb(nrgcols),axes = FALSE, xlab="", ylab = "")
    axis(1, at = 1:nrgcols, labels = round(seq(-mm.cor,mm.cor,len=nrgcols),1), las = 1, cex.axis = 0.5, col.axis = 1, cex=.1, cex.axis=font.size)
    box(col="black")

    #14: title of the whole plot
    plot.new()
    par(mar=c(0,3,3,3),mgp=c(1,0,0))
    title(main = list(main, cex=1.8, col="Black", font=2))
  }
  if (gene.locator){
    ### In Class1
    rr <- cor(t(x[match(sel, rownames(x)),ind1]))
    hc <- hclust(Dist(rr, method=dist.method),"ave")
    ord1<-hc$order
    rr.ord <- rr[ord1, ord1]
    out1 <- geneLocator(rr.ord,main=paste(main,"---Class", lab1, sep=""))
    
    ### In Class2
    rr <- cor(t(x[match(sel, rownames(x)),ind2]))
    hc <- hclust(Dist(rr, method=dist.method),"ave")
    ord2<-hc$order
    rr.ord <- rr[ord2, ord2]
    out2 <- geneLocator(rr.ord,main=paste(main,"---Class", lab2, sep=""))

    cat("You have chosen", length(out1), "genes in", paste("Class", lab1, sep=""), "\n")
    cat("You have chosen", length(out2), "genes in", paste("Class", lab2, sep=""), "\n")
    list(out1=out1, out2=out2)
  }
  }

