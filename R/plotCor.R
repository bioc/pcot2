"plotCor" <-
function(x, inputP, av1, av2, av3, range.x, range.dif, nrgcols = 15, labels = FALSE, labcols = 1,main, add.name, font.size, ...){
    par(tcl=0.4)
    n <- ncol(x)
    if (!is.null(inputP)) layout(matrix(c(1:4,4,5,5,5,6,7),5,2), heights=c(1.3,1,3, 5,5), widths=c(12,1)) else if (add.name) layout(matrix(c(1:4,4,5,5,5,6,7),5,2), heights=c(1.5,1,2.5,5,5), widths=c(12,1)) else  layout(matrix(c(1:4,4,5,5,5,6,7),5,2), heights=c(1.5,1,1,5,5), widths=c(12,1)) 

    mm.cor <- max(x)
    
    #1
    if (!is.null(inputP)) par(mar=c(1,3,1.5,1),mgp=c(1,0,0)) else par(mar=c(1,3,2.5,1),mgp=c(1,0,0))
    image(1:n, 1:1, as.matrix(av1), col = wb(nrgcols), axes = FALSE, xlab = "", ylab = "A(a)", zlim=range.x)
    title(main = main)
    box(col="black")
    #2
    par(mar=c(1,3,0,1),mgp=c(1,0,0))
    image(1:n, 1:1, as.matrix(av2), col = wb(nrgcols), axes = FALSE, xlab = "", ylab = "A(b)", zlim=range.x)
    box(col="black")
    #3
    if (!is.null(inputP)){
      pvalues <- round(inputP, digits=2)
      ids <- names(inputP)
      if (add.name){
        par(mar=c(7,3,2,1),mgp=c(1,0.5,0))
        image(1:n, 3:3, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "B", zlim=range.dif)
        axis(1, at = 1:n, labels = labels, las = 2, cex.axis = 0.5, col.axis = labcols, cex=.1, cex.axis=font.size)} else {
          par(mar=c(1,3,2,1),mgp=c(1,0.5,0))
          image(1:n, 3:3, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "B", zlim=range.dif)}
      axis(3, at = 1:n, labels = pvalues[match(labels, ids)], las = 2, cex.axis = 0.5, col.axis = labcols, cex=.1, cex.axis=font.size)  #adding limma pvalues
      ii <- which(pvalues[match(labels, ids)]<=0.05)
      axis(3, at = ii, labels = pvalues[match(labels, ids)][ii], las = 2, cex.axis = 0.5, col.axis = "red", cex=.1, col="red", cex.axis=font.size)  #sig pvalues
    } else if (add.name){
        par(mar=c(7,3,0,1),mgp=c(1,0.5,0))
        image(1:n, 1:1, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "B", zlim=range.dif)
        axis(1, at = 1:n, labels = labels, las = 2, cex.axis = 0.5, col.axis = labcols, cex=.1, cex.axis=font.size)
      } else {
          par(mar=c(1,3,0,1),mgp=c(1,0.5,0))
          image(1:n, 1:1, as.matrix(av3), col = rg(nrgcols), axes = FALSE, xlab = "", ylab = "B", zlim=range.dif)
        }
    box(col="black")
    #4
    par(mar=c(1,3,0,1),mgp=c(1,0,0))
    image(1:n, 1:n, x[,n:1], col = rb(nrgcols), axes = FALSE, xlab = "", ylab = "C", zlim=c(-mm.cor, mm.cor))
    box(col="black")

    ### color scale panel ###
    #5
    par(mar=c(1,2,2.5,3),mgp=c(1,0.5,0))
    image(1:1, 1:nrgcols, t(as.matrix(round(seq(min(range.x), max(range.x), len=nrgcols),1))), ylab = "A", col = wb(nrgcols), axes = FALSE)
    title(main = "Color Key")
    axis(4, at = 1:nrgcols, labels = round(seq(min(range.x), max(range.x), len=nrgcols),1), las = 2, cex.axis = 0.5, col.axis = labcols, cex=.1, cex.axis=font.size)
    box(col="black")
    #6
    par(mar=c(1,2,1,3),mgp=c(1,0.5,0))
    image(1:1, 1:nrgcols, t(as.matrix(round(seq(min(range.dif),max(range.dif), len=nrgcols),1))), ylab = "B", col = rg(nrgcols), axes = FALSE)
    axis(4, at = 1:nrgcols, labels = round(seq(min(range.dif),max(range.dif), len=nrgcols),1), las = 2, cex.axis = 0.5, col.axis = labcols, cex=.1, cex.axis=font.size)
    box(col="black")
    #7
    par(mar=c(2,2,1,3),mgp=c(1,0.5,0))
    image(1:1, 1:nrgcols, t(as.matrix(round(seq(-mm.cor,mm.cor, len=nrgcols),1))), col = rb(nrgcols), ylab = "C", xlab="", axes = FALSE)
    axis(4, at = 1:nrgcols, labels = round(seq(-mm.cor,mm.cor,len=nrgcols),1), las = 2, cex.axis = 0.5, col.axis = labcols, cex=.1, cex.axis=font.size)
    box(col="black")
    
  }

