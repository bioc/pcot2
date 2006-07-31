###################################################
### chunk number 1: 
###################################################
library(pcot2)
library(multtest)
library(hu6800)  
set.seed(1234567)


###################################################
### chunk number 2: 
###################################################
data(golub)
rownames(golub) <- golub.gnames[,3]
colnames(golub) <- golub.cl


###################################################
### chunk number 3: 
###################################################
golub.cl


###################################################
### chunk number 4: 
###################################################
KEGG.list <- as.list(hu6800PATH)
imat <- getImat(golub, KEGG.list, ms=10) 
colnames(imat) <- paste("KEGG", colnames(imat), sep="")
dim(imat)


###################################################
### chunk number 5: 
###################################################
results <- pcot2(golub, golub.cl, imat, iter=10)


###################################################
### chunk number 6: 
###################################################
results$res.sig
results$res.all


###################################################
### chunk number 7: 
###################################################
sel <- c("04620","04120")
pvalue <- c(0.001, 0.72)
library(KEGG)
pname <- unlist(mget(sel, env=KEGGPATHID2NAME))
main <- paste("KEGG", sel, ": ", pname, ": ", "P=", pvalue, sep="")
for(i in 1:length(sel)){
fname <- paste("corplot2-KEGG",sel[i] , ".jpg", sep="")
jpeg(fname, width=1600, height=1200, quality=100)
selgene <- rownames(imat)[imat[,match(paste("KEGG",sel,sep="")[i],colnames(imat))]==1]
corplot2(golub, selgene, golub.cl, main=main[i])
dev.off()
}


###################################################
### chunk number 8: 
###################################################
pathlist <- as.list(hu6800PATH)
pathlist <- pathlist[match(rownames(golub), names(pathlist))]
ids <- unlist(mget(names(pathlist), env=hu6800SYMBOL))
#### transform data matrix only ####
newdata <- aveProbe(x=golub, ids=ids)$newx
#### transform both data and imat ####
output <- aveProbe(x=golub, imat=imat, ids=ids)
newdata <- output$newx
newimat <- output$newimat
newimat <- newimat[,apply(newimat, 2, sum)>=10]
dim(newdata)
dim(newimat)


