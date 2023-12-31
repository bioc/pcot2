\name{corplot2}
\alias{corplot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce a plot for jointly visualizing unpooled correlation
  information and expression data for selected genes}
\description{This plot is used for looking at unpooled inter-gene correlation within a
pre-defined group of genes, in conjunction with information about
differences in expression activity between classes.}
\usage{
corplot2(x, sel, cla = NULL, inputP = NULL, main, gene.locator = FALSE, add.name = TRUE, font.size = 1, dist.method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix with no missing values; Each row represents a gene and each column represents a sample. }
  \item{sel}{ A vector of selected gene identifiers. }
  \item{cla}{ Class labels representing two distinct experimental conditions (e.g., normal and disease). }
 \item{inputP}{ This option allows users to input p-values for each
    gene (e.g., if produced by another software package).}
  \item{main}{ A title for the plot. }
  \item{gene.locator}{ This option allows users to click of the plot to
    identify groups of genes.  Clicking twice on the diagonal of the
    plot returns the identifiers of genes between the points clicked.}
  \item{add.name}{Specifies whether gene identifiers should be printed
    on the plot.}
  \item{font.size}{Adjusts the size of gene names printed on the plot. }
  \item{dist.method}{Specifies the method for calculating
    inter-gene distance (used when ordering the rows and columns of the
    correlation plot).  The available distance methods are
    "euclidean", "maximum", "manhattan", "canberra", "binary",
    "pearson","correlation" or "spearman". For additional details see the
    \code{amap} package and the help documentation for the \code{Dist} function.}
}

\author{ Sarah Song and Mik Black }

\seealso{\code{\link{pcot2}},\code{\link{corplot}},\code{\link{aveProbe}}}


\examples{
library(multtest)
library(hu6800.db)  
data(golub)
rownames(golub) <- golub.gnames[,3]
colnames(golub) <- golub.cl
KEGG.list <- as.list(hu6800PATH)
imat <- getImat(golub, KEGG.list, ms=10) 
colnames(imat) <- paste("KEGG", colnames(imat), sep="")
sel <- c("04620","04120")
pvalue <- c(0.001, 0.72)
library(KEGG.db)
pname <- unlist(mget(sel, env=KEGGPATHID2NAME))
main <- paste("KEGG", sel, ": ", pname, ": ", "P=", pvalue, sep="")
for(i in 1:length(sel)){
fname <- paste("corplot2-KEGG",sel[i] , ".jpg", sep="")
jpeg(fname, width=1600, height=1200, quality=100)
selgene <- rownames(imat)[imat[,match(paste("KEGG",sel,sep="")[i],colnames(imat))]==1]
corplot2(golub, selgene, golub.cl, main=main[i])
dev.off()
}

}
\keyword{ hplot }
\keyword{ iplot }
