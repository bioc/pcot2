\name{aveProbe}
\alias{aveProbe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transform Affymetrix data so that unique genes with multiple
  probes are represented by a single expression value on each array.}
\description{
  In Affymetrix gene expression data, a unique gene can often link to multiple
probe sets, with such genes then having a greater influence on the
analysis (particularly if the gene is differentially expressed).  To
overcome this problem the median is taken across all probes sets which
represent a unique gene.   
}
\usage{
aveProbe(x, imat = NULL, ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix with no missing values; Each row represents a gene and each column represents a sample. }
  \item{imat}{ A matrix indicating presence or absence of genes in the
    gene sets. The indicator matrix contains rows representing gene
    identifiers of genes present in the expression data and columns
    representing group (gene set) names. }
  \item{ids}{ A vector of identifiers (e.g., UniGene or LocusLink
    identifiers) representing unique genes which match to the probe ids
    in the expression data. }
}

\value{
  \item{newx}{A data matrix with rows representing the input
    identifiers and columns representing samples.}
  \item{newimat}{A new imat (indicator matrix) with rows representing
    the unique gene identifiers and columns representing gene sets. }
}
\references{S.Song and M.Black (2006) PCOT2: an R package for assessing
  expression changes in groups of related genes. Submitted to Bioinformatics.}
\author{ Sarah Song and Mik Black }

\seealso{\code{\link{pcot2}},\code{\link{corplot}},\code{\link{corplot2}}}

\examples{
library(multtest)
library(hu6800)  
data(golub)
rownames(golub) <- golub.gnames[,3]
colnames(golub) <- golub.cl
KEGG.list <- as.list(hu6800PATH)
imat <- getImat(golub, KEGG.list, ms=10) 
colnames(imat) <- paste("KEGG", colnames(imat), sep="")


pathlist <- as.list(hu6800PATH)
pathlist <- pathlist[match(rownames(golub), names(pathlist))]
ids <- unlist(mget(names(pathlist), env=hu6800SYMBOL))
#### transform data matrix only ####
newdat <- aveProbe(x=golub, ids=ids)$newx
#### transform both data and imat ####
output <- aveProbe(x=golub, imat=imat, ids=ids)
newdat <- output$newx
newimat <- output$newimat
newimat <- newimat[,apply(newimat, 2, sum)>=10]
}

\keyword{ htest }

