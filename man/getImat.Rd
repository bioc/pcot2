\name{getImat}
\alias{getImat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate an indicator matrix}
\description{
This function is used to generate an indicator matrix as an input to
the pcot2 function. The gene category indicator matrix indicates
presence or absence of genes in pre-defined gene sets
(e.g., gene pathways). The indicator matrix contains rows representing
gene identifiers of genes present in the expression data and columns
representing pre-defined group names. A value of 1 indicates the
presence of a gene and 0 indicates the absence for the gene in a
particular group.  
}
\usage{
getImat(x, pathlist, ms = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix with no missing values; Each row represents a gene
    and each column represents a sample.}
  \item{pathlist}{ A list of gene sets. }
  \item{ms}{ The minimum gene set size. Gene sets containing less
    than this number of genes will be excluded from the analysis.}
}

\value{
  An indicator matrix is returned. The matrix value is 1 (gene in) or 0 (gene out)}

\author{ Sarah Song and Mik Black }

\seealso{\code{\link{pcot2}},\code{\link{corplot}},\code{\link{corplot2}},\code{\link{aveProbe}}}
\examples{
library(multtest)
library(hu6800)  
data(golub)
rownames(golub) <- golub.gnames[,3]
colnames(golub) <- golub.cl
KEGG.list <- as.list(hu6800PATH)
imat <- getImat(golub, KEGG.list, ms=10) 
}

\keyword{htest}
