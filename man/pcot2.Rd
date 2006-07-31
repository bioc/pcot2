\name{pcot2}
\alias{pcot2}
\title{Principle Coordinates and Hotelling's T-Square}
\description{
The \code{pcot2} function implements the PCOT2 testing method, which is a
two-stage permutation-based approach for testing changes in activity in
pre-specified gene sets.
}
\usage{
pcot2(emat, class = NULL, imat, permu = "ByColumn", iter = 1000, alpha = 0.05, adjP.method = "BY", var.equal = TRUE, ncomp = 2, dist.method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{emat}{ A gene expression matrix with no missing values; Each row represents a gene and each column represents a sample. }
  \item{class}{ Class labels representing two distinct experimental conditions (e.g., normal and disease). }
  \item{imat}{ The gene category indicator matrix indicates
presence or absence of genes in pre-defined gene sets
(e.g., gene pathways). The indicator matrix contains rows representing
gene identifiers of genes present in the expression data and columns
representing pre-defined group names. A value of 1 indicates the
presence of a gene and 0 indicates the absence for the gene in a
particular group.}
  \item{permu}{ Specifies whether genes or samples are permuted.  By
    default, permutations are performed by sample ("ByColumn").}
  \item{iter}{ The number indicates how many permutations will be
    performed in the analysis. }
  \item{alpha}{ alpha determines the significance threshold for the permutation p-values. }
  \item{adjP.method}{ Specifies that p-values be adjusted by one of
    the following methods: "bonferroni", "holm", "hochberg", "hommel",
    "BH" (Benjamini and Hochberg), or "BY"  (Benjamini and Yekutieli).}
  \item{var.equal}{ Specifies the use of either a pooled estimate of
    correlation for the two classes or an unpooled estimate for
    calculating each T-squared statistic. By default, the pooled estimate is used.}
  \item{ncomp}{  The dimensionality to which the data matrix is reduced
    via principle coordinates. The default dimensionality is set as
    \code{ncomp=2}.}
   \item{dist.method}{Specifies the method for calculating
    distance in the PCO procedure.  The available distance methods are
    "euclidean", "maximum", "manhattan", "canberra", "binary",
    "pearson","correlation" or "spearman". For additional details see the
    \code{amap} package and the help documentation for the \code{Dist} function.}
}
\details{
  The raw permutation p-values are adjusted for multiple testing by a call to 'p.adjust'.
}
\value{
  \item{res.all}{A data frame which prints information for all pathways}
  \item{res.sig}{A data frame which prints information for significant pathways at a given alpha level}
  \item{comparison}{Print the contrast used in the analysis}
  ...
}
\references{S.Song and M.Black (2006) PCOT2: an R package for assessing
  expression changes in groups of related genes. Submitted to Bioinformatics.}
\author{ Sarah Song and Mik Black }

\seealso{\code{\link{corplot}},\code{\link{corplot2}},\code{\link{aveProbe}}}

\examples{
ns <- 40  ## 40 samples
cla <- rep(c("Trt","Ctr"),each=ns/2)
ngene <- 10  ## 10 genes per group 
npath <- 10  ## 10 groups

nreal <- 3  ## alter groups ##
nnull <- npath-nreal   ## null groups ##
pname <- c(paste("RealP",1:nreal, sep=""), paste("NullP",1:nnull, sep=""))

## Three main inputs in the function ##
## [1] Simulate (gene) expression matrix (emat) ##
rmv <- function(mn, covm, nr, nc){
   sigma <- diag(nr)
   sigma[sigma==0] <- covm
   x1 <- rmvnorm(nc/2, mean=mn, sigma=sigma)
   x0 <- rmvnorm(nc/2, mean=rep(0,nr), sigma=sigma)
   mat <- t(rbind(x1,x0))
  return(mat)
}

covm <- 0.9  ##covariance 
ct <- c(6,8,10)  ##mean

library(mvtnorm)
emat <- c()
for (i in 1:nreal) emat <- rbind(emat, rmv(rep(ct[i],ngene),covm=covm, ngene, ns))  # for alt pathways
for (i in 1:(npath-nreal)) emat <- rbind(emat, rmv(mn=rep(0,ngene),covm=covm, nr=ngene, nc=ns))
dimnames(emat) <- list(paste("Gene", 1:(ngene*npath),sep=""), cla)

## [2] class label ##
cla

## [3] indicator matrix (row: genes and col: pathways)
imat <- kronecker(diag(npath),rep(1,ngene))
dimnames(imat) <- list(paste("Gene",1:(ngene*npath), sep=""), pname)

results.pcot2 <- pcot2(emat, cla, imat)
results.pcot2$res.sig
results.pcot2$res.all

}
\keyword{methods}

