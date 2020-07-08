\name{fixNonMonotHca}
\alias{fixNonMonotHca}
\title{fixNonMonotHca}
\description{Solving non-monotonous heights in HCA clustering.
}
\usage{fixNonMonotHca(hca, method = c("eps", "halfway"), eps = NULL)}
\arguments{
  \item{hca}{an object of class \code{\link[stats]{hclust}}}
  \item{method}{replacement method, one of 'eps', and
'halfway'.
The 'eps' method makes the non-monotonous heights to be monotonous
by shifting them just above the last monotonous height (\code{eps}
above it). The resulting heights are close to reality,
eventhough the denrogram is not so visually appealing.
The 'halfway' method replaces the non-monotonous heights with
heights lying halfway between neighbouring monotonous heights. This
gives nicer dendrogram at the expense of making the heights too
artificial. See examples.}
  \item{eps}{parameter of the 'eps' method, if \code{NULL}, it
defaults to half of the minimal difference in monotonous heights, or
\code{1}, if there are no consecutive monotonous heights.}
}
\details{\code{\link{fixNonMonotHca}} replaces non-monotonous heights in a
dendrogram representing a hierarchical clustering with artificial
values making the heights monotonous.
However, please note that non-monotonous heights appear naturally in
dendrograms resulting from some HCA methods including
\code{link{mhclust}}, such that altering the heights result in
artificial dendrograms that do NOT represent the structure of the
clustered data. Please consider using HCA methods producing
monotonous before attempting to remove the non-monotonous heights,
e.g. for the purpose of graphical presentation of the dendrogram.}
\value{An object of class \code{link[stats]{hclust}} having the
\code{height} component fixed.}

\author{Tomas Sieger, Karel Fiser}




\examples{
# simple example
d<-cbind(1:3,1:3)
hd<-mhclust(d)
# original dendrogram
print(hd$height)
hdFixed<-fixNonMonotHca(hd)
# dendrogram with resolved non-monotonous heights
print(hdFixed$height)

# another example
set.seed(1)
x<-cbind(runif(20),runif(20))
hx<-hclust(dist(x)^2,'cen')
hx1<-fixNonMonotHca(hx,method='eps')
hx2<-fixNonMonotHca(hx,method='halfway')
opar<-par(mfrow=c(1,3))
plot(hx,main='Original')
plot(hx1,main='Fixed by the "eps" method')
plot(hx2,main='Fixed by the "halfway" method')
par(opar)
}
