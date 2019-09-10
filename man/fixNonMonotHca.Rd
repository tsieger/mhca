\name{fixNonMonotHca}
\alias{fixNonMonotHca}
\title{fixNonMonotHca}
\description{Solving non-monotonous heights in HCA clustering.
}
\usage{fixNonMonotHca(hca)}
\arguments{
  \item{hca}{an object of class \code{\link[stats]{hclust}}}
}
\details{\code{\link{fixNonMonotHca}} removes non-monotonous heights in a
dendrogram representing a hierarchical clustering. he non-monotonous
heights get replaced with artificial values lying halfway between
neighbouring monotonous heights.
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
d<-cbind(1:3,1:3)
hd<-mhclust(d)
# original dendrogram
print(hd$height)

hdFixed<-fixNonMonotHca(hd)
# dendrogram with resolved non-monotonous heights
print(hdFixed$height)
}
