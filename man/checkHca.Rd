\name{checkHca}
\alias{checkHca}
\title{Check HCA consistency.}
\description{This function checks that the results of HCA is consistent. It
prints info on key components.}
\usage{checkHca(h, verb = TRUE, dbg = 0)}
\arguments{
  \item{h}{hierarchical clustering (e.g. the result of
\code{\link[stats]{hclust}} or \code{\link[mhca]{mhclust}}).}
  \item{verb}{verbosity level}
  \item{dbg}{debug level}
}

\value{A boolean value indicating whether \code{h} is valid or not.}

\author{Tomas Sieger, Karel Fiser}




\examples{
hGood<-hclust(dist(1:3))
checkHca(hGood)

hBad<-list(merge=rbind(c(-1,-2),c(-1,2)),height=1:2,order=1:3,method='fake-method',dist.method='fake-dist.method')
checkHca(hBad)
}
