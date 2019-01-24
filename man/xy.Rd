\name{xy}
\alias{xy}
\docType{data}
\title{
mhca demo data
}
\description{
Data to demonstrate Mahalanobis hierarchical cluster analysis. 
The consist of a 150x2 numeric matrix containing three 2D elliptical
clusters, each comprising of 50 observations.
}
\format{
A 150x2 numeric matrix of 150 observations and 2 features.
}
\usage{
data(xy)
}
\examples{
if (interactive()) {
    library(mhca)
    data(xy)
    mh<-mhclust(xy,thresh=.3)
    plot(xy[,1],xy[,2],asp=1,col=cutree(mh,2))
}
}
\keyword{datasets}
