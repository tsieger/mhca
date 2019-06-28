\name{cutreeApriori}
\alias{cutreeApriori}
\title{Cut off the non-apriori portion of a tree resulting from}
\description{This function restricts a tree resulting from \code{mhclust} to the
portion representing the mergings of the apriori clusters.}
\usage{cutreeApriori(h, verb = 0)}
\arguments{
  \item{h}{an object of class *hclust*, resulting from \code{mhclust}.}
  \item{verb}{verbosity level}
}

\value{An object of class *hclust*. The object corresponds to the
original tree \code{h} having the non-apriori portion removed.
See \code{\link{mhclust}} \code{\link{hclust}} for details.}

\author{Tomas Sieger, Karel Fiser}




\examples{
# demo data to cluster
x<-c(1,2,4,10,13,20,24)
x<-cbind(x,x)
h<-mhclust(x,g=c(1,1,1,2,2,3,3))
h2<-cutreeApriori(h)
opar<-par(mfrow=c(1,2))
plot(h)
plot(h2)
par(opar)
}
