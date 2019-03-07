\name{mhclust}
\alias{mhclust}
\title{mhca}
\description{Hierarchical clustering with Mahalanobis inter-cluster distances.
}
\usage{mhclust(x, thresh = 0.5, scale = FALSE, quick = FALSE, g = NULL, 
    normalize = FALSE, verb = 0, useR = FALSE)}
\arguments{
  \item{x}{data.frame or matrix containing coordinates of elementary
observations. Rows correspond to observations, columns correspond
to dimensions of the feature space.}
  \item{thresh}{real number in the interval of (0,1) defining
the minimal relative size of cluster (relative to the number of
observations) whose distance to other clusters will be computed as a
pure Mahalanobis distance. The distance from smaller clusters will
be computed as a mixture of Mahalanobis and Euclidean distances,
with the contribution of the Mahalanobis distance being proportional
to the cluster size. This threshold balances the uncertainty in the
cluster shape as estimated by the covariance matrix of its members.}
  \item{scale}{boolean. Should we transform observations by the
`\code{\link{scale}}' function?}
  \item{quick}{boolean. If \code{TRUE}, inter-cluster
distances will be computed using centroids only. If \code{FALSE},
all observations contained in the clusters will be used.}
  \item{g}{Optional assignment of samples to apriori clusters that
should get formed (in a hierarchical fashion) before any other
merging takes place. By default, there are no apriori clusters, and
clustering starts from individual observations. If `g' is supplied,
a numeric vector of length corresponding to the number of rows of
`x' is expected, holding the index of apriori cluster each sample is
member of. 0's can appear in the vector meaning that the
corresponding sample is not member of any apriori cluster. See the
`apriori' demo for an example.}
  \item{normalize}{boolean. If \code{TRUE}, cluster size
will be ignored when computing Mahalanobis distance from the
cluster. If \code{FALSE}, once all clusters are of at least the
\code{thresh} relative size, both cluster shape and size will
affect inter-cluster distance.}
  \item{verb}{level of verbosity, the greater the more detailed
info, defaults to 0 (no info).}
  \item{useR}{if TRUE, R implementation gets used, otherwise, C
implementation gets used}
}
\details{This is 'mahalanobis-average' hierarchical clustering similar to
`hclust' with advanced merging strategy. The shape of clusters is
considered when computing inter-cluster distances.

The distance between two clusters `c1' and `c2' is the mean of
the distances of members of `c1' to the cluster `c2' and the distances
of members of `c2' to the cluster `c1'. The formula is:

dist(c1,c2) = 1/2 * (
mean ( dist ( c1_i, c2 ) ) +
mean ( dist ( c2_i, c1 ) )
),

where `c1',`c2' denote clusters, `c1_i' iterates over members of `c1',
`c2_i' iterates over members of `c2', and

\eqn{dist(x,c) = sqrt ( ( x - \hat{c} )' * cov(c)^-1 * ( x - \hat{c} ) )}

where `c' is a cluster, `x' is an observation column vector,
\eqn{\hat{c}} is the center of cluster `c' (mean of the observations
contained in `c').

The distance between an individual observation `x' and a cluster
`c' is a mixture of the Mahalanobis and Euclidean distances from
`c' to `x', weighted by the relative cluster size (see the
`thresh' parameter).
}
\value{An object of class *hclust*. The object is a list with
components:

'merge': an n-1 by 2 matrix. The `i'-th row describes the two
clusters merged at the `i'-th step of the clustering. If an
element `j' is negative, then observation `-j' was merged at this
stage.  If `j' is positive, the merge was with the cluster
formed at the (earlier) stage `j' of the algorithm. Negative
entries thus indicate agglomerations of singletons, and
positive entries indicate agglomerations of non-singletons.

'height': a set of `n'-1 non-decreasing real values, the
clustering heights - the distances between the clusters merged
at each step.

'order': a vector giving the permutation of the original
observations suitable for plotting, in the sense that a cluster
plot using this ordering and matrix 'merge' will not have
crossings of the branches.

'labels': labels of each of the objects being clustered.

'method': 'mahalanobis-average'

'dist.method': 'euclidean'
A list with components:

'merge': an n-1 by 2 matrix. The `i'-th row describes the two
clusters merged at the `i'-th step of the clustering. If an
element `j' is negative, then observation `-j' was merged at this
stage.  If `j' is positive, the merge was with the cluster
formed at the (earlier) stage `j' of the algorithm. Negative
entries thus indicate agglomerations of singletons, and
positive entries indicate agglomerations of non-singletons.

'height': a set of `n'-1 non-decreasing real values, the
clustering heights - the distances between the clusters merged
at each step.

'order': a vector giving the permutation of the original
observations suitable for plotting, in the sense that a cluster
plot using this ordering and matrix 'merge' will not have
crossings of the branches.

'labels': labels of each of the objects being clustered.

'method': 'mahalanobis-average'

'dist.method': 'euclidean'}
\references{Fiser K., Sieger T., Schumich A., Wood B., Irving J.,
Mejstrikova E., Dworzak MN. _Detection and monitoring of normal and
leukemic cell populations with hierarchical clustering of flow
cytometry data_. Cytometry A. 2012 Jan;81(1):25-34.
doi:10.1002/cyto.a.21148.
}
\author{Tomas Sieger, Karel Fiser}



\seealso{\code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}}.
}
\examples{
opar<-par(mfrow=c(2,2))

k<-3
n<-nrow(xy)

# classical HCA
h<-hclust(dist(xy))

# Mahalanobis HCA
mh<-mhclust(xy,thresh=.3)

ch<-cutree(h,k=k)
cmh<-cutree(mh,k=k)

# feature space plots with 3 top clusters
plot(xy[,1],xy[,2],asp=1,col=ch,main='HCA')
plot(xy[,1],xy[,2],asp=1,col=cmh,main='Mahalanobis HCA')

# HCA dendrogram
plot(h,hang=0,labels=FALSE,main='Dendrogram of HCA')
y<-min(h$height)-diff(range(h$height))/20
text(1:n,y,(1:n)[h$order],col=ch[h$order],srt=90)

# MHCA dendrogram
plot(mh,labels=FALSE,main='Dendrogram of MHCA')
y<-min(mh$height)-diff(range(mh$height))/10
text(1:n,y,(1:n)[mh$order],col=cmh[mh$order],srt=90)

par(opar)
}
