\name{mhclust}
\alias{mhclust}
\title{mhca}
\description{Hierarchical clustering with Mahalanobis inter-cluster distances.
}
\usage{mhclust(x, thresh = 0.5, 
    scale = FALSE, quick = FALSE, 
    g = NULL, gIntra = TRUE, 
    gParallel = TRUE, 
    normalize = FALSE, 
    subthreshHandling = c("mahal", 
        "mahal0", "euclid", 
        "euclidMahal"), 
    warn = TRUE, verb = 0, 
    verbRecursive = verb, 
    useR = FALSE, nFull = nrow(as.matrix(x)))}
\arguments{
  \item{x}{data.frame or matrix containing coordinates of elementary
observations. Rows correspond to observations, columns correspond
to dimensions of the feature space.}
  \item{thresh}{real number in the interval of (0,1) defining
the minimal relative size of cluster (relative to the number of
observations) whose distance to other clusters will be computed as a
pure Mahalanobis distance. The distance from smaller clusters will
be computed as a mixture of Mahalanobis and Euclidean distances,
with the contribution of the Mahalanobis distance being larger
for larger clusters, in general. However, the exact handling of
subthreshold clusters is controlled by the \code{subthreshHandling}
argument (see below). This threshold balances the uncertainty in the
cluster shape as estimated by the covariance matrix of its members.}
  \item{scale}{boolean. Should we transform observations by the
`\code{\link{scale}}' function?}
  \item{quick}{boolean. If \code{TRUE}, inter-cluster
distances will be computed using centroids only. If \code{FALSE},
all observations contained in the clusters will be used.}
  \item{g}{Optional assignment of samples to apriori clusters. By
default, there are no apriori clusters, and clustering starts from
individual observations. If \code{g} is supplied, the clustering
starts from apriori clusters (and the structure within the apriori
clusters depends on the value of the \code{gIntra} argument).
Intuitively, apriori clusters group observations that are known
(from principle) to form compact clusters, whose internal structure
is not of interest from the perspective of the whole hierarchical
clustering. Apriori clusters are typically small: they are much
smaller compared to the clusters whose relative distances
are to be computed using pure Mahalanobis distance. That said, it
makes not much sense to define an apriori cluster of size close to
the number of observations multiplied by \code{thresh} (we raise a
warning in that case).
\code{g} is expected to be a numeric vector of length corresponding
to the number of rows of \code{x}, holding the index of apriori
cluster each sample is member of. \code{0}'s can appear in the
vector meaning that the corresponding sample is not member of any
apriori cluster. Run \code{demo(apriori)} for an example.}
  \item{gIntra}{boolean. If \code{TRUE}, the intrinsic structure of
the apriori clusters (defined using the \code{g} argument) gets
computed before the apriori clusters get merged with their neighbours.
If \code{FALSE},
the intrinsic structure of the apriori clusters is not of interest
(the corresponding elements of the result are not defined and should
be ignored) and clustering starts directly from the level of the
apriori clusters. Run \code{demo(aprioriNonIntra)}' for an example.}
  \item{gParallel}{boolean. If \code{TRUE}, and \code{gIntra} is
also \code{TRUE}, apriori clusters get formed in parallel using the
\code{doMC} package. If \code{FALSE}, apriori clusters get formed
sequentially.}
  \item{normalize}{boolean. If \code{TRUE}, cluster size
will be ignored when computing Mahalanobis distance from the
cluster. If \code{FALSE}, once all clusters are of at least the
\code{thresh} relative size, both cluster shape and size will
affect inter-cluster distance.}
  \item{subthreshHandling}{method
how subthreshold clusters distort the space around them - i.e. how
their size and shape gets reflected in distance computation.
The idea is that subthreshold clusters should not use the
fully-blown Mahalanobis distance, as they are too small to reliably
estimate the covariance matrix for them. However, resigning to pure
Euclidean distance for subthreshold clusters could be too harsh.

The Mahalanobis distance is not used for clusters of two samples
only, as the covariance is degenerated in that case.

The \code{mahal} method pushes the covariance matrices of small
clusters towards a sphere, such that those clusters do not distort
the space around them too much.

The \code{mahal0} method is similar to \code{mahal}, but it
somehow ignores the \code{scale} of the data contained in subthreshold
clusters, which can lead to spurious results for clusters of data
living on too large / small scales. This option is mostly present
only for backward compatibility, and the \code{mahall} method is
considered as areplacement for \code{mahall0}.

The \code{euclidMahal} method enforces spherical (unit) covariances
of all subthreshold clusters (i.e. Euclidean distances are used to
compute distances from them), but clusters above the threshold use
the Mahalanobis distance, eventhough there still could exist
subthreshold clusters using the Euclidean distance. This disparity
could demonstrate in merging large super-threshold elliptical
clusters before subthreshold clusters form another super-threshold
cluster, which could be deemed non-intuitive.

The \code{euclid} method enforces spherical (unit) covariances of all
clusters until there are not any subthreshold clusters. This option
usually leads to formation of compact superthreshold clusters, but
completely ignores the possible intrinsic structure of subthreshold
clusters.

See \code{subthreshX} demos to learn more.}
  \item{warn}{boolean. If \code{FALSE}, warnings about
non-monotonous heights will be suppressed.}
  \item{verb}{level of verbosity, the greater the more detailed
info, defaults to 0 (no info).}
  \item{verbRecursive}{level of verbosity of the recursive calls}
  \item{useR}{if TRUE, R implementation gets used, otherwise, C
implementation gets used. It is not guarranteed that these two
implementations yield exactly same results for the same input data
when the inter-cluster distances are not unique. This is due to
technical differences (the finding of the minimum inter-cluster
distance differs between C and R).}
  \item{nFull}{number of observations; this equals
the number of rows of \code{x} on primary call, but on recursive
calls over a subset of \code{x} it refers to the original full
data set (internal parameter)}
}
\details{This is 'mahalanobis-average' hierarchical clustering similar to
\code{\link[stats]{hclust}} with advanced merging strategy. The shape of clusters
is considered when computing inter-cluster distances.

The distance between two clusters `c1' and `c2' is the mean of
the distances of members of `c1' to the cluster `c2' and the distances
of members of `c2' to the cluster `c1'. The formula is:

dist(c1,c2) = 1/2 * (
mean ( dist ( c1_i, c2 ) ) +
mean ( dist ( c2_i, c1 ) )
),

where `c1',`c2' denote clusters, `c1_i' iterates over members of `c1',
`c2_i' iterates over members of `c2', and

\eqn{dist(x,c) = sqrt ( ( \code{x} - \hat{c} )' * cov(c)^-1 * ( \code{x} - \hat{c} ) )}

where `c' is a cluster, `x' is an observation column vector,
\eqn{\hat{c}} is the center of cluster `c' (mean of the observations
contained in `c').

The distance between an individual observation `x' and a cluster
`c' is some mixture of the Mahalanobis and Euclidean distances from
`c' to `x', depending on the relative cluster size (see the
`thresh' argument) and the method of handling of clusters of
subthreshold size (see the `subthreshHandling' argument).
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

'g': if non-trivial apriori clusters were supplied in the 'g'
argument, this holds the 'g' assignment of samples to apriori
clusters.

'n.height.apriori': if non-trivial apriori clusters were
supplied in the 'g' argument, this component holds the number
of mergings within the apriori clusters. I.e., the first
'n.height.apriori' entries of 'height' and 'merge' describe the
mergings within the apriori clusters, while the subsequent
entries describe the mergings of and above the apriori
clusters.

'g.labels': if non-trivial apriori clusters were
supplied in the 'g' argument, this is a numeric vector of
length 'n.height.apriori', which for each merging within apriori
clusters holds the id of the apriori cluster the merging
appears within.

If apriori clusters were supplied, \code{\link{cutreeApriori}}
can be used to cut off the within-apriori mergings from this
tree object to focus on the mergings of the apriori clusters.}
\references{Fiser K., Sieger T., Schumich A., Wood B., Irving J.,
Mejstrikova E., Dworzak MN. _Detection and monitoring of normal and
leukemic cell populations with hierarchical clustering of flow
cytometry data_. Cytometry A. 2012 Jan;81(1):25-34.
doi:10.1002/cyto.a.21148.
}
\author{Tomas Sieger, Karel Fiser}



\seealso{\code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}},
\code{\link{cutreeApriori}}.
}
\examples{
opar<-par(mfrow=c(2,2))

data(xy)
# the desired number of clusters
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
