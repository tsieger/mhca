mhclust<-structure(function # mhca
### Hierarchical clustering with Mahalanobis inter-cluster distances.
###
##details<<
##   This is 'mahalanobis-average' hierarchical clustering similar to
##   `hclust' with advanced merging strategy. The shape of clusters is
##   considered when computing inter-cluster distances.
##
##   The distance between two clusters `c1' and `c2' is the mean of
##   the distances of members of `c1' to the cluster `c2' and the distances
##   of members of `c2' to the cluster `c1'. The formula is:
##
##      dist(c1,c2) = 1/2 * (
##                      mean ( dist ( c1_i, c2 ) ) +
##                      mean ( dist ( c2_i, c1 ) )
##                    ),
##
##   where `c1',`c2' denote clusters, `c1_i' iterates over members of `c1',
##   `c2_i' iterates over members of `c2', and
##
##      \eqn{dist(x,c) = sqrt ( ( x - \hat{c} )' * cov(c)^-1 * ( x - \hat{c} ) )}
##
##   where `c' is a cluster, `x' is an observation column vector,
##   \eqn{\hat{c}} is the center of cluster `c' (mean of the observations
##   contained in `c').
##
##   The distance between an individual observation `x' and a cluster
##   `c' is a mixture of the Mahalanobis and Euclidean distances from
##   `c' to `x', weighted by the relative cluster size (see the
##   `thresh' parameter).
##
##references<< Fiser K., Sieger T., Schumich A., Wood B., Irving J.,
## Mejstrikova E., Dworzak MN. _Detection and monitoring of normal and
## leukemic cell populations with hierarchical clustering of flow
## cytometry data_. Cytometry A. 2012 Jan;81(1):25-34.
## doi:10.1002/cyto.a.21148.
##
##seealso<< \code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}}.
##
(x, ##<< data.frame or matrix containing coordinates of elementary
## observations. Rows correspond to observations, columns correspond
## to dimensions of the feature space.
thresh = 0.5, ##<< real number in the interval of (0,1) defining
## the minimal relative size of cluster (relative to the number of
## observations) whose distance to other clusters will be computed as a
## pure Mahalanobis distance. The distance from smaller clusters will
## be computed as a mixture of Mahalanobis and Euclidean distances,
## with the contribution of the Mahalanobis distance being proportional
## to the cluster size. This threshold balances the uncertainty in the
## cluster shape as estimated by the covariance matrix of its members.
scale = FALSE, ##<< boolean. Should we transform observations by the
## `\code{\link{scale}}' function?
quick = FALSE, ##<< boolean. If \code{TRUE}, inter-cluster
## distances will be computed using centroids only. If \code{FALSE},
## all observations contained in the clusters will be used.
g = NULL, ##<< Optional assignment of samples to apriori clusters that
## should get formed (in a hierarchical fashion) before any other
## merging takes place. By default, there are no apriori clusters, and
## clustering starts from individual observations. If `g' is supplied,
## a numeric vector of length corresponding to the number of rows of
## `x' is expected, holding the index of apriori cluster each sample is
## member of. 0's can appear in the vector meaning that the
## corresponding sample is not member of any apriori cluster. See the
## `apriori' demo for an example.
normalize = FALSE, ##<< boolean. If \code{TRUE}, cluster size
## will be ignored when computing Mahalanobis distance from the
## cluster. If \code{FALSE}, once all clusters are of at least the
## \code{thresh} relative size, both cluster shape and size will
## affect inter-cluster distance.
verb = 0, ##<< level of verbosity, the greater the more detailed
## info, defaults to 0 (no info).
useR = FALSE ##<< if TRUE, R implementation gets used, otherwise, C
## implementation gets used
) {

    # convert to matrix
    if (!is.matrix(x)) x<-as.matrix(x)

    n<-nrow(x)
    if (n==1) stop('please, supply at least 2 observations')
    if (ncol(x)==1) stop('please, supply at least 2-dimensional data')
    if (!is.null(g)) {
        # convert factor to numeric
        if (!is.numeric(g)) stop('g must be numeric')
        g<-as.integer(g)
        if (min(g)<0 || max(g)>n) stop('g must consist of indices of apriori clusters in the range of 1 to n')
        # compute the cardinality of each member of g
        gt<-table(g)
        if (verb>1) printWithName(gt)
        gMergingCount<-sum(gt[gt>1]-1)
        gtLevels<-as.numeric(names(gt))
        # indices into the distance matrix where pairs of points g reside
        gDistIdx<-rep(0L,n*(n-1)/2)
        gDistIdxCnt<-0
        for (gti in 1:length(gt)) {
            i<-which(g==gtLevels[gti])
            if (length(i)==n) {
                # there is only one huge apriori cluster - no need to do anything special
                g<-NULL
            } else if (length(i)>1) {
                tmp<-getDistGroupIdx(n,i)
                gDistIdx[gDistIdxCnt+(1:length(tmp))]<-tmp
                gDistIdxCnt<-gDistIdxCnt+length(tmp)
            } else {
                # set the apriori cluster to 0 for apriori clusters of only 1 sample
                g[i]<-0
            }
        }
        gDistIdx<-gDistIdx[1:gDistIdxCnt]
    } else {
      gMergingCount<-0L
      gDistIdx<-NULL
    }

    # scale coordinates if necessary
    if (scale) x<-scale(x)

    # implementation switch
    if (useR) {
        if (verb) cat('using R implementation\n')
        rv<-mhclust_Rimpl(x,thresh,scale,quick,normalize,g,gMergingCount,gDistIdx,verb)
    } else {
        if (verb) cat('using C implementation\n')
        rv<-mhclust_Cimpl(x,thresh,scale,quick,normalize,g,gMergingCount,gDistIdx,verb)
    }

    # sort clusters in rows of `merge'
    mn<-pmin(rv$merge[,1],rv$merge[,2])
    mx<-pmax(rv$merge[,1],rv$merge[,2])
    rv$merge[,1]<-mn
    rv$merge[,2]<-mx

    if (any(diff(rv$height)<0)) warning('clusters form non-monotonic tree')

    # follow the style of hclust: observations get represented by negative numbers, clusters by positive numbers
    rv$merge[rv$merge<=n]<- -rv$merge[rv$merge<=n]
    rv$merge[rv$merge>n]<-rv$merge[rv$merge>n]-n

    # assignment of leafs in dendrogram to original observations
    ordering<-computeLeafOrder(rv$merge)
    # the inverse permutation gives the assignment of original observations to the leafs in the dendrogram
    ordering<-order(ordering)

    tmp=list(merge=rv$merge,height=rv$height,order=ordering,labels=rownames(x),
        method='mahalanobis-average',dist.method='euclidean')
    class(tmp)<-'hclust'

    return(tmp)
    ### An object of class *hclust*. The object is a list with
    ### components:
    ###
    ### 'merge': an n-1 by 2 matrix. The `i'-th row describes the two
    ### clusters merged at the `i'-th step of the clustering. If an
    ### element `j' is negative, then observation `-j' was merged at this
    ### stage.  If `j' is positive, the merge was with the cluster
    ### formed at the (earlier) stage `j' of the algorithm. Negative
    ### entries thus indicate agglomerations of singletons, and
    ### positive entries indicate agglomerations of non-singletons.
    ###
    ### 'height': a set of `n'-1 non-decreasing real values, the
    ### clustering heights - the distances between the clusters merged
    ### at each step.
    ###
    ### 'order': a vector giving the permutation of the original
    ### observations suitable for plotting, in the sense that a cluster
    ### plot using this ordering and matrix 'merge' will not have
    ### crossings of the branches.
    ###
    ### 'labels': labels of each of the objects being clustered.
    ###
    ### 'method': 'mahalanobis-average'
    ###
    ### 'dist.method': 'euclidean'
    ### A list with components:
    ###
    ### 'merge': an n-1 by 2 matrix. The `i'-th row describes the two
    ### clusters merged at the `i'-th step of the clustering. If an
    ### element `j' is negative, then observation `-j' was merged at this
    ### stage.  If `j' is positive, the merge was with the cluster
    ### formed at the (earlier) stage `j' of the algorithm. Negative
    ### entries thus indicate agglomerations of singletons, and
    ### positive entries indicate agglomerations of non-singletons.
    ###
    ### 'height': a set of `n'-1 non-decreasing real values, the
    ### clustering heights - the distances between the clusters merged
    ### at each step.
    ###
    ### 'order': a vector giving the permutation of the original
    ### observations suitable for plotting, in the sense that a cluster
    ### plot using this ordering and matrix 'merge' will not have
    ### crossings of the branches.
    ###
    ### 'labels': labels of each of the objects being clustered.
    ###
    ### 'method': 'mahalanobis-average'
    ###
    ### 'dist.method': 'euclidean'
},ex=function() {
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
})
