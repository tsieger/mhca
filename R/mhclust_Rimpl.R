mhclust_Rimpl<-structure(function # mhca
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
##   `c' is some mixture of the Mahalanobis and Euclidean distances from
##   `c' to `x', depending on the relative cluster size (see the
##   `thresh' argument) and the method of handling of clusters of
##   subthreshold size (see the `subthreshHandling' argument).
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
thresh, ##<< real number in the interval of (0,1) defining
## the minimal relative size of cluster (relative to the number of
## observations) whose distance to other clusters will be computed as a
## pure Mahalanobis distance. The distance from smaller clusters will
## be computed as a mixture of Mahalanobis and Euclidean distances,
## with the contribution of the Mahalanobis distance being larger
## for larger clusters, in general. However, the exact handling of
## subthreshold clusters is controlled by the \code{subthreshHandling}
## argument (see below). This threshold balances the uncertainty in the
## cluster shape as estimated by the covariance matrix of its members.
scale, ##<< boolean. Should we transform observations by the
## `\code{\link{scale}}' function?
quick, ##<< boolean. If \code{TRUE}, inter-cluster
## distances will be computed using centroids only. If \code{FALSE},
## all observations contained in the clusters will be used.
normalize, ##<< boolean. If \code{TRUE}, cluster size
## will be ignored when computing Mahalanobis distance from the
## cluster. If \code{FALSE}, once all clusters are of at least the
## \code{thresh} relative size, both cluster shape and size will
## affect inter-cluster distance.
g, ##<< Optional assignment of samples to apriori clusters that
## should get formed (in a hierarchical fashion) before any other
## merging takes place. If NULL, there are no apriori clusters, and
## clustering starts from individual observations. If `g' is not NULL,
## a numeric vector of length corresponding to the number of rows of
## `x' is expected, holding the index of apriori cluster each sample is
## member of. 0's can appear in the vector meaning that the
## corresponding sample is not member of any apriori cluster. See the
## `aprioriMhca' demo for an example.
gMergingCount, ##<< number of merging operations to be performed when
## forming the apriori clusters (this is directly computed from `g' and
## appears here only to eliminate code duplication).
subthreshHandling = c('mahal','mahal0','euclid','euclidMahal'), ##<< method
## how subthreshold clusters distort the space around them - i.e. how
## their size and shape gets reflected in distance computation.
## The idea is that subthreshold clusters should not use the
## fully-blown Mahalanobis distance, as they are too small to reliably
## estimate the covariance matrix for them. However, resigning to pure
## Euclidean distance for subthreshold clusters could be too harsh.
##
## The Mahalanobis distance is not used for clusters of two samples
## only, as the covariance is degenerated in that case.
##
## The \code{mahal} method pushes the covariance matrices of small
## clusters towards a sphere, such that those clusters do not distort
## the space around them too much.
##
## The \code{mahal0} method is similar to \code{mahal}, but it
## somehow ignores the scale of the data contained in subthreshold
## clusters, which can lead to spurious results for clusters of data
## living on too large / small scales. This option is mostly present
## only for backward compatibility, and the \code{mahall} method is
## considered as areplacement for \code{mahall0}.
##
## The \code{euclidMahal} method enforces spherical (unit) covariances
## of all subthreshold clusters (i.e. Euclidean distances are used to
## compute distances from them), but clusters above the threshold use
## the Mahalanobis distance, eventhough there still could exist
## subthreshold clusters using the Euclidean distance. This disparity
## could demonstrate in merging large super-threshold elliptical
## clusters before subthreshold clusters form another super-threshold
## cluster, which could be deemed non-intuitive.
##
## The \code{euclid} method enforces spherical (unit) covariances of all
## clusters until there are not any subthreshold clusters. This option
## usually leads to formation of compact superthreshold clusters, but
## completely ignores the possible intrinsic structure of subthreshold
## clusters.
##
## See \code{subthreshX} demos to learn more.
verb, ##<< level of verbosity, the greater the more detailed
## info, defaults to 0 (no info).
.nFull = NULL, ##<< number of observations; this equals
## the number of rows of \code{x} on primary call, but it refers to
## the original full data set on recursive calls over a subset of
## \code{x} (this is an internal parameter)
.nLeft = NULL, ##<< number of clusters and observations left for
## clustering after apriori clusters have been processed (internal
## parameter used when clustering pre-clustered apriori clusters)
## parameter)
.distX = NULL, ##<< distance matrix (internal parameter used when
## clustering pre-clustered apriori clusters)
.centroid = NULL, ##<< centroid of current clusters/observations to
## cluster (internal parameter used when clustering pre-clustered
## apriori clusters)
.members = NULL, ##<< members of the current clusters/observations to
## cluster (internal parameter used when clustering pre-clustered
## apriori clusters)
.invcov = NULL, ##<< inverse of the covariance matrices of the current
## clusters/observations to cluster (internal parameter used when
## clustering pre-clustered apriori clusters)
.invcovNmf = NULL, ##<< normalization constants of the inverse of the
## covariance matrices (internal parameter used when clustering
## pre-clustered apriori clusters)
.weightFactor = NULL, ##<< weight factors of the current
## clusters/observations (internal parameter used when clustering
## pre-clustered apriori clusters)
.clusterId = NULL, ##<< IDs of the current clusters/observations
## (internal parameter used when clustering pre-clustered apriori
## clusters)
.clusterSize, ##<< the size of the current clusters/observations
## (internal parameter used when clustering pre-clustered apriori
## clusters)
.merging = NULL, ##<< the `merge' matrix representing the clustering of
## the apriori clusters (internal parameter used when clustering
## pre-clustered apriori clusters)
.height = NULL ##<< the `height' vector representing the clustering of
## the apriori clusters (internal parameter used when clustering
## pre-clustered apriori clusters)
) {

    if (verb>2) {
        printWithName(.nFull)
    }
    # sanity check: the optional '.*' arguments starting from '.nLeft' should all be present or all be NULL
    tmp<-c(
        is.null(.nLeft),
        is.null(.distX),
        is.null(.centroid),
        is.null(.members),
        is.null(.invcov),
        is.null(.invcovNmf),
        is.null(.weightFactor),
        is.null(.clusterId),
        is.null(.clusterSize),
        is.null(.merging),
        is.null(.height))
    if (verb>1) printWithName(tmp)
    if (!is.null(.nLeft)) {
        if (verb>2) {
            printWithName(.nLeft)
            printWithName(.centroid)
            printWithName(.members)
            printWithName(.invcov)
            printWithName(.invcovNmf)
            printWithName(.weightFactor)
            printWithName(.clusterId)
            printWithName(.clusterSize)
            printWithName(.merging)
            printWithName(.height)
        }
        # check the consistency of the .* arguments
        stopifnot(all(tmp) || all(!tmp))
        stopifnot(.nLeft==nrow(.centroid))
        stopifnot(.nLeft==length(.members))
        stopifnot(.nLeft==length(.invcov))
        stopifnot(.nLeft==length(.invcovNmf))
        stopifnot(.nLeft==length(.weightFactor))
        stopifnot(.nLeft==length(.clusterId))
        stopifnot(.nLeft==length(.clusterSize))
        stopifnot(nrow(x)-.nLeft==nrow(.merging))
        stopifnot(nrow(x)-.nLeft==length(.height))
    }

    printClusterMembers<-function (m,addParentheses=TRUE) {
        x<-paste(m,collapse=' ')
        if (addParentheses) {
            x<-paste0('(',x,')')
        }
        x
    }
    printMembers<-function (members) {
        for (i in 1:length(members)) {
          cat(paste0('members[',i,']: ',printClusterMembers(members[[i]],addParentheses=FALSE),'\n'))
        }
    }
    printVectorMembers<-function (x) {
        xn<-deparse(substitute(x))
        for (i in 1:length(x)) {
          cat(paste0(xn,'[',i,']: ',x[i],'\n'))
        }
    }

    dbg<-verb

    # number of elementary points subject to clustering (merging)
    pointCount<-nrow(x)
    # number of points in the full (original) data set; on recursive calls over a subset of data
    # pointCount refers to the subset, while fullPointCount refers to the original data
    if (!is.null(.nFull)) {
        fullPointCount<-.nFull
    } else {
        fullPointCount<-pointCount
    }
    if (dbg>1) printWithName(pointCount)
    if (dbg>1) printWithName(fullPointCount)
    # number of dimensions of the feature space
    spaceDim<-ncol(x)
    if (dbg>1) printWithName(spaceDim)
    if (!is.null(.distX)) {
        distX<-.distX
    } else {
        # compute distance matrix
        distX<-dist(x)
    }
    # merging: the output of the clustering - the structure of clusters
    merging<-matrix(0,pointCount-1,3)
    if (!is.null(.merging) && !is.null(.height)) {
        merging[1:(pointCount-.nLeft),1:2]<-.merging
        merging[1:(pointCount-.nLeft),3]<-.height
    }
    # how may clusters (not merged together yet) remain; in the
    # beginning, all points are clusters, in the end only one huge
    # cluster exists
    if (!is.null(.nLeft)) {
        clusterCount<-.nLeft
    } else {
        clusterCount<-pointCount
    }
    # number of elementary points in each cluster
    if (!is.null(.clusterSize)) {
        clusterSize<-.clusterSize
    } else {
        clusterSize<-rep(1,clusterCount)
    }
    # clusters being made (by merging two smaller clusters) are
    # assigned unique IDs, but reside in data structured indexed by
    # index of one of its subclusters - thus we need to map the
    # 1:clusterCount space into IDs of current clusters
    if (!is.null(.clusterId)) {
        clusterId<-.clusterId
    } else {
        clusterId<-1:clusterCount
    }

    # if inverse of covariance matrix can't be computed, use
    # this surrogate
    fakeCov<-fakeInvCov<-diag(spaceDim)

    # members (elementary observations) of each cluster
    if (!is.null(.members)) {
        members<-.members
    } else {
        members<-vector('list',clusterCount)
        members[1:clusterCount]<-1:clusterCount
    }
    # inverse of covariance matrix of members of given cluster
    # (representing the shape of the clusters)
    if (!is.null(.invcov)) {
        invcov<-.invcov
    } else {
        invcov<-rep(list(fakeInvCov),clusterCount)
    }
    # normalizing factor for each cluster (computed from determinant of
    # the inverse of the covariance matrix) making the N-dim volume of
    # clusters equal to 1 if `invcov[[i]]' gets divided by `invcovNmf[i]'
    if (!is.null(.invcovNmf)) {
        invcovNmf<-.invcovNmf
    } else {
        invcovNmf<-rep(1,clusterCount)
    }
    # centroids of cluster
    if (!is.null(.centroid)) {
        centroid<-.centroid
    } else {
        centroid<-x
    }
    # proportional size of each cluster
    if (!is.null(.weightFactor)) {
        weightFactor<-.weightFactor
    } else {
        weightFactor<-rep(0,clusterCount)
    }
    # number of clusters whose relative size is at least
    # thresh
    fullMahalClusterCount<-sum(clusterSize>=thresh*fullPointCount)
    # have all clusters reached the thresh
    # relative size and have we, therefore, switched into "full
    # Mahalanobis" mode?
    switchedToFullMahal<-fullMahalClusterCount==clusterCount

    if (dbg>1) printWithName(thresh)
    if (dbg>1) printWithName(normalize)
    if (dbg>1) printWithName(quick)

    # merge two closest clusters at each step `s'
    for (s in mySeq(pointCount-clusterCount+1,pointCount-1)) {
        if (dbg>2) {
            cat(sprintf('\n====================== step %d ============================\n',s-1))
            printMembers(members)
            printWithName(invcov)
            #printWithName(invcovNmf)
            printVectorMembers(invcovNmf)
            printWithName(centroid)
            printVectorMembers(weightFactor)
            printWithName(clusterCount)
            printVectorMembers(clusterId)
            printVectorMembers(clusterSize)
            printWithName(fullMahalClusterCount)

            cat('distX\n')
            tmp<-matrix(NA,clusterCount,clusterCount)
            tmp[lower.tri(tmp)]<-distX
            print(as.dist(tmp))
        }

        # find two mutually nearest clusters
        k<-which.min(distX)
        v<-distX[k]
        if (dbg>3) cat(sprintf('found minimum %g at %d\n',v,k))

        # (i,j) @ k = n*(i-1) - i*(i-1)/2 + j-i
        # k = -(i^2)/2 + i*(n+1/2-1) -n+j
        # 0 = -(i^2)/2 + i*(n-1/2) -n+j-k
        # 0 = (i^2)/2 - i*(n-1/2) + n+k-j          {1}
        #  now try to substitute "i" for "j" to make the equation {1} to contain "i" only:
        # j = i+1+j0, j0>=0
        # 0 = (i^2)/2 - i*(n-1/2) + n+k-(i+1+j0)
        # 0 = (i^2)/2 - i*(n-1/2+1) + n+k-(1+j0)
        # 0 = (i^2)/2 - i*(n+1/2) + n+k-1-j0       {2}
        # D = (n+1/2)^2 - 4*(n+k-1-j0)/2 = (n-1/2)^2-2*(k-1-j0) >= (n-1/2)^2-2*(k-1) as j0>=0
        #  now solve the equation {2} for "i":
        # i = floor( (n+1/2) - sqrt((n-1/2)^2-2*(k-1)) )
        #  and compute j from {1}:
        # j = (i^2)/2 - i*(n-1/2) + n+k
        #   = k-i*(n-1/2-i/2)+n
        #   = k-i*(n-i/2)+i/2+n
        #   = k-(i-1)*(n-i/2)+i

        # we are merging clusters `i' and `j' into a new one
        i<-as.integer(floor(clusterCount+1/2-sqrt(clusterCount^2-clusterCount+1/4-2*(k-1))))
        j<-as.integer(k-(i-1)*(clusterCount-i/2)+i)
        if (dbg>1) cat(sprintf('c1=%d, c2=%d\n',i,j))

        if (dbg>0) {
            if (dbg>1 || clusterSize[i] + clusterSize[j]<=25) {
                cat(sprintf('Cluster %d: depth %g, merged clusters %d and %d (%s and %s).\n',
                    s+pointCount,v,clusterId[i],clusterId[j],printClusterMembers(members[[i]]),printClusterMembers(members[[j]])))
            } else {
                cat(sprintf('Cluster %d: depth %g, merged clusters %d and %d.\n',
                    s+pointCount,v,clusterId[i],clusterId[j]))
            }
        }

        merging[s,1:3]<-c(clusterId[i],clusterId[j],v)

        # the other clusters just need to be updated in respect to
        # distance to the newly created cluster
        otherClusters1<-as.integer(mySeq(1,i-1))
        otherClusters2<-as.integer(mySeq(i+1,j-1))
        otherClusters3<-as.integer(mySeq(j+1,clusterCount))
        otherClusters<-c(otherClusters1,otherClusters2,otherClusters3)
        if (dbg>2) cat(sprintf("clusterCount %d\n",clusterCount))
        if (dbg>2) printWithName(otherClusters)

        # get all samples constituting the merged clusters
        xij<-x[c(members[[i]],members[[j]]),,drop=FALSE]
      
        # compute the weight factor controlling the Mahalanobis-Euclidean
        # balance (to be applied when measuring distances relatively to
        # the merged cluster)
        if (thresh > 0) {
            wf1<-min( 1, (clusterSize[i] + clusterSize[j]) / ( fullPointCount * thresh ) )
        } else {
            wf1<-0
        }
        # if the cluster being born consists of two samples only, don't try to compute
        # its covariance matrix, but use the spherical covariance matrix instead
        if (clusterSize[i] + clusterSize[j]==2) {
            wf1<-0
            covXij<-fakeCov
        } else {
            # try to fit an ellipsoid to the merged cluster - try to
            # compute the inverse of covariance matrix
            covXij<-cov(xij)
        }
        # update fullMahalClusterCount if necessary
        if (wf1==1) {
            if (weightFactor[i] + weightFactor[j] == 2) fullMahalClusterCount<-fullMahalClusterCount - 1
            else if (weightFactor[i] < 1 && weightFactor[j] < 1) fullMahalClusterCount<-fullMahalClusterCount + 1
        }
        weightFactor[i]<-wf1
        if (dbg>1) printWithName(wf1)

        if (dbg>2) printWithName(xij)
        if (dbg>1) printWithName(covXij)
        # if the cluster consists of few members only, its shape should
        # not be took too serious: the covariance could be made closer
        # to circle, or even replaced by spherical covariance
        if (subthreshHandling=='mahal' || subthreshHandling=='mahal0' ||
            (subthreshHandling=='euclid' || subthreshHandling=='euclidMahal') && wf1==1) { # the condition is verbose, but error-prone
            # shift the covariance towards a sphere

            # determine scaling factor
            if (subthreshHandling=='mahal0') {
                # NOTE: this is NOT scale-independent, i.e. this affects covXij with large elements
                # much less than covXij with small elements!
                mf<-1
                if (dbg>2) cat('mf fixed to 1\n')
            } else {
                # scale the unit covariance matrix by the determinant
                # of the empirical covariance matrix
                mf<-det(covXij)^(1/spaceDim)
                # in case of singular covXij, det is 0 and mf also, but we want to set mf to 1
                if (mf==0) {
                    if (dbg>2) cat('fallback: setting mf to 1\n')
                    mf<-1
                }
            }
            if (dbg>2) printWithName(mf)
            # scale the covariance
            covXij<-wf1 * covXij + (1-wf1) * mf * fakeInvCov
            if (dbg>2) printWithName(covXij)

            # try to compute the Cholesky decomposition of the covariance matrix
            # if it fails, fall back to the unit matrix
            c.cholDecomp<-tryCatch(chol(covXij),error=function(e) fakeInvCov)
            invcov_merged<-chol2inv(c.cholDecomp)
            if (dbg>2) printWithName(invcov_merged)

            # normalization factor is given by the spaceDim-th root of
            # the determinant of the inverse of the covariance matrix,
            # which is the inverse of the determinant of the covariance matrix,
            # which is the square of the product of diagonal entries of its
            # Cholesky decomposition
            icNmf<-(1/prod(diag(c.cholDecomp)))^(2/spaceDim)
            if (dbg>2) printWithName(icNmf)
        } else {
            # force the covariance matrix to be the unit matrix
            invcov_merged<-fakeInvCov
            if (dbg>2) printWithName(invcov_merged)
            icNmf<-1
            if (dbg>2) printWithName(icNmf)
        }

        # compute a new center of the merged cluster
        centroid[i,]<-(clusterSize[i]*centroid[i,,drop=FALSE] + clusterSize[j]*centroid[j,,drop=FALSE])/
            (clusterSize[i] + clusterSize[j])
        if (dbg>1) printWithName(centroid[i,,drop=FALSE])

        # update distX if we haven't reached the point at which we switch to the full
        # Mahalanobis style or if we haven't clustered all the samples in the apriori clusters
        # (otherwise we recompute the whole distX later on)
        if (fullMahalClusterCount < clusterCount-1 || # we haven't reached the point of switch or
          switchedToFullMahal || # we have already switched or
          pointCount-(clusterCount-1) < gMergingCount) { # we've not clustered all the samples in the apriori clusters
          # ( pointCount-(clusterCount-1) is the number of samples clustered so far)
            ic1<-invcov_merged
            icNmf1<-icNmf
            if (dbg>2) printWithName(ic1)
            if (normalize || fullMahalClusterCount < clusterCount-1) { # note: there is no need to check whether
                # we've clustered all the samples in the apriori clusters, as apriori clusters should be small and
                # the "fullMahalClusterCount < clusterCount-1" condition should thus be trivially satisfied when
                # still clustering apriori clusters
                if (dbg>2) cat(sprintf(' normalizing (normalize %d, clusters with full Mahalanobis = %d, clusters =  %d)\n',normalize,fullMahalClusterCount,clusterCount))
                if (fullMahalClusterCount < clusterCount-1 && subthreshHandling=='euclid') {
                    # enforce unit covariance even if this cluster is large enough,
                    # but some of the others are not
                    ic1<-fakeInvCov
                    if (dbg>2) cat(' using fakeInvCov\n')
                } else {
                    ic1<-ic1 / icNmf1
                }
                if (dbg>2) printWithName(ic1)
            }
            for (ii in seq(along=otherClusters)) {
                if (dbg>2) printWithName(ii)
                if (dbg>2) printWithName(otherClusters[ii])

                iRelDistXIdx<-c(otherClusters1*(clusterCount-(otherClusters1+1)/2)-clusterCount+i,
                    i*(clusterCount-(i+1)/2)-clusterCount+otherClusters2,
                    i*(clusterCount-(i+1)/2)-clusterCount+otherClusters3)
                if (dbg>4) {
                    printWithName(iRelDistXIdx)
                    printWithName(ii)
                    printWithName(iRelDistXIdx[ii])
                }

                # if we are still clustering samples of the apriori clusters, the distance
                # between clusters in disctinct apriori clusters are, by definition, set to Inf
                if (!is.null(g) && g[members[[i]][1]]!=g[members[[otherClusters[ii]]][1]] && # clusters in distinct apriori clusters
                     pointCount-(clusterCount-1) < gMergingCount) { # and we are not done with clustering the samples from prior clusters
                    # the distance between clusters in different apriori clusters is defined to be Inf
                    distX[iRelDistXIdx[ii]]<-Inf
                } else {
                    # compute the distance from the newly merged cluster i+j to cluster otherClusters(ii)
                    if (quick) {
                        xc1<-centroid[otherClusters[ii],,drop=FALSE]
                    } else {
                        xc1<-x[members[[otherClusters[ii]]],,drop=FALSE]
                    }
                    if (dbg>3) printWithName(xc1)
                    if (dbg>3) printWithName(centroid[i,,drop=FALSE])
                    xc1<-xc1-matrix(centroid[i,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE)
                    if (dbg>3) printWithName(xc1)
                    # distMaha1 holds a vector of squares of mahalanobis distances:
                    # mean Mahalanobis distance from i+j to some other cluster
                    distMaha1<-mean(sqrt(rowSums((xc1%*%ic1)*xc1)))
                    if (dbg>2) printWithName(distMaha1)

                    # compute the distance from cluster otherClusters(ii) to the newly merged cluster i+j
                    if (quick) {
                        xc2<-centroid[i,,drop=FALSE]
                    } else {
                        xc2<-xij
                    }
                    if (dbg>3) printWithName(xc2)
                    if (dbg>3) printWithName(centroid[otherClusters[ii],,drop=FALSE])
                    xc2<-xc2-matrix(centroid[otherClusters[ii],,drop=FALSE],nrow(xc2),ncol(xc2),byrow=TRUE)
                    if (dbg>3) printWithName(xc2)
                    #if (dbg>3) printWithName(centroid[otherClusters[ii],])
                    #if (dbg>3) printWithName(invcov[[otherClusters[ii]]])
                    ic2<-invcov[[otherClusters[ii]]]
                    if (dbg>3) printWithName(ic2)
                    if (normalize || fullMahalClusterCount < clusterCount-1) {
                        if (fullMahalClusterCount < clusterCount-1 && subthreshHandling=='euclid') {
                            # enforce unit covariance even if this cluster is large enough,
                            # but some of the others are not
                            ic2<-fakeInvCov
                        } else {
                            ic2<-ic2 / invcovNmf[otherClusters[ii]]
                        }
                        if (dbg>3) printWithName(ic2)
                    }
                    # mean Mahalanobis distance
                    distMaha2<-mean(sqrt(rowSums((xc2%*%ic2)*xc2)))
                    if (dbg>2) printWithName(distMaha2)

                    # merge the clusterId(i) <-> clusterId(j) distances
                    distX[iRelDistXIdx[ii]]<-mean(c(distMaha1,distMaha2))
                }
                if (dbg>2) cat(sprintf('Dist from %d=%s to %d=%s: %g.\n',clusterId[otherClusters[ii]],printClusterMembers(members[[otherClusters[ii]]]),
                  s+pointCount,printClusterMembers(c(members[[i]], members[[j]])),distX[iRelDistXIdx[ii]]))
            }
        }
        members[[i]]<-c(members[[i]],members[[j]])
        # clusters clusterId[i] and clusterId[j] merged, remove
        # clusterId[j]-related info, put info about the newly created
        # cluster at position occupied by clusterId[i] previously
        members<-members[-j]
        invcov<-invcov[-j]
        invcov[[i]]<-invcov_merged
        invcovNmf<-invcovNmf[-j]
        invcovNmf[i]<-icNmf
        centroid<-centroid[-j,]
        weightFactor<-weightFactor[-j]
        jRelDistXIdx<-c(otherClusters1*(clusterCount-(otherClusters1+1)/2)-clusterCount+j,
            otherClusters2*(clusterCount-(otherClusters2+1)/2)-clusterCount+j,
            j*(clusterCount-(j+1)/2)-clusterCount+otherClusters3,
            i*(clusterCount-(i+1)/2)-clusterCount+j)
        distX<-distX[-jRelDistXIdx]

        # update clusterCount, clusterSize, clusterId
        clusterCount<-clusterCount-1
        clusterSize[i]<-clusterSize[i]+clusterSize[j]
        clusterSize<-clusterSize[-j]
        clusterId[i]<-pointCount+s
        clusterId<-clusterId[-j]

        # recompute distX in case we reached the point at which we switch to the full Mahalanobis style
        # and we've also clustered all the samples in the apriori clusters
        if (fullMahalClusterCount == clusterCount && # we reached the switch point and
          !switchedToFullMahal && # we haven't swicthed yet and
          gMergingCount <= pointCount-clusterCount) { # we have clustered all samples in the apriori clusters
            switchedToFullMahal<-TRUE
            if (dbg>1) cat('Recomputing all distances.\n')
            idx<-1

            for (i1 in mySeq(1,clusterCount-1)) {
                if (quick) {
                    xc1Orig<-centroid[i1,,drop=FALSE]
                } else {
                    xc1Orig<-x[members[[i1]],,drop=FALSE]
                }
                if (dbg>2) printWithName(xc1Orig)

                ic2<-invcov[[i1]]
                if (dbg>2) printWithName(ic2)
                if (normalize) {
                    ic2<-ic2 / invcovNmf[i1]
                    if (dbg>2) printWithName(ic2)
                }

                for (i2 in mySeq(i1+1,clusterCount)) {
                    if (dbg>1) cat(paste0('recomputing dist from ',i1,' to ',i2,'\n'))

                    # compute the distance from cluster (i1) to (i2)
                    xc1<-xc1Orig-matrix(centroid[i2,,drop=FALSE],nrow(xc1Orig),ncol(xc1Orig),byrow=TRUE)
                    if (dbg>2) printWithName(xc1)
                    ic1<-invcov[[i2]]
                    if (dbg>2) printWithName(ic1)
                    if (normalize) {
                        ic1<-ic1 / invcovNmf[i2]
                        if (dbg>2) printWithName(ic1)
                    }
                    # mean Mahalanobis distance
                    distMaha1<-mean(sqrt(rowSums((xc1%*%ic1)*xc1)))
                    if (dbg>1) printWithName(distMaha1)

                    # compute the distance from cluster (i2) to (i1)
                    if (quick) {
                        xc2<-centroid[i2,,drop=FALSE]
                    } else {
                        xc2<-x[members[[i2]],,drop=FALSE]
                    }
                    if (dbg>2) printWithName(xc2)
                    xc2<-xc2-matrix(centroid[i1,,drop=FALSE],nrow(xc2),ncol(xc2),byrow=TRUE)
                    if (dbg>2) printWithName(xc2)
                    # mean Mahalanobis distance
                    distMaha2<-mean(sqrt(rowSums((xc2%*%ic2)*xc2)))
                    if (dbg>1) printWithName(distMaha2)

                    # merge the clusterId(i) <-> clusterId(j) distances
                    distX[idx]<-mean(c(distMaha1,distMaha2))
                    idx<-idx + 1
                }
            }
        }
    }

    return(list(merge=merging[,1:2,drop=FALSE],height=merging[,3]))
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
},ex=function() {
# TODO
})
