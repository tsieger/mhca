mhclust<-structure(function # mhca
### Hierarchical clustering with Mahalanobis inter-cluster distances.
###
##details<<
##   This is 'mahalanobis-average' hierarchical clustering similar to
##   \code{\link{hclust}} with advanced merging strategy. The shape of clusters
##   is considered when computing inter-cluster distances.
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
##seealso<< \code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}},
## \code{\link{cutreeApriori}}.
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
g = NULL, ##<< Optional assignment of samples to apriori clusters. By
## default, there are no apriori clusters, and clustering starts from
## individual observations. If \code{g} is supplied, the clustering
## process depends on the value of the \code{gIntra} argument.
## \code{g} is expected to be a numeric vector of length corresponding
## to the number of rows of \code{x}, holding the index of apriori
## cluster each sample is member of. \code{0}'s can appear in the
## vector meaning that the corresponding sample is not member of any
## apriori cluster. Run \code{demo(apriori)} for an example.
gIntra = TRUE, ##<< boolean. If \code{TRUE}, the intrinsic structure of
## the apriori clusters (defined using the \code{g} argument) gets
## computed before the apriori clusters get clustered. If \code{FALSE},
## the intrinsic structure of the apriori clusters is not of interest
## (the corresponding elements of the result are not defined and should
## be ignored) and clustering starts directly from the level of the
## apriori clusters. Run \code{demo(aprioriNonIntra)}' for an example.
gParallel = TRUE, ##<< boolean. If \code{TRUE}, and \code{gIntra} is
## also \code{TRUE}, apriori clusters get formed in parallel using the
## \code{doMC} package. If \code{FALSE}, apriori clusters get formed
## sequentially.
normalize = FALSE, ##<< boolean. If \code{TRUE}, cluster size
## will be ignored when computing Mahalanobis distance from the
## cluster. If \code{FALSE}, once all clusters are of at least the
## \code{thresh} relative size, both cluster shape and size will
## affect inter-cluster distance.
verb = 0, ##<< level of verbosity, the greater the more detailed
## info, defaults to 0 (no info).
verbRecursive = verb, ##<< level of verbosity of the recursive calls
useR = FALSE, ##<< if TRUE, R implementation gets used, otherwise, C
## implementation gets used. It is not guarranteed that these two
## implementations yield exactly same results for the same input data
## when the inter-cluster distances are not unique. This is due to
## technical differences (the finding of the minimum inter-cluster
## distance differs between C and R).
nFull = nrow(as.matrix(x)) ##<< number of observations; this equals
## the number of rows of \code{x} on primary call, but on recursive
## calls over a subset of \code{x} it refers to the original full
## data set (internal parameter)
) {

    # convert to matrix
    if (!is.matrix(x)) x<-as.matrix(x)

    # scale coordinates if necessary
    if (scale) x<-scale(x)

    n<-nrow(x)
    if (verb>1) {
        printWithName(n)
        printWithName(nFull)
        printWithName(thresh)
    }
    if (n==1) stop('please, supply at least 2 observations')
    if (ncol(x)==1) stop('please, supply at least 2-dimensional data')
    if (!is.null(g)) {
        # convert factor to numeric
        if (!is.numeric(g)) stop('g must be numeric')
        g<-as.integer(g)
        if (length(g)!=n) stop('g must be a vector of ',n,' elements')
        if (min(g)<0 || max(g)>n) stop('g must consist of indices of apriori clusters in the range of 1 to ',n)
        # compute the cardinality of each non zero member of g
        gt<-table(g[g>0])
        if (verb>1) printWithName(gt)
        # check for pathological cases
        if (length(gt)==0 || # no non zero entries in g
            gt[1]==n || # a single apriori cluster of all observations
            all(gt==1)) { # each observation forms its own apriori cluster
            g<-NULL
        }
    }
    if (!is.null(g)) {
        # there are some non-trivial apriori clusters
        # Example:
        #     g = 1,1,1,2,2,2,3,0
        #     gt = [1:] 3, [2:] 3  (plus [3:] 1, which gets removed below)
        #     n = 8
        #     nApriori = 6
        #     gMergingCount = 4
        #     nLeft = 4
        #     IDs of subclusters in apirori clusters: 7,8,9,10
        #     clusterIds = 9,10 (n+gMergingCount-length(gt))
        # restrict to clusters of at least 2 observations
        gt<-gt[gt>1]
        # number of observations in the apriori clusters
        nApriori<-sum(gt)
        # number of mergings in apriori clusters
        nHeightApriori<-gMergingCount<-sum(gt-1L)
        # number of clusters left to be merged after apriori clusters have been merged
        nLeft<-n-gMergingCount
        gtLevels<-as.integer(names(gt))

        # cluster apriori clusters first, then merge the clusterings,
        # and cluster the apriori clusters subsequently

        # identity matrix - a surrogate for a correlation matrix in degenerated cases
        fakeInvCov<-diag(ncol(x))
        # vector of IDs of apriori clusters, corresponding to accumulated 'gh' and 'gm'
        gg<-integer(gMergingCount)
        # 'gg' counting index
        gIdx<-0L
        # indices of the last observation in each apriori cluster in 'gg'
        ggLastIndices<-integer(length(gt))
        # 'ggLastIndices' counting index 
        ggLastIndicesIdx<-0L
        # 'height' accumulated over recursive calls
        gh<-integer(gMergingCount)
        # 'merge' accumulated over recursive calls
        gm<-matrix(0L,gMergingCount,2)
        # "g labels" of observations accumulated over recursive calls
        # i.e. the ids of group each observation is member of
        gl<-integer(gMergingCount)
        # members of apriori clusters
        members<-vector('list',nLeft)
        # centroids of apriori clusters
        centroid<-matrix(NA,nLeft,ncol(x))
        # weight factors of apriori clusters
        weightFactor<-numeric(nLeft)
        # square roots of determinants of the covariance matrices of the apriori clusters
        detsSqrt<-numeric(nLeft)
        # covariance matrices of the apriori clusters
        invcov<-vector('list',nLeft)

        # step 1: perform recursive MHCA over apriori subclusters
        #
        # iterate over apriori clusters
        if (gIntra && gParallel && require(doMC)) {
            if (verb) cat(paste0('parallelizing MHCA\n'))
            registerDoMC()
            mhs<-foreach(gti=1:length(gt)) %dopar% {
                gi<-gtLevels[gti]
                i<-which(g==gi)
                iLen<-length(i)
                iLen1<-iLen-1L
                xx<-x[i,]
                mh<-mhclust(xx,thresh=thresh,scale=FALSE,quick=quick,g=NULL,normalize=normalize,verb=verbRecursive,useR=useR,nFull=nFull)
            }
        } else {
            mhs<-NULL
        }
        for (gti in 1:length(gt)) {
            gi<-gtLevels[gti]
            i<-which(g==gi)
            if (verb>2) printWithName(i)
            iLen<-length(i)
            iLen1<-iLen-1L
            xx<-x[i,]
            if (gIntra) {
                if (!is.null(mhs)) {
                    mh<-mhs[[gti]]
                } else {
                    if (verb) cat(paste0('> recursive call (',gti,'/',length(gt),') to cluster ',iLen,' observations\n'))
                    mh<-mhclust(xx,thresh=thresh,scale=FALSE,quick=quick,g=NULL,normalize=normalize,verb=verbRecursive,useR=useR,nFull=nFull)
                    if (verb>1) cat('< returned from the recursive call\n')
                }
            } else {
                # make up a fake, but valid merge structure:
                # -1, -2
                # -3, 1
                # -4, 2
                # etc.
                tmp.merge<-c(-1,-2)
                if (iLen>2) {
                   tmp.merge<-rbind(tmp.merge,cbind(1:(iLen-2),-(3:iLen)))
                }
                tmp.height<-rep(0,iLen-1)
                mh<-list(merge=tmp.merge,height=tmp.height)
            }
            if (verb>2) printWithName(mh)
            if (verb>2) printWithName(mh$height)
            if (verb>2) printWithName(mh$merge)

            # accumulate ID of the current apriori cluster
            gg[gIdx+(1:iLen1)]<-gi
            # accumulate the index of the last member of the current apriori cluster
            ggLastIndices[ggLastIndicesIdx+1]<-gIdx+iLen1
            ggLastIndicesIdx<-ggLastIndicesIdx+1

            members[[gti]]<-i
            centroid[gti,]<-colMeans(xx)
            if (thresh > 0) {
                wf1<-min(1,iLen/(nFull*thresh))
            } else {
                wf1<-0
            }
            weightFactor[gti]<-wf1
            covXx<-cov(xx)
            covXx<-wf1*covXx+(1-wf1)*fakeInvCov
            c.cholDecomp<-tryCatch(chol(covXx),error=function(e)fakeInvCov)
            detsSqrt[gti]<-(1/prod(diag(c.cholDecomp)))^(2/ncol(x))
            invcov[[gti]]<-chol2inv(c.cholDecomp)

            # accumulate the 'height' from the recursive call
            gh[gIdx+(1:iLen1)]<-mh$height
            # accumulate the 'merge' from the recursive call;
            # transform observation indices in 'gm' from local
            # (specific to each apriori cluster) into global
            mh$merge[mh$merge<0]<--i[-mh$merge[mh$merge<0]]
            #  (the clusters formed will get transformed later)
            gm[gIdx+(1:iLen1),]<-mh$merge
            gl[gIdx+(1:iLen1)]<-gi
            gIdx<-gIdx+iLen1
        }
        tmp<-which(g==0)
        if (length(tmp)>0) {
            for (i in 1:length(tmp)) {
                members[length(gt)+i]<-tmp[i]
                centroid[length(gt)+i,]<-x[tmp[i],,drop=FALSE]
                weightFactor[length(gt)+i]<-0.0
                detsSqrt[length(gt)+i]<-1.0
                invcov[[length(gt)+i]]<-fakeInvCov
            }
        }
        clusterSize<-sapply(members,length)
        if (verb>2) printWithName(gg)
        if (verb>2) printWithName(gh)
        if (verb>2) printWithName(gm)

        # step 2: merge the results of MHCA over apriori subclusters
        if (verb) cat('merging the results of MHCA over apriori subclusters\n')
        # we can't simply sort by height, as MHCA is non-monotonic, in general,
        # and naive sorting would produce incorrect dendrogram, in which not-yet-created
        # branches get used before they are formed; instead, we must merge
        # the individual dendrograms, respecting the order of branches in each dendrogram
        gho<-orderStrata(gh,gg)
        # 'gho' determines how to order 'gh' to come sorted increasingly
        if (verb>2) printWithName(gho)
        # 'ghoRank' gives the order of the subclusters of the apriori clusters
        ghoRank<-integer(length(gho))
        ghoRank[gho]<-seq(along=gho)
        if (verb>2) printWithName(ghoRank)
        # 'ggLastIndices' hold the indices of the last member of each apriori cluster,
        # 'gho[ggLastIndices]' then hold the indices of the apriori clusters
        # as formed by HCA, sorted according the their heights
        clusterId<-ghoRank[ggLastIndices]+n
        if (verb>2) printWithName(ggLastIndices)

        # the second pass through apriori clusters
        if (verb>2) cat('merging pre-pre tx\n')
        if (verb>2) printWithName(gm[gho,])
        gIdx<-0L
        for (gti in 1:length(gt)) {
            if (verb>1) cat(paste0(gti,'/',length(gt),'\n'))
            if (verb>2) printWithName(gti)
            gi<-gtLevels[gti]
            i<-which(g==gi)
            if (verb>2) printWithName(i)
            iLen1<-length(i)-1L
            # transform cluster indices in 'gm' from local (specific
            # to each apriori cluster) into global
            tmp<-gm[gIdx+(1:iLen1),]
            mapping<-ghoRank[gg==gi]
            if (verb>2) printWithName(mapping)
            tmp[tmp>0]<-mapping[tmp[tmp>0]]
            gm[gIdx+(1:iLen1),]<-tmp
            gIdx<-gIdx+iLen1
        }

        # sort the accumulated data according to the increasing height of the HCA structure in each apriori cluster
        height<-gh[gho]
        merging<-gm[gho,]
        gLabels<-gl[gho]
        # convert merging to the internal style: all indices are positive
        if (verb>2) cat('merging pre tx\n')
        if (verb>2) printWithName(merging)
        merging[merging>0]<-merging[merging>0]+nApriori
        merging[merging<0]<--merging[merging<0]
        if (verb>2) cat('merging post tx\n')
        if (verb>2) printWithName(merging)
        if (verb>1) cat('merged data\n')
        if (verb>1) printWithName(clusterId)
        if (verb>1) printWithName(clusterSize)
        if (verb>1) printWithName(members)
        if (verb>1) printWithName(centroid)
        if (verb>2) printWithName(weightFactor)
        if (verb>2) printWithName(detsSqrt)
        if (verb>2) printWithName(invcov)
        if (verb>1) printWithName(merging)
        if (verb>1) printWithName(height)
        # initialize distance matrix
        if (verb) cat('computing distances between apriori clusters\n')
        distX<-computeMahalDistMat(distX=NULL,nLeft,x,centroid,members,invcov,detsSqrt,normalize,quick,dbg=verb)
        gMergingCount<-0L
    } else {
        # g not specified
        gMergingCount<-0L
        nLeft<-distX<-centroid<-members<-invcov<-detsSqrt<-weightFactor<-clusterId<-clusterSize<-merging<-height<-NULL
    }

    # implementation switch
    if (useR) {
        if (verb) cat('using R implementation\n')
        rv<-mhclust_Rimpl(x,thresh,scale,quick,normalize,g,gMergingCount,verb,nFull,nLeft,distX,centroid,members,invcov,detsSqrt,weightFactor,clusterId,clusterSize,merging,height)
    } else {
        if (verb) cat('using C implementation\n')
        rv<-mhclust_Cimpl(x,thresh,scale,quick,normalize,g,gMergingCount,verb,nFull,nLeft,distX,centroid,members,invcov,detsSqrt,weightFactor,clusterId,clusterSize,merging,height)
    }

    if (verb>2) cat('rv$merge pre tx\n')
    if (verb>2) printWithName(rv$merge)
    if (verb>2) printWithName(rv$height)

    # sort clusters in rows of `merge'
    mn<-pmin(rv$merge[,1],rv$merge[,2])
    mx<-pmax(rv$merge[,1],rv$merge[,2])
    rv$merge[,1]<-mn
    rv$merge[,2]<-mx

    if (any(diff(rv$height)<0) && nFull==n) warning('clusters form non-monotonic tree')

    # follow the style of hclust: observations get represented by negative numbers, clusters by positive numbers
    rv$merge[rv$merge<=n]<- -rv$merge[rv$merge<=n]
    rv$merge[rv$merge>n]<-rv$merge[rv$merge>n]-n
    if (verb>2) cat('rv$merge post tx\n')
    if (verb>2) printWithName(rv$merge)

    # assignment of leafs in dendrogram to original observations
    ordering<-computeLeafOrder(rv$merge)
    # the inverse permutation gives the assignment of original observations to the leafs in the dendrogram
    ordering<-order(ordering)

    tmp<-list(merge=rv$merge,height=rv$height,order=ordering,labels=rownames(x),
        method='mahalanobis-average',dist.method='euclidean')
    if (!is.null(g)) {
        tmp<-c(tmp,list(n.height.apriori=nHeightApriori,g.labels=gLabels))
    }
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
    ###
    ### 'n.height.apriori': if non-trivial apriori clusters were
    ### supplied in the 'g' argument, this component holds the number
    ### of mergings within the apriori clusters. I.e., the first
    ### 'n.height.apriori' entries of 'height' and 'merge' describe the
    ### mergings within the apriori clusters, while the subsequent
    ### entries describe the mergings of and above the apriori
    ### clusters.
    ###
    ### 'g.labels': if non-trivial apriori clusters were
    ### supplied in the 'g' argument, this is a numeric vector of
    ### length 'n.height.apriori', which for each merging within apriori
    ### clusters holds the id of the apriori cluster the merging
    ### appears within.
    ### 
    ### If apriori clusters were supplied, \code{\link{cutreeApriori}}
    ### can be used to cut off the within-apriori mergings from this
    ### tree object to focus on the mergings of the apriori clusters.

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
