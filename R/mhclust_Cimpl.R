#' @useDynLib mhca mhclust_
mhclust_c<-function(X,DistX,Merging,Height,Threshold,Quick,Normalize,G,GMergingCount,SubthreshHandlingId,Verb,
    .NFull,.NLeft,.Centroid,.Members,.Invcov,.DetsSqrt,.WeightFactor,.ClusterId,.ClusterSize,.membersPoolSize){

    .Call(mhclust_,X,DistX,Merging,Height,Threshold,Quick,Normalize,G,GMergingCount,SubthreshHandlingId,Verb,
        .NFull,.NLeft,.Centroid,.Members,.Invcov,.DetsSqrt,.WeightFactor,.ClusterId,.ClusterSize,.membersPoolSize)
}

mhclust_Cimpl<-function(x,thresh,scale,quick,normalize,g,gMergingCount,subthreshHandling,verb,
    .nFull,.nLeft,.distX,.centroid,.members,.invcov,.invcovNmf,.weightFactor,.clusterId,.clusterSize,.merging,.height) {

    SUBTHRESHOLD_METHOD_MAHAL<-0
    SUBTHRESHOLD_METHOD_MAHAL0<-1
    SUBTHRESHOLD_METHOD_EUCLID<-2
    SUBTHRESHOLD_METHOD_EUCLID_MAHAL<-3

    if (subthreshHandling=='mahal') {
        subthreshHandlingId<-SUBTHRESHOLD_METHOD_MAHAL
    } else if (subthreshHandling=='mahal0') {
        subthreshHandlingId<-SUBTHRESHOLD_METHOD_MAHAL0
    } else if (subthreshHandling=='euclid') {
        subthreshHandlingId<-SUBTHRESHOLD_METHOD_EUCLID
    } else if (subthreshHandling=='euclidMahal') {
        subthreshHandlingId<-SUBTHRESHOLD_METHOD_EUCLID_MAHAL
    } else {
        stop(paste0('unsupported subthreshHandling of \'',subthreshHandling,'\''))
    }

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

    n<-nrow(x)
    if (!is.null(.distX)) {
        d<-.distX
    } else {
        if (verb) cat('computing distance matrix\n')
        d<-dist(x)
    }
    merge<-matrix(0L,n-1,2)
    if (!is.null(.merging)) {
        merge[1:(n-.nLeft),]<-as.integer(.merging)
    }
    height<-rep(0,n-1)
    if (!is.null(.height)) {
        height[1:(n-.nLeft)]<-.height
    }

    if (!is.null(.nFull)) .nFull<-as.integer(.nFull)
    if (!is.null(.nLeft)) {
        .nLeft<-as.integer(.nLeft)
        .clusterId<-as.integer(.clusterId)
        .clusterSize<-as.integer(.clusterSize)
        # the maximum size of pool needed to hold all the members of clusters formed during the clustering
        # e.g. in case of four clusters, we need to hold at most:
        # n1, n2, n3, n4, n1+n2, n1+n2+n3, n1+n2+n3+n4,
        # where n1 >= n2 >= n3 >= n4,
        # i.e. sum(Ni) + cumsum(Ni) - n1, where Ni=c(n1,n2,n3,n4)
        .membersPoolSize<-sum(.clusterSize)+sum(cumsum(sort(.clusterSize,decreasing=TRUE)))-max(.clusterSize)
    } else {
        .membersPoolSize<-NULL
    }
    mhclust_c(x,d,merge,height,thresh,quick,normalize,as.integer(g),as.integer(gMergingCount),as.integer(subthreshHandlingId),as.integer(verb),
        .nFull,.nLeft,.centroid,.members,.invcov,.invcovNmf,.weightFactor,.clusterId,.clusterSize,.membersPoolSize)

    return(list(merge=merge,height=height))
}
