#' @useDynLib mhca mhclust_
mhclust_c<-function(X, DistX, Merging, Height, Threshold, Quick, Normalize, G, GMergingCount, Verb) {
    .Call(mhclust_, X, DistX, Merging, Height, Threshold, Quick, Normalize, G, GMergingCount, Verb)
}

# x data matrix
mhclust_Cimpl<-function(x, thresh, scale, quick, normalize, g, gMergingCount, gDistIdx, verb) {
    n<-nrow(x)
    d<-dist(x)
    if (!is.null(gDistIdx)) {
      # force distances between distinct apriori clusters to be much greater
      # than the largest distance (we can't use Inf for technical reasons)
      d[-gDistIdx]<-1000*max(d)
    }
    merge<-matrix(0L,n-1,2)
    height<-rep(0,n-1)

    mhclust_c(x,d,merge,height,thresh,quick,normalize,as.integer(g),as.integer(gMergingCount),as.integer(verb))

    return(list(merge=merge,height=height))
}
