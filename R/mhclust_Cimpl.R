#' @useDynLib mhca mhclust_
mhclust_c<-function(X, DistX, Merging, Height, Threshold, Quick, Normalize, Verb) {
    .Call(mhclust_, X, DistX, Merging, Height, Threshold, Quick, Normalize, Verb)
}

# x data matrix
mhclust_Cimpl<-function(x, thresh, scale, quick, normalize, verb) {
    n<-nrow(x)
    d<-dist(x)
    merge<-matrix(0L,n-1,2)
    height<-rep(0,n-1)

    mhclust_c(x,d,merge,height,thresh,quick,normalize,as.integer(verb))

    return(list(merge=merge,height=height))
}
