fixNonMonotHca<-structure(function
### Solving non-monotonous heights in HCA clustering.
###
##details<<
## \code{\link{fixNonMonotHca}} removes non-monotonous heights in a
## dendrogram representing a hierarchical clustering. he non-monotonous
## heights get replaced with artificial values lying halfway between
## neighbouring monotonous heights.
## However, please note that non-monotonous heights appear naturally in
## dendrograms resulting from some HCA methods including
## \code{link{mhclust}}, such that altering the heights result in
## artificial dendrograms that do NOT represent the structure of the
## clustered data. Please consider using HCA methods producing
## monotonous before attempting to remove the non-monotonous heights,
## e.g. for the purpose of graphical presentation of the dendrogram.
(hca ##<< an object of class \code{\link[stats]{hclust}}
) {
    h<-hca$height
    n<-length(h)
    dh<-diff(h)
    idx<-which(dh<=0)
    idx_len<-length(idx)
    idx_i<-1
    while (idx_i <= idx_len) {
        i<-idx[idx_i]+1
        # i is the FIRST bad value index
        j<-i+1
        while (j<=n && h[j]<=h[i-1]) j<-j+1
        # j is AFTER the last bad value
        if (j>n) {
            if (i>2) {
                df<-h[i-1]-h[i-2]
            } else {
                df<-1
            }
        } else {
            df<-(h[j]-h[i-1])/(j-i+1)
        }
        h[i:(j-1)]<-h[i-1] + df*(1:(j-i))
        idx_i<-idx_i+1
        while (idx_i <= idx_len && idx[idx_i] <= j-1) idx_i<-idx_i+1
    }
    hca$height<-h
    return(hca)
    ### An object of class \code{link[stats]{hclust}} having the
    ### \code{height} component fixed.
},ex=function() {
    d<-cbind(1:3,1:3)
    hd<-mhclust(d)
    # original dendrogram
    print(hd$height)

    hdFixed<-fixNonMonotHca(hd)
    # dendrogram with resolved non-monotonous heights
    print(hdFixed$height)
})
