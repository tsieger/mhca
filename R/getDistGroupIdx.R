getDistGroupIdx<-structure(
function # Get indices into the distance matrix where members of a group reside.
##description<<
## This function computes the indices of entries in the distance matrix
## where distances between members of a given sample set reside.
(n, ##<< number of observations the distance matrix is computed from
idx ##<< indices of observations whose corresponding indices into the
## distance matrix are to be computed
) {
    if (!is.numeric(idx)) stop('group idx must be numeric')
    if (any(idx<1 || idx>n)) stop('group idx must consist of observation indices in the range of 1 to n')
    m<-length(idx)
    if (m<=1) {
        res<-c()
    } else {
        res<-rep(0L,m*(m-1)/2)
        resCnt<-0
        for (i in 1:(m-1)) {
            idxFrom<-idx[i]
            idxTo<-idx[(i+1):m]
            res[resCnt+(1:(m-i))] <- n*(idxFrom-1) - idxFrom*(idxFrom-1)/2 + idxTo-idxFrom
            resCnt<-resCnt+m-i
        }
    }
    res<-as.integer(res)
    return(res)
    ### An integer vector of indices into the distance matrix.

},ex=function() {
    getDistGroupIdx(6,1:3)
    getDistGroupIdx(6,5:6)
})
