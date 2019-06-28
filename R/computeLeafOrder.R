computeLeafOrder<-structure(
function # Order leafs in a HCA dendrogram. 
##description<<
## This function orders leafs in a dendrogram resulting from a
## hierarchical cluster analysis such that the branches in the
## dendrogram do not cross each other.
(merge, ##<< an `n-1` by 2 matrix describing the merges made during
## HCA, usually the `merge` component of the return value from
## \code{\link[mhca]{mhclust}} or \code{\link[stats]{hclust}}.
verb=0 ##<< verbosity level
) {
    
    # m = # of clusters
    m<-dim(merge)[1]
    if (verb) printWithName(m)
    # n = # of leafs
    n<-m+1
    if (verb) printWithName(n)

    # compute 'order'
    lifo<-rep(0,m)
    lifoIdx<-0
    ordering<-rep(0,n)
    # running leaf index being assigned to leafs in order
    oi<-1
    # push the top merge
    if (verb) cat(sprintf('pushing %d\n',m))
    lifoIdx<-lifoIdx+1
    lifo[lifoIdx]<-m
    while (lifoIdx>0) {
        clstr<-lifo[lifoIdx]
        lifoIdx<-lifoIdx-1
        if (verb) cat(sprintf('popping %d\n',clstr))
        if (clstr<0) {
            # leaf found -> assign it an index
            if (verb) cat(sprintf('  -> %d\n',oi))
            ordering[-clstr]<-oi
            oi<-oi+1
        } else {
            # process both subclusters
            for (i in 2:1) {
                subclstr<-merge[clstr,i]
                if (verb) cat(sprintf('pushing %d\n',subclstr))
                lifoIdx<-lifoIdx+1
                lifo[lifoIdx]<-subclstr
            }
        }
    }
    if (verb) printWithName(ordering)
    return(ordering)
    ### A vector of leafs ordering.
},ex=function() {
# internal function
})
