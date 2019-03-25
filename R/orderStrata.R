orderStrata<-structure(function # Order a vector, preserving the ordering in strata.
##details<<
## This function merges subvectors (strata) of a vector and returns the ordering
## of the vector elements. The ordering of elements in the strata is preserved,
## even if they are NOT sorted.
(x, ##<< a vector of values
g ##<< a vector defining the strata
) {
    n<-length(x)
    if (length(g)!=n) {
        stop('"g" and "x" must be of the same length')
    }
    if (n==0) return(numeric(0))

    ug<-unique(g)
    ng<-length(ug)
    if (ng==1) {
        return(1:n)
    }

    # Sequentially merge stratum #1 with stratum #2, then
    # merge stratum #1+2 with stratum #3, etc. This takes
    # O(ng*n) time.
    # TODO: this could be optimized if needed.
    i1<-which(g==ug[1])
    x1<-x[i1]
    n1<-length(i1)
    ig<-2
    repeat {
        i2<-which(g==ug[ig])
        n2<-length(i2)
        x2<-x[i2]
        o<-orderVectors(x1,x2)
        i1<-c(i1,i2)[o]
        if (ig==ng) {
            break
        } else {
            n1<-n1+n2
            x1<-c(x1,x2)[o]
            ig<-ig+1
        }
    }
    return(i1)
},ex=function() {
    orderStrata(c(),c())
    orderStrata(1,1)
    orderStrata(c(1,3,5,2,4,6),c(1,1,1,2,2,2))
    orderStrata(c(1,3,5,4,6,2),c(1,1,1,2,2,2))
    orderStrata(c(1,2,3,4,5,6),c(1,2,3,1,2,3))
    orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2))
    orderStrata(c(1,1,3,4,2,2,4,3,5),c(1,2,1,2,1,2,1,2,2))
})
