orderVectors<-structure(function # Order elements of two vectors, preserving the within-vector ordering.
##details<<
## This function merges two almost-increasingly-sorted vectors in
## attempt to produce one increasingly sorted vector. If the vectors
## are sorted, the result is sorted perfectly. If, however, the vectors
## are not sorted, the relative ordering of their elements is preserved
## even in the output vector, which is thus not sorted.
## The algorithm simply goes through elements of the input vectors and
## iteratively puts the smaller of the currently considered elements in
## the resulting vector. For example, the vectors of
## \code{x=c(1, 4, 2)} and \code{y=c(3, 5)} get merged to form the
## vector of \code{sxy=c(1, 3, 4, 2, 5)}. The ordering
## \code{o=c(1, 4, 2, 3, 5)} gets returned, which gives
## \code{sxy = c(x,y)[o]}.
(x, ##<< first vector to merge
y ##<< second vector to merge
) {
    nx<-length(x)
    ny<-length(y)
    n<-nx+ny

    if (nx==0) {
        return(seq(along=y))
    }
    if (ny==0) {
        return(seq(along=x))
    }

    rv<-numeric(n)
    ix<-iy<-1
    for (i in 1:n) {
        if (x[ix]<=y[iy]) {
            rv[i]<-ix
            if (ix==nx) {
                rv[(i+1):n]<-nx+(iy:ny)
                break
            }
            ix<-ix+1
        } else {
            rv[i]<-iy+nx
            if (iy==ny) {
                rv[(i+1):n]<-(ix:nx)
                break
            }
            iy<-iy+1
        }
    }
    return(rv)
    ### The order \code{o} such that \code{c(x,y)[o]} forms the merged
    ### vector, which is "almost sorted" (see Details).
},ex=function() {
    orderVectors(c(),c())
    orderVectors(c(),1:3)
    orderVectors(1:3,c())
    orderVectors(1:3,2:4)
    orderVectors(c(1,3,2),c(1.1,1.9,4))
    orderVectors(c(1,4,2),c(3,5))
})
