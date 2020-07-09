fixNonMonotHca<-structure(function
### Solving non-monotonous heights in HCA clustering.
###
##details<<
## \code{\link{fixNonMonotHca}} replaces non-monotonous heights in a
## dendrogram representing a hierarchical clustering with artificial
## values making the heights monotonous.
## However, please note that non-monotonous heights appear naturally in
## dendrograms resulting from some HCA methods including
## \code{link{mhclust}}, such that altering the heights result in
## artificial dendrograms that do NOT represent the structure of the
## clustered data. Please consider using HCA methods producing
## monotonous before attempting to remove the non-monotonous heights,
## e.g. for the purpose of graphical presentation of the dendrogram.
(hca, ##<< an object of class \code{\link[stats]{hclust}}
method=c('eps','halfway'), ##<< replacement method, one of 'eps', and
## 'halfway'.
## The 'eps' method makes the non-monotonous heights to be monotonous
## by shifting them just above the last monotonous height (\code{eps}
## above it). The resulting heights are close to reality,
## eventhough the denrogram is not so visually appealing.
## The 'halfway' method replaces the non-monotonous heights with
## heights lying halfway between neighbouring monotonous heights. This
## gives nicer dendrogram at the expense of making the heights too
## artificial. See examples.
eps = NULL ##<< parameter of the 'eps' method, if \code{NULL}, it
## defaults to half of the minimal positive difference in monotonous
## heights, or \code{1}, if there are no consecutive monotonous
## heights.
) {
    method<-match.arg(method)

    if (method=='eps') {
        h<-hca$height
        n<-length(h)
        dh<-diff(h)

        # compute eps, if not given
        if (is.null(eps)) {
            if (sum(dh>0)>0) {
                eps<-min(dh[dh>0])/2
            } else {
                eps<-1
            }
        }

        if (any(dh<0)) {
            lastH<-0
            for (i in 1:n) {
                if (h[i]<lastH) {
                    lastH<-h[i]<-lastH+eps
                } else {
                    lastH<-h[i]
                }
            }
        }
        hca$height<-h
    } else if (method=='halfway') {
        h<-hca$height
        n<-length(h)
        dh<-diff(h)
        idx<-which(dh<0)
        idx_len<-length(idx)
        idx_i<-1
        while (idx_i <= idx_len) {
            i<-idx[idx_i]+1
            # i is the FIRST "bad" value index in a chunk of "bad" values
            j<-i+1
            while (j<=n && h[j]<h[i-1]) j<-j+1
            # j is AFTER the last "bad" value in a chunk of "bad" values
            if (j>n) {
                if (i>2) {
                    df<-h[i-1]-h[i-2]
                } else {
                    df<-1
                }
            } else {
                df<-(h[j]-h[i-1])/(j-i+1)
            }
            # replace heights by uniform values between the neighbouring "good" values
            h[i:(j-1)]<-h[i-1] + df*(1:(j-i))
            idx_i<-idx_i+1
            while (idx_i <= idx_len && idx[idx_i] <= j-1) idx_i<-idx_i+1
        }
        hca$height<-h
    } else {
        stop('unknown method "',method,'"',sep='')
    }
    return(hca)
    ### An object of class \code{link[stats]{hclust}} having the
    ### \code{height} component fixed.
},ex=function() {
    # simple example
    d<-cbind(1:3,1:3)
    hd<-mhclust(d)
    # original dendrogram
    print(hd$height)
    hdFixed<-fixNonMonotHca(hd)
    # dendrogram with resolved non-monotonous heights
    print(hdFixed$height)

    # another example
    set.seed(1)
    x<-cbind(runif(20),runif(20))
    hx<-hclust(dist(x)^2,'cen')
    hx1<-fixNonMonotHca(hx,method='eps')
    hx2<-fixNonMonotHca(hx,method='halfway')
    opar<-par(mfrow=c(1,3))
    plot(hx,main='Original')
    plot(hx1,main='Fixed by the "eps" method')
    plot(hx2,main='Fixed by the "halfway" method')
    par(opar)
})
