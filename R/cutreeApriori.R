cutreeApriori<-structure(
function # Cut off the non-apriori portion of a tree resulting from
## \code{mhclust}.
##description<<
## This function restricts a tree resulting from \code{mhclust} to the
## portion representing the mergings of the apriori clusters.
(h, ##<< an object of class *hclust*, resulting from \code{mhclust}.
verb=0 ##<< verbosity level
) {
    if (!is.null(h$n.height.apriori) && !is.null(h$g.labels) && h$n.height.apriori>0) {
        n<-length(h$height)+1L
        if (verb) printWithName(n)
        m<-h$n.height.apriori
        if (verb) printWithName(m)

        # compute the cut 'merge':
        # cut off "m" entries from the bottom
        h$merge<-h$merge[-(1:m),,drop=FALSE]
        if (verb>1) printWithName(h$merge)
        # some clusters are made of observations within aprirori
        # clusters, others from from apriori clusters themselves;
        # find indices of clusters made from apriori clusters:
        idx<-h$merge>m
        if (verb>1) printWithName(idx)
        # shift the ids of these clusters to start at the index of 1
        h$merge[idx]<-h$merge[idx]-m
        # observations in the cut tree correspond to apriori clusters,
        # so generate "labels of observations" in the cut tree by
        # taking the ids of observations within the apriori clusters
        # (these ids come in the g.labels), and pick those referred
        # from the "merge" above the 'n.height.apriori' limit:
        if (verb>1) printWithName(h$g.labels)
        if (verb>1) printWithName(t(h$merge)[t(!idx)])
        tmp<-rep(0L,m)
        tmp[rank(h$merge[!idx])]<-h$g.labels[h$merge[!idx]]
        if (verb>1) printWithName(tmp)
        h$labels<-paste(h$g.labels[t(h$merge)[t(!idx)]])
        if (verb>1) printWithName(h$labels)
        h$labels<-paste(tmp[1:(n-m)])
        if (verb>1) printWithName(h$labels)
        # take the apriori clusters as elementary observations: convert
        # to the range of 1 to k and invert them:
        h$merge[!idx]<--rank(h$merge[!idx])

        # compute the cut 'height':
        h$height<-h$height[-(1:m)]

        # compute the cut 'order':
        h$order<-computeLeafOrder(h$merge)

        # remove entries specific to apriori clusters, which are now
        # not contained in the tree:
        h$n.height.apriori<-NULL
        h$g.labels<-NULL
    }
    return(h)
    ### An object of class *hclust*. The object corresponds to the
    ### original tree \code{h} having the non-apriori portion removed.
    ### See \code{\link{mhclust}} \code{\link{hclust}} for details.
},ex=function() {
    # demo data to cluster
    x<-c(1,2,4,10,13,20,24)
    x<-cbind(x,x)
    h<-mhclust(x,g=c(1,1,1,2,2,3,3))
    h<-mhclust(x,g=4-c(1,1,1,2,2,3,3))
    h<-mhclust(x,g=c(1,1,1,2,2,2,2))
    h2<-cutreeApriori(h)
    opar<-par(mfrow=c(1,2))
    plot(h)
    plot(h2)
    par(opar)
})
