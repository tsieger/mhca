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
        # # of observations
        n<-length(h$height)+1L
        if (verb) printWithName(n)
        # # of apriori merge operations
        m<-h$n.height.apriori
        if (verb) printWithName(m)

        # compute the 'merge' after the cut:
        # cut off "m" entries from the bottom
        h$merge<-h$merge[-(1:m),,drop=FALSE]
        if (verb>1) printWithName(h$merge)

        # Note: some clusters are made of apriori clusters, others from
        # singletons (i.e. apriori clusters of only a single
        # observation), and others from single observations not part of
        # any apriori cluster.

        # find indices of clusters (entries in 'h$merge') made from
        # formed apriori clusters (these entries will be positive in
        # the resulting 'merge':
        idx<-h$merge>m
        if (verb>1) printWithName(idx)
        # shift the ids of these clusters to start at the index of 1
        h$merge[idx]<-h$merge[idx]-m
        if (verb>1) printWithName(h$merge)
        # the rest of the 'merge' entries will be negative - they will
        # denote apriori clusters as elementry observations on the
        # resulting 'merge'

        if (verb>2) printWithName(h$g)
        if (verb>2) printWithName(h$g.labels)
        # merge entries that correspond to observations in the 'merge'
        # being computed
        i<-c(h$merge[!idx])
        if (verb>1) printWithName(i)
        # observations in 'i' correspond to apriori clusters,
        # so generate "labels of observations" in the cut tree by
        # taking the labels of observations within the apriori
        # clusters; the labels of formed apriori clusters come in
        # 'g.labels'), labels of singleton can be recalled from 'g'.
        # labels of apriori clusters (incl. singleton apriori
        # clusters):
        id<-rep(0L,length(i))
        # assign labels of merged non-singleton apriori clusters (which
        # are represented as positive numbers), these come in 'g.labels'
        id[i>0]<-h$g.labels[i[i>0]]
        # elementary observations and singleton apriori clusters are
        # represented as negative numbers - refer to the original
        # assignment of observations to apriori clusters in 'g'
        id[i<0]<-h$g[-i[i<0]]
        if (verb>1) .pn(id)
        # sort the labels to correspond to 'h$merge'
        if (verb>1) .pn(order(i))
        h$labels<-paste(id)[order(i)]
        if (verb>1) printWithName(h$labels)

        # populate the new 'merge' with numbers -1..-k referring to
        # "new observations" corresponding to merged apriori clusters
        # (incl. singleton apriori clusters) and also to some original
        # elementary observations not part of any apriori cluster
        h$merge[!idx]<--rank(i)

        # compute the cut 'height':
        h$height<-h$height[-(1:m)]

        # compute the cut 'order':
        ordering<-computeLeafOrder(h$merge)
        # the inverse permutation gives the assignment of observations
        # to the leafs in the dendrogram
        h$order<-order(ordering)

        # remove entries specific to apriori clusters, which are now
        # not contained in the tree:
        h$n.height.apriori<-NULL
        h$g.labels<-NULL
    }
    if (verb>1) printWithName(h$merge)
    return(h)
    ### An object of class *hclust*. The object corresponds to the
    ### original tree \code{h} having the non-apriori portion removed.
    ### See \code{\link{mhclust}} \code{\link{hclust}} for details.
},ex=function() {
    # demo data to cluster
    x<-c(1,2,4,10,13,20,24)
    x<-cbind(x,x)
    h<-mhclust(x,g=c(1,1,1,2,2,3,3))
    h2<-cutreeApriori(h)
    opar<-par(mfrow=c(1,2))
    plot(h)
    plot(h2)
    par(opar)
})
