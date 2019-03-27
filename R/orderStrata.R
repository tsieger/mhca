orderStrata<-structure(function # Order a vector, preserving the ordering in strata.
##details<<
## This function merges subvectors (strata) of a vector and returns the ordering
## of the vector elements. The ordering of elements in the strata is preserved,
## even if they are NOT sorted.
(x, ##<< a vector of values
g, ##<< a vector defining the strata
method = c('log','parallel','naive') ##<< method to use, the 'log' outperfmorms the other two
) {
    method<-match.arg(method)

    n<-length(x)
    if (length(g)!=n) {
        stop('"g" and "x" must be of the same length')
    }
    if (n==0) return(integer(0))
    if (!is.numeric(g)) stop('g must be numeric')

    ug<-sort(unique(g))
    ng<-length(ug)
    if (ng==1) {
        return(1:n)
    }
    if (ng==2) {
        i1<-which(g==ug[1])
        i2<-setdiff(1:n,i1)
        o<-orderVectors(x[i1],x[i2])
        o<-c(i1,i2)[o]
        return(o)
    }

    # ng>2
    if (method=='parallel') {
        # Merge individual strata in parallel.
        # This is a heavy-weight solution, slower compared to the 'log'
        # method.
        o<-integer(n)
        indices<-vector('list',ng)
        # vector of the first elements from each stratum
        heads<-numeric(ng)
        for (i in seq(along=ug)) {
            tmp<-which(g==ug[i])
            indices[[i]]<-tmp
            heads[i]<-x[tmp[1]]
        }
        # ordering of heads
        headsOrder<-order(heads)
        # sorted heads, the first element is the first to come in the output vector
        headsSorted<-heads[headsOrder]
        for (i in 1:n) {
            ##printWithName(headsSorted)
            # index of stratum whose first element comes first in headsSorted
            ig<-headsOrder[1]
            ##catnl(paste0('found minimum ',headsSorted[1],' in stratum ',ig))
            # indices of this stratum
            tmp<-indices[[ig]]
            # output the first index from this stratum
            o[i]<-tmp[1]
            ###printWithName(o[i])
            if (length(tmp)>1) {
                indices[[ig]]<-tmp[-1]
                ##printWithName(indices[[ig]])
                # merge the second element of indices[headsOrder[1]] into headsSorted
                # element to insert into heads
                e<-x[tmp[2]]
                ###printWithName(e)
                ###printWithName(headsSorted[2])
                if (e<=headsSorted[2]) {
                    # 'e' comes in the first slot
                    ##printWithName(headsSorted)
                    ##catnl(paste0('inserting ',e,' in the pos 1'))
                    headsSorted[1]<-e
                    # we're done
                } else {
                    # find the position after which to insert 'e'
                    p1<-2L
                    p2<-min(3L,length(headsSorted))
                    while (p1!=p2) {
                        ##printWithName(p1)
                        ##printWithName(p2)
                        if (e>headsSorted[p2]) {
                            # 'e' should come after p2
                            # move the p1..p2 window forward and expand it
                            tmp<-p2
                            p2<-min(p2+2*(p2-p1),length(headsSorted))
                            p1<-tmp
                        } else {
                            # e <= heads[p2]
                            # 'e' should come after p1 and before p2
                            p3<-p1+ceiling((p2-p1)/2)
                            if (e>headsSorted[p3]) {
                                p1<-p3
                            } else {
                                p2<-p3-1L
                            }
                        }
                    }
                    # if inserting before an element with the same value, the
                    # stratum with lower ID takes precedence
                    if (p1<length(headsSorted) && e==headsSorted[p1+1]) {
                        if (headsOrder[p1+1]<ig) {
                            p1<-p1-1L
                        }
                    }
                    # insert 'e' after p1
                    ##printWithName(headsSorted)
                    ##catnl(paste0('inserting ',e,' after element ',headsSorted[p1],' (pos ',p1,')'))
                    if (p1<length(headsSorted)) {
                        headsSorted<-c(headsSorted[2:p1],e,headsSorted[(p1+1):length(headsSorted)])
                        headsOrder<-c(headsOrder[2:p1],ig,headsOrder[(p1+1):length(headsOrder)])
                    } else {
                        headsSorted<-c(headsSorted[2:p1],e)
                        headsOrder<-c(headsOrder[2:p1],ig)
                    }
                }
            } else {
                # this stratum is exhausted
                ##catnl(paste0('stratum ',ig,' exhausted'))
                headsOrder<-headsOrder[-1]
                headsSorted<-headsSorted[-1]
                if (length(headsSorted)==1) {
                    ig<-headsOrder[1]
                    ##catnl(paste0('single stratum ',ig,' remains'))
                    o[(i+1):n]<-indices[[ig]]
                    break
                }
            }
        }
        return(o)
    } else if (method=='log') {
        # Merge the 1st stratum with the 2nd one, the 3rd stratum with
        # the 4th one, etc., and iteratively merge the merged strata in
        # subsequently.
        # This performs well across different choices of the number of
        # strata and the number of observations in them.
        indices<-vector('list',ng)
        for (i in seq(along=ug)) {
            indices[[i]]<-which(g==ug[i])
        }
        while (ng>1) {
            newIndices<-vector('list',ceiling(ng/2))
            i<-1L
            # merge indices consequtive lists of indices
            while (2L*i<=ng) {
                o<-orderVectors(x[indices[[2L*i-1L]]],x[indices[[2L*i]]])
                newIndices[[i]]<-c(indices[[2L*i-1L]],indices[[2L*i]])[o]
                i<-i+1L
            }
            # append the last list, if any
            if (2L*(i-1L)+1L==ng) {
                newIndices[[i]]<-indices[[ng]]
            }
            indices<-newIndices
            ng<-as.integer(ceiling(ng/2))
        }
        return(newIndices[[1]])
    } else if (method=='naive') {
        # Sequentially merge stratum #1 with stratum #2, then
        # merge strat #1+#2 with stratum #3, etc. This takes
        # O(ng*n) time and is too slow for large numbers of strata.
        i1<-which(g==ug[1])
        x1<-x[i1]
        n1<-length(i1)
        ig<-2L
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
                ig<-ig+1L
            }
        }
        return(i1)
    } else {
        stop('unknown method "',method,'"')
    }

},ex=function() {
    orderStrata(c(),c())
    orderStrata(1,1)
    orderStrata(c(1,3,5,2,4,6),c(1,1,1,2,2,2))
    orderStrata(c(1,3,5,4,6,2),c(1,1,1,2,2,2))
    orderStrata(c(1,2,3,4,5,6),c(1,2,3,1,2,3))
    orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2))
    orderStrata(c(1,1,3,4,2,2,4,3,5),c(1,2,1,2,1,2,1,2,2))
})
