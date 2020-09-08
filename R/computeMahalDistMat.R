# Compute a matrix of between-cluster Mahalanobis distances.
computeMahalDistMat<-function(distX=NULL,clusterCount,x,centroid,members,invcov,invcovNmf,normalize,subthreshHandling,quick,dbg) {
    if (is.null(distX)) {
        distX<-rep(0,1,clusterCount*(clusterCount-1)/2)
    }
    fakeInvCov<-diag(ncol(x))
    idx<-1
    for (i1 in mySeq(1,clusterCount-1)) {
        if (quick) {
            xc1Orig<-centroid[i1,,drop=FALSE]
        } else {
            xc1Orig<-x[members[[i1]],,drop=FALSE]
        }
        if (dbg>2) printWithName(xc1Orig)

        ic2<-invcov[[i1]]
        if (dbg>2) printWithName(ic2)
        if (normalize) {
            if (subthreshHandling=='euclid') {
                ic2<-fakeInvCov
            } else {
                ic2<-ic2 / invcovNmf[i1]
            }
            if (dbg>2) printWithName(ic2)
        }

        for (i2 in mySeq(i1+1,clusterCount)) {
            if (dbg>1) cat(paste0('recomputing dist from ',i1,' to ',i2,'\n'))

            # compute the distance from cluster (i1) to (i2)
            xc1<-xc1Orig-matrix(centroid[i2,,drop=FALSE],nrow(xc1Orig),ncol(xc1Orig),byrow=TRUE)
            if (dbg>2) printWithName(xc1)
            ic1<-invcov[[i2]]
            if (dbg>2) printWithName(ic1)
            if (normalize) {
                if (subthreshHandling=='euclid') {
                    ic1<-fakeInvCov
                } else {
                    ic1<-ic1 / invcovNmf[i2]
                }
                if (dbg>2) printWithName(ic1)
            }
            # mean Mahalanobis distance
            distMaha1<-mean(sqrt(rowSums((xc1%*%ic1)*xc1)))
            if (dbg>1) printWithName(distMaha1)

            # compute the distance from cluster (i2) to (i1)
            if (quick) {
                xc2<-centroid[i2,,drop=FALSE]
            } else {
                xc2<-x[members[[i2]],,drop=FALSE]
            }
            if (dbg>2) printWithName(xc2)
            xc2<-xc2-matrix(centroid[i1,,drop=FALSE],nrow(xc2),ncol(xc2),byrow=TRUE)
            if (dbg>2) printWithName(xc2)
            # mean Mahalanobis distance
            distMaha2<-mean(sqrt(rowSums((xc2%*%ic2)*xc2)))
            if (dbg>1) printWithName(distMaha2)

            # merge the clusterId(i) <-> clusterId(j) distances
            distX[idx]<-mean(c(distMaha1,distMaha2))
            idx<-idx + 1
        }
    }

    return(distX)
}
