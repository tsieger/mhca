# Test of issue #2: bad covariance matrix for superthreshold clusters with subthreshHandling set to 'euclid'

opar<-par(mfcol=c(2,2),mar=c(3,2,2,1)+0.1)

xy<-rbind(
    # cluster 1
    cbind(-10.5, 0),
    cbind(-10.5, .1),
    cbind(-10.5, .2),
    cbind(-10, 0),
    cbind(-10, .1),
    cbind(-10, .2),
    cbind(-9.5, 0),
    cbind(-9.5,.1),
    cbind(-9.5,.2),

    # cluster 2
    cbind(-.5,0),
    cbind(-.5,.1),
    cbind(-.5,.2),
    cbind(0,0),
    cbind(0,.1),
    cbind(0,.2),
    cbind(.5,0),
    cbind(.5,.1),
    cbind(.5,.2),

    # cluster 3
    cbind(8,-.5),
    cbind(8.1,-.5),
    cbind(8.2,-.5),
    cbind(8,0),
    cbind(8.1,0),
    cbind(8.2,0),
    cbind(8,.5),
    cbind(8.1,.5),
    cbind(8.2,.5))
# each of the 3 top-level clusters consist of 3 superthreshold clusters

useR<-FALSE
useR<-TRUE

# Mahalanobis HCA with normalization off
thresh<-3/27 # 3 points form a full clusters
mh.normOff<-mhclust(xy,thresh=thresh,useR=useR,subthreshHandling='euclid',normalize=FALSE)
mh.normOn<-mhclust(xy,thresh=thresh,useR=useR,subthreshHandling='euclid',normalize=TRUE)
cmh.normOff<-cutree(mh.normOff,k=2)
cmh.normOn <-cutree(mh.normOn, k=2)

plotDendro<-function(h,ch,main) {
    n<-length(h$height)+1
    plot(h,labels=FALSE,main=main)
    y<-min(h$height)-diff(range(h$height))/10
    text(1:n,y,(1:n)[h$order],col=ch[h$order],srt=90)
}
plot(xy[,1],xy[,2],asp=1,col=cmh.normOff,main=c('normalize=FALSE','the two left-most clusters should come together'));
plotDendro(mh.normOff,cmh.normOff,main='normalize=FALSE')
plot(xy[,1],xy[,2],asp=1,col=cmh.normOn,main=c('normalize=TRUE','the two left-most clusters should come together'));
plotDendro(mh.normOn,cmh.normOn,main='normalize=TRUE')

par(opar)

