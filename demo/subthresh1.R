# Demo on treatment of clusters with subthreshold size: three elongated clusters.

if (exists('saveMe')&&saveMe) png('subthresh1.png',1400,1000)
opar<-par(mfcol=c(3,8),mar=c(3,2,2,1)+0.1)

xy2<-xy
# desired number of clusters
k<-3

n<-nrow(xy2)

useR<-FALSE
#useR<-TRUE

# Mahalanobis HCA
thresh<-1/3
   mh1.vm<-mhclust(  1*xy2,thresh=thresh,useR=useR,subthreshHandling='mahal')
  mh10.vm<-mhclust( 10*xy2,thresh=thresh,useR=useR,subthreshHandling='mahal')
 mh100.vm<-mhclust(100*xy2,thresh=thresh,useR=useR,subthreshHandling='mahal')
  mh1.det<-mhclust(  1*xy2,thresh=thresh,useR=useR,subthreshHandling='mahal0')
 mh10.det<-mhclust( 10*xy2,thresh=thresh,useR=useR,subthreshHandling='mahal0')
mh100.det<-mhclust(100*xy2,thresh=thresh,useR=useR,subthreshHandling='mahal0')
  mh1.eu<-mhclust(  1*xy2,thresh=thresh,useR=useR,subthreshHandling='euclid')
 mh10.eu<-mhclust( 10*xy2,thresh=thresh,useR=useR,subthreshHandling='euclid')
mh100.eu<-mhclust(100*xy2,thresh=thresh,useR=useR,subthreshHandling='euclid')
  mh1.em<-mhclust(  1*xy2,thresh=thresh,useR=useR,subthreshHandling='euclidMahal')
 mh10.em<-mhclust( 10*xy2,thresh=thresh,useR=useR,subthreshHandling='euclidMahal')
mh100.em<-mhclust(100*xy2,thresh=thresh,useR=useR,subthreshHandling='euclidMahal')

  cmh1.det<-cutree(  mh1.det,k=k)
 cmh10.det<-cutree( mh10.det,k=k)
cmh100.det<-cutree(mh100.det,k=k)
  cmh1.vm<-cutree(  mh1.vm,k=k)
 cmh10.vm<-cutree( mh10.vm,k=k)
cmh100.vm<-cutree(mh100.vm,k=k)
  cmh1.eu<-cutree(  mh1.eu,k=k)
 cmh10.eu<-cutree( mh10.eu,k=k)
cmh100.eu<-cutree(mh100.eu,k=k)
  cmh1.em<-cutree(  mh1.em,k=k)
 cmh10.em<-cutree( mh10.em,k=k)
cmh100.em<-cutree(mh100.em,k=k)

plotDendro<-function(h,ch,main) {
    plot(h,labels=FALSE,main=main)
    y<-min(h$height)-diff(range(h$height))/10
    text(1:n,y,(1:n)[h$order],col=ch[h$order],srt=90)
}
addMark<-function(m) {
    mtext(paste0('1',m),line=.3,adj=0,font=2,cex=1.5)
}

#
# subthreshHandling='mahal'
#
# feature space plots with 3 top clusters
plot(  1*xy2[,1],  1*xy2[,2],asp=1,col=  cmh1.vm,main=  '1*xy');addMark('A')
plot( 10*xy2[,1], 10*xy2[,2],asp=1,col= cmh10.vm,main= '10*xy');addMark('B')
plot(100*xy2[,1],100*xy2[,2],asp=1,col=cmh100.vm,main='100*xy');addMark('C')
# MHCA dendrograms
plotDendro(  mh1.vm,  cmh1.vm,main='subthreshHandling=\'mahal\'')
plotDendro( mh10.vm, cmh10.vm,main='subthreshHandling=\'mahal\'')
plotDendro(mh100.vm,cmh100.vm,main='subthreshHandling=\'mahal\'')

#
# subthreshHandling='mahal0
#
# feature space plots with 3 top clusters
plot(  1*xy2[,1],  1*xy2[,2],asp=1,col=  cmh1.det,main=  '1*xy');addMark('D')
plot( 10*xy2[,1], 10*xy2[,2],asp=1,col= cmh10.det,main= '10*xy');addMark('E')
plot(100*xy2[,1],100*xy2[,2],asp=1,col=cmh100.det,main='100*xy');addMark('F')
# MHCA dendrograms
plotDendro(  mh1.det,  cmh1.det,main='subthreshHandling=\'mahal0\'')
plotDendro( mh10.det, cmh10.det,main='subthreshHandling=\'mahal0\'')
plotDendro(mh100.det,cmh100.det,main='subthreshHandling=\'mahal0\'')

#
# subthreshHandling='euclid'
#
# feature space plots with 3 top clusters
plot(  1*xy2[,1],  1*xy2[,2],asp=1,col=  cmh1.eu,main=  '1*xy');addMark('G')
plot( 10*xy2[,1], 10*xy2[,2],asp=1,col= cmh10.eu,main= '10*xy');addMark('H')
plot(100*xy2[,1],100*xy2[,2],asp=1,col=cmh100.eu,main='100*xy');addMark('I')
# MHCA dendrograms
plotDendro(  mh1.eu,  cmh1.eu,main='subthreshHandling=\'euclid\'')
plotDendro( mh10.eu, cmh10.eu,main='subthreshHandling=\'euclid\'')
plotDendro(mh100.eu,cmh100.eu,main='subthreshHandling=\'euclid\'')

#
# subthreshHandling='euclidMahal'
#
# feature space plots with 3 top clusters
plot(  1*xy2[,1],  1*xy2[,2],asp=1,col=  cmh1.em,main=  '1*xy');addMark('J')
plot( 10*xy2[,1], 10*xy2[,2],asp=1,col= cmh10.em,main= '10*xy');addMark('K')
plot(100*xy2[,1],100*xy2[,2],asp=1,col=cmh100.em,main='100*xy');addMark('L')
# MHCA dendrograms
plotDendro(  mh1.em,  cmh1.em,main='subthreshHandling=\'euclidMahal\'')
plotDendro( mh10.em, cmh10.em,main='subthreshHandling=\'euclidMahal\'')
plotDendro(mh100.em,cmh100.em,main='subthreshHandling=\'euclidMahal\'')

par(opar)
if (exists('saveMe')&&saveMe) dev.off()
