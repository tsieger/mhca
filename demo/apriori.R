## Mahalanobis HCA with apriori clusters demo.
##

# prepare demo data
set.seed(1)
x0n<-100
# a unit uniform circle
x0<-data.frame(x=runif(2*x0n,-1,1),y=runif(2*x0n,-1,1))
x0<-x0[x0$x^2+x0$y^2<1,]
x0<-x0[1:min(nrow(x0),x0n),]
# a horizontal ellipse
x1<-x0
x1$x<-x1$x+1.5
x1$y<-x1$y/10
# a vertical ellipse
x2<-x0
x2$x<-x2$x/10
# another vertical  ellipse
x3<-x0
x3$x<-x3$x/10
x3$y<-x3$y-3.5
# random noise
x4<-data.frame(x=rnorm(x0n)-1,y=rnorm(x0n)-1)
# combine the data
x<-rbind(x1,x2,x3,x4)
rownames(x)<-1:nrow(x)
n<-nrow(x)

# number of top-most clusters to take
k<-4

# MHCA
mh<-mhclust(x,thresh=.2)
cmh<-cutree(mh,k=k)
# MHCA with apriori clusters
mhg<-mhclust(x,thresh=.25,g=rep(c(1:3,0),c(nrow(x1),nrow(x2),nrow(x3),nrow(x4))))
# note the "thresh=.25" argument: wel-formed clusters are expected to be formed
# by 25% of observations
cmhg<-cutree(mhg,k=k)

# make plots
opar<-par(mfcol=c(2,2))
# MHCA
# feature space plots with 4 top clusters
plot(x[,1],x[,2],asp=1,col=cmh,main='Mahalanobis HCA (MHCA)')
# dendrogram
plot(mh,labels=FALSE,main='MHCA')
y<-min(mh$height)-diff(range(mh$height))/10
text(1:n,y,(1:n)[mh$order],col=cmh[mh$order],srt=90)
# MHCA with apriori clusters
# feature space plots with 4 top clusters
plot(x[,1],x[,2],asp=1,col=cmhg,main='MHCA with apriori clusters')
# dendrogram
plot(mhg,labels=FALSE,main='MHCA with apriori clusters')
y<-min(mhg$height)-diff(range(mhg$height))/10
text(1:n,y,(1:n)[mhg$order],col=cmhg[mhg$order],srt=90)
par(opar)
