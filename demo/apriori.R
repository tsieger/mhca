## Mahalanobis HCA with apriori clusters demo.
##

# prepare demo data
set.seed(1)
x0n<-100
# a unit uniform circle
tmpRadius<-sqrt(runif(x0n))
tmpAngle<-runif(x0n)*2*pi
x0<-data.frame(x=tmpRadius*cos(tmpAngle),y=tmpRadius*sin(tmpAngle))
# shift
a<-1.5
# scaling
b<-8
# left vertical ellipse
x1<-x0
x1$x<-x1$x/b-a
x1$y<-x1$y+a/2
# horizontal ellipse
x2<-x0
x2$x<-x2$x
x2$y<-x2$y/b+a/2
# upper vertical ellipse
x3<-x0
x3$x<-x3$x/b
x3$y<-x3$y+a
# lower vertical ellipse
x4<-x0
x4$x<-x4$x/b
x4$y<-x4$y-a
# combine the data
x<-rbind(x1,x2,x3,x4)
rownames(x)<-1:nrow(x)
n<-nrow(x)

# number of top-most clusters to take
k<-4

# MHCA (w/out apriori clusters)
mh<-mhclust(x,thresh=.2)
cmh<-cutree(mh,k=k)
# MHCA with apriori clusters
mhg<-mhclust(x,thresh=.25,g=rep(c(1:4),c(nrow(x1),nrow(x2),nrow(x3),nrow(x4))))
# note the "thresh=.25" argument: well-formed clusters are expected to be formed
# by 25% of observations
cmhg<-cutree(mhg,k=k)
# heights of and above the apriori clusters:
mhg$height[-(1:mhg$n.height.apriori)]

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
plot(x[,1],x[,2],asp=1,col=cmhg,main='MHCA with 4 apriori clusters')
# dendrogram
plot(mhg,labels=FALSE,main='MHCA with 4 apriori clusters')
y<-min(mhg$height)-diff(range(mhg$height))/10
text(1:n,y,(1:n)[mhg$order],col=cmhg[mhg$order],srt=90)
par(opar)
