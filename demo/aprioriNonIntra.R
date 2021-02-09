## Mahalanobis HCA with apriori clusters demo showing how to start
## clustering from the level of the apriori clusters, ignoring the
## intrinsic structure of the apriori clusters.

# cut apriori mergings off the tree using 'cutreeApriori' and show the
# resulting restricted tree?
SHOW_RESTRICTED_TREE<-FALSE

# prepare demo data
set.seed(1)
n<-20
# a unit uniform circle
tmpRadius<-sqrt(runif(n))
tmpAngle<-runif(n)*2*pi
x0<-data.frame(x=tmpRadius*cos(tmpAngle),y=tmpRadius*sin(tmpAngle))
# left circle
x1<-x0
x1$x<-x1$x-2
# right circle
x2<-x0
x2$x<-x2$x+2
# upper circle
x3<-x0
x3$y<-x3$y+2
# combine the data
x<-rbind(x1,x2,x3)
rownames(x)<-1:nrow(x)
n<-nrow(x)

# number of top-most clusters to take
k<-3

# MHCA with apriori clusters, computing the internal structure of the
# apriori clusters
mhgIn<-mhclust(x,thresh=1/3,g=rep(c(1:3),c(nrow(x1),nrow(x2),nrow(x3))))
# note the "thresh=.33" argument: well-formed clusters are expected to be formed
# by 33% of observations
cmhgIn<-cutree(mhgIn,k=k)
# heights of and above the apriori clusters:
mhgIn$height[-(1:mhgIn$n.height.apriori)]

# MHCA with apriori clusters, ignoring the internal structure of the
# apriori clusters
mhgEx<-mhclust(x,thresh=1/3,g=rep(c(1:3),c(nrow(x1),nrow(x2),nrow(x3))),gIntra=FALSE)
# note the "thresh=.33" argument: well-formed clusters are expected to be formed
# by 33% of observations
cmhgEx<-cutree(mhgEx,k=k)
# heights of and above the apriori clusters:
mhgEx$height[-(1:mhgEx$n.height.apriori)]

# make plots
opar<-par(mfcol=c(2+SHOW_RESTRICTED_TREE,2))
# MHCA with apriori clusters, computing the internal structure of the
# apriori clusters
# feature space plots with 2 top clusters
plot(x[,1],x[,2],asp=1,col=cmhgIn,main='Mahalanobis HCA (MHCA)')
if (SHOW_RESTRICTED_TREE) {
    text(mean(x1[,1]),mean(x1[,2]),'1',cex=3,col=scales::alpha('blue',.3))
    text(mean(x2[,1]),mean(x2[,2]),'2',cex=3,col=scales::alpha('blue',.3))
    text(mean(x3[,1]),mean(x3[,2]),'3',cex=3,col=scales::alpha('blue',.3))
}
# dendrogram
plot(mhgIn,labels=FALSE,main='Structure of apriori clusters computed')
y<-min(mhgIn$height)-diff(range(mhgIn$height))/10
text(1:n,y,(1:n)[mhgIn$order],col=cmhgIn[mhgIn$order],srt=90)
if (SHOW_RESTRICTED_TREE) {
    mhgInCut<-cutreeApriori(mhgIn)
    plot(mhgInCut,main='Apriori structure cut off')
}
# MHCA with apriori clusters, ignoring the internal structure of the
# apriori clusters
plot(x[,1],x[,2],asp=1,col=cmhgEx,main='Mahalanobis HCA (MHCA)')
if (SHOW_RESTRICTED_TREE) {
    text(mean(x1[,1]),mean(x1[,2]),'1',cex=3,col=scales::alpha('blue',.3))
    text(mean(x2[,1]),mean(x2[,2]),'2',cex=3,col=scales::alpha('blue',.3))
    text(mean(x3[,1]),mean(x3[,2]),'3',cex=3,col=scales::alpha('blue',.3))
}
# dendrogram
plot(mhgEx,labels=FALSE,main='Structure of apriori clusters ignorable')
y<-min(mhgEx$height)-diff(range(mhgEx$height))/10
text(1:n,y,(1:n)[mhgEx$order],col=cmhgEx[mhgEx$order],srt=90)
if (SHOW_RESTRICTED_TREE) {
    mhgExCut<-cutreeApriori(mhgEx)
    plot(mhgExCut,main='Apriori structure cut off')
}
par(opar)


#######################
# timing on larger data
#######################
# number of total samples
n<-1000
# number of apriori groups
ng<-10
y<-cbind(rnorm(n),rnorm(n),rnorm(n),rnorm(n),rnorm(n))
g<-rep(1:ng,length.out=n)
# computing the internal structure of apriori clusters:
system.time(mhcgIn<-mhclust(y,thresh=1/ng,g=g,gIntra=TRUE))
# ignoring the internal structure of apriori clusters:
system.time(mhcgEx<-mhclust(y,thresh=1/ng,g=g,gIntra=FALSE))
