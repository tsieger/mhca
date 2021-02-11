## Demo on the effect of the 'quick' argument.
## This demo shows how the choice of the 'quick' argument can change
## the merging process.

# prepare demo data
set.seed(1)
# number of samples in each apriori cluster
n1<-100
# number of apriori clusters
m<-3
# 2D data
x<-matrix(NA,n1*m,2)
x[(1-1)*n1+(1:n1),1]<-rnorm(n1,1,.1)
x[(1-1)*n1+(1:n1),2]<-rnorm(n1,0, 1)

x[(2-1)*n1+(1:n1),1]<-rnorm(n1,2,.1)
x[(2-1)*n1+(1:n1),2]<-rnorm(n1,0, 1)

x[(3-1)*n1+(1:n1),1]<-rnorm(n1,3,.1)
x[(3-1)*n1+(1:n1),2]<-rnorm(n1,0,.1)

# assignment of observation into apriori clusters
apriori<-rep(1:m,each=n1)

plotHca<-function(mh,main) {
    plot(mh,labels=FALSE,main=main)
    y<-min(mh$height)-diff(range(mh$height))/10
    cmh<-cutree(mh,k=m)
    text(1:(n1*m),y,(1:(n1*m))[mh$order],col=cmh[mh$order]+1,srt=90)
}

# MHCA with apriori clusters
mh1<-mhclust(x,thresh=1/m,g=apriori,quick=FALSE)
mh2<-mhclust(x,thresh=1/m,g=apriori,quick=TRUE)

layout(rbind(1:3))
plot(x[,1],x[,2],asp=1,col=apriori+1,main=paste0('Mahalanobis HCA (MHCA) (n=',n1*m,')'))
plotHca(mh1,main='quick=FALSE')
plotHca(mh2,main='quick=TRUE')
layout(1)
