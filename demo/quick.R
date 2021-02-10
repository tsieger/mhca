## Demo on the effect of the 'quick' argument.
## The 'quick' argument controls whether all observations are
## considered when computing a distance to cluster (quick=FALSE, the
## default), or whether centroids only get considered (quick=TRUE).
## The choice of quick=TRUE leads to a slight improvement in runtime,
## but the result can be compromised.

# prepare demo data
set.seed(1)
# number of samples in each apriori cluster
n1<-100
# number of apriori clusters
m<-6
# 2D data
x<-matrix(NA,n1*m,2)
x[(1-1)*n1+(1:n1),1]<-rnorm(n1,-8,1)
x[(1-1)*n1+(1:n1),2]<-rnorm(n1, 1,.2)
x[(2-1)*n1+(1:n1),1]<-rnorm(n1,-8,.2)
x[(2-1)*n1+(1:n1),2]<-rnorm(n1,-1.5,1)

x[(3-1)*n1+(1:n1),1]<-rnorm(n1, 0,.2)
x[(3-1)*n1+(1:n1),2]<-rnorm(n1, 3,1)
x[(4-1)*n1+(1:n1),1]<-rnorm(n1, 0,.2)
x[(4-1)*n1+(1:n1),2]<-rnorm(n1,-3,1)

x[(5-1)*n1+(1:n1),1]<-rnorm(n1, 8,1)
x[(5-1)*n1+(1:n1),2]<-rnorm(n1, 1,.2)
x[(6-1)*n1+(1:n1),1]<-rnorm(n1, 8,1)
x[(6-1)*n1+(1:n1),2]<-rnorm(n1,-1,.2)

# assignment of observation into apriori clusters
apriori<-rep(1:m,each=n1)

plotHca<-function(mh,main) {
    plot(mh,labels=FALSE,main=main)
    y<-min(mh$height)-diff(range(mh$height))/10
    cmh<-cutree(mh,k=m)
    text(1:(n1*m),y,(1:(n1*m))[mh$order],col=cmh[mh$order]+1,srt=90)
}

# MHCA with apriori clusters
tm1<-system.time(mhg1<-mhclust(x,thresh=1/m,g=apriori)) # note the "thresh=1/m" argument:
# well-formed clusters are expected to be formed by `1/m-th of observations, i.e. by `n1' observation
tm2<-system.time(mhg2<-mhclust(x,thresh=1/m,g=apriori,quick=TRUE))
tm3<-system.time(mhg3<-mhclust(x,thresh=1/m))
tm4<-system.time(mhg4<-mhclust(x,thresh=1/m,quick=TRUE))

layout(rbind(c(0,2,3,4),c(1,2,3,4),c(1,5,6,7),c(0,5,6,7)))
plot(x[,1],x[,2],asp=1,col=apriori+1,main=paste0('Mahalanobis HCA (MHCA) (n=',n1*m,')'))
plotHca(mhg1,main=c('apriori, quick=FALSE',paste0('time ',round(tm1[3],3),'s')))
plotHca(mhg2,main=c('apriori, quick=TRUE',paste0('time ',round(tm2[3],3),'s')))
plotHca(fixNonMonotHca(mhg2),main='apriori, quick=TRUE, non-monot. fixed')
plotHca(mhg3,main=c('no apriori, quick=FALSE',paste0('time ',round(tm3[3],3),'s')))
plotHca(mhg4,main=c('no apriori, quick=TRUE',paste0('time ',round(tm4[3],3),'s')))
plotHca(fixNonMonotHca(mhg4),main='no apriori, quick=TRUE, non-monot. fixed')
layout(1)

# Note that when quick=TRUE, the top-most cluster is weird: the distance
# between the cluster made up from the pair of cyan and blue clusters
# and the other cluster made up from the rest of clusters is small, as it
# got computed from centroids only, which are close to each other.
# The default choice quick=FALSE leads to a more intuitive
# result, as all observations in the clusters get considered in
# distance computations.
#
# Also note that using the fixNonMonotHca() function to make the
# dendrogram more visually appealing is an ad-hoc heuristics, which
# can't be considered as a proper solution to the problem.
