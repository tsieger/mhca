## Mahalanobis HCA with "big" data (n1=2000) split / not split into apriori clusters.
## Grouping observations in apriori clusters can speed up MHCA considerably.

# prepare demo data
set.seed(1)
# number of samples in each apriori cluster
n1<-100
# number of apriori clusters
m<-20
x<-matrix(NA,n1*m,2)
for (i in 1:m) {
    x[(i-1)*n1+(1:n1),1]<-rnorm(n1)/runif(1,.2,5)+runif(1,-15,15)
    x[(i-1)*n1+(1:n1),2]<-rnorm(n1)/runif(1,.2,5)+runif(1,-15,15)
}
apriori<-rep(1:m,each=n1)

# MHCA (w/out apriori clusters)
st<-system.time(mh<-mhclust(x,thresh=1/m))
# note the "thresh=1/m" argument: wel-formed clusters are expected to be formed
# by `1/m'-th of observations, i.e. by `n1' observation
cmh<-cutree(mhg,k=m)
# MHCA with apriori clusters
stg<-system.time(mhg<-mhclust(x,thresh=1/m,g=apriori))
# note the "thresh=1/m" argument: well-formed clusters are expected to be formed
# by `1/m-th of observations, i.e. by `n1' observation
cmhg<-cutree(mhg,k=m)

# make plots
opar<-par(mfcol=c(2,2))
# MHCA
# feature space plot
plot(x[,1],x[,2],asp=1,col=apriori+1,main=paste0('Mahalanobis HCA (MHCA) (n=',n1*m,')'))
# dendrogram
plot(mh,labels=FALSE,main=paste('MHCA took',round(st[3],2),'s'))
y<-min(mh$height)-diff(range(mh$height))/10
text(1:(n1*m),y,(1:(n1*m))[mh$order],col=cmh[mh$order],srt=90)
# MHCA with apriori clusters
# feature space plot
plot(x[,1],x[,2],asp=1,col=apriori+1,main=paste0('MHCA with apriori clusters (n=',n1*m,')'))
# dendrogram
plot(mhg,labels=FALSE,main=paste('apriori MHCA took',round(stg[3],2),'s'))
y<-min(mhg$height)-diff(range(mhg$height))/10
text(1:(n1*m),y,(1:(n1*m))[mhg$order],col=cmhg[mhg$order],srt=90)
par(opar)
