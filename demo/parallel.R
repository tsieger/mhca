## Mahalanobis HCA with parallel clustering of apriori clusters.

# prepare demo data
set.seed(1)
# number of samples in each apriori cluster
n1<-1000
# number of apriori clusters
m<-20
x<-matrix(NA,n1*m,2)
for (i in 1:m) {
    x[(i-1)*n1+(1:n1),1]<-rnorm(n1)/runif(1,.2,5)+runif(1,-15,15)
    x[(i-1)*n1+(1:n1),2]<-rnorm(n1)/runif(1,.2,5)+runif(1,-15,15)
}
apriori<-rep(1:m,each=n1)

# disabling parallel procesing:
system.time(mhgSeq<-mhclust(x,thresh=1/m,g=apriori,gParallel=FALSE))

# parallel procesing enabled by default:
system.time(mhgPar<-mhclust(x,thresh=1/m,g=apriori,gParallel=TRUE))
