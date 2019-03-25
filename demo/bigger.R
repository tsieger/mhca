## Mahalanobis HCA with 1,000,000 observations in 1,000 apriori clusters.

# prepare demo data
set.seed(1)
# number of samples in each apriori cluster
n1<-1000
# number of apriori clusters
m<-1000
x<-matrix(NA,n1*m,2)
for (i in 1:m) {
    x[(i-1)*n1+(1:n1),1]<-rnorm(n1)/runif(1,.2,5)+runif(1,-15,15)
    x[(i-1)*n1+(1:n1),2]<-rnorm(n1)/runif(1,.2,5)+runif(1,-15,15)
}
apriori<-rep(1:m,each=n1)

# MHCA with apriori clusters
system.time(mhg<-mhclust(x,thresh=1/m,g=apriori))
# note the "thresh=1/m" argument: well-formed clusters are expected to be formed
# by `1/m-th of observations, i.e. by `n1' observation

