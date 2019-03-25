context("mhclust: apriori clusters")

library(mhca)

# prepare test data
set.seed(1)
n<-100
x<-data.frame(x=c(rnorm(n,-1),rnorm(n,1)),y=c(rnorm(n),rnorm(n)))
g<-rep(1:2,each=n)

test_that("R impl.: recusive vs. non-recursive", {
    expect_equal(mhclust(x,g=g,useR=TRUE)$height,
        mhclust(x,g=g)$height,useR=TRUE,gRecursive=FALSE)
})

test_that("C impl.: recusive vs. non-recursive", {
    expect_equal(mhclust(x,g=g)$height,
        mhclust(x,g=g)$height,gRecursive=FALSE)
})

# permute the observations and g and try again
p<-sample(1:(2*n))
x<-x[p,]
g<-g[p]

test_that("R impl.: recusive vs. non-recursive", {
    expect_equal(mhclust(x,g=g,useR=TRUE)$height,
        mhclust(x,g=g)$height,useR=TRUE,gRecursive=FALSE)
})

test_that("C impl.: recusive vs. non-recursive", {
    expect_equal(mhclust(x,g=g)$height,
        mhclust(x,g=g)$height,gRecursive=FALSE)
})
