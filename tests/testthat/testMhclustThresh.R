context("mhclust: thresh influence")

library(mhca)

test_that("default thresh", {
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)))$height,
        c(0.1,0.2496894,0.4,0.6872512,3.3333333),tolerance=1e-6)
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),useR=TRUE)$height,
        c(0.1,0.2496894,0.4,0.6872512,3.3333333),tolerance=1e-6)
})
test_that("low thresh", {
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),thresh=.1)$height,
        c(0.1,0.25,0.4,0.7,3.333333),tolerance=1e-6)
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),thresh=.1,useR=TRUE)$height,
        c(0.1,0.25,0.4,0.7,3.333333),tolerance=1e-6)
})
test_that("low thresh", {
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),thresh=1)$height,
        c(0.1,0.249922,0.4,0.696585,3.248385),tolerance=1e-6)
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),thresh=1,useR=TRUE)$height,
        c(0.1,0.249922,0.4,0.696585,3.248385),tolerance=1e-6)
})
