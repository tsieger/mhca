context("mhclust: large apriori clusters get detected")

test_that("large apriori clusters get detected", {
    x<-cbind(1:10,1:10)
    expect_warning(
        mhclust(cbind(1:10,1:10),g=c(rep(1,6),rep(2,4)),thresh=.5))
    expect_silent(
        mhclust(cbind(1:10,1:10),g=c(rep(1,6),rep(2,4)),thresh=.6))
})
