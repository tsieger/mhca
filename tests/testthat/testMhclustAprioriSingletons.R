context("mhclust: apriori clusters with incl. singletons")
# (test of #1)

library(mhca)

test_that("R impl.: a singleton apriori cluster", {
    expect_equal(
        mhclust(data.frame(a=c(1,2,3), b=c(2,3,4)), g=c(1,3,3),useR=TRUE)$merge,
        rbind(c(-2,-3),c(-1,1))
    )
})

test_that("C impl.: a singleton apriori cluster", {
    expect_equal(
        mhclust(data.frame(a=c(1,2,3), b=c(2,3,4)), g=c(1,3,3),useR=FALSE)$merge,
        rbind(c(-2,-3),c(-1,1))
    )
})

test_that("only singleton apriori clusters", {
    expect_equal(
        mhclust(data.frame(a=c(1,2,3), b=c(2,3,5)), g=c(1,2,3))$merge,
        rbind(c(-1,-2),c(-3,1))
    )
})

test_that("only one non-singleton apriori cluster", {
    expect_equal(
        mhclust(data.frame(a=c(1,2,3,4,5), b=c(2,3,4,5,6)), g=c(1,2,4,4,5))$merge,
        rbind(
            c(-3,-4),
            c(-1,-2),
            c( 1, 2),
            c(-5, 3)
        )
    )
})

