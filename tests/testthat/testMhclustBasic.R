context("mhclust: basic functionality")

library(mhca)

test_that("case #1: 4 observations", {
    expect_equal(mhclust(cbind(rep(1,4),c(1,1.1,4,4.5)))$merge,
        rbind(
            c(-1,-2),
            c(-3,-4),
            c(1,2)
            ))
    expect_equal(mhclust(cbind(rep(1,4),c(1,1.1,4,4.5)),useR=TRUE)$merge,
        rbind(
            c(-1,-2),
            c(-3,-4),
            c(1,2)
            ))
    expect_equal(mhclust(cbind(rep(1,4),c(1,1.1,4,4.5)))$height,
        c(0.1,0.5,3.2))
    expect_equal(mhclust(cbind(rep(1,4),c(1,1.1,4,4.5)),useR=TRUE)$height,
        c(0.1,0.5,3.2))
})

test_that("case #2: 6 observations", {
    # case #2
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)))$merge,
        rbind(
            c(-1,-2),
            c(-3, 1),
            c(-5,-6),
            c(-4, 3),
            c( 2, 4)))
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),useR=TRUE)$merge,
        rbind(
            c(-1,-2),
            c(-3, 1),
            c(-5,-6),
            c(-4, 3),
            c( 2, 4)))
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)))$height,
        c(0.1, 0.25, 0.4, 0.7, 10/3),tolerance=1e-6)
    expect_equal(mhclust(cbind(rep(1,6),c(1,1.1,1.3,4,4.5,4.9)),useR=TRUE)$height,
        c(0.1, 0.25, 0.4, 0.7, 10/3),tolerance=1e-6)
})
