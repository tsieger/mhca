context("mhclust: cutreeApriori")

library(mhca)

test_that("3 full apriori clusters", {
    x<-c(1,2,4,10,13,20,24)
    x<-cbind(x,x)
    h<-mhclust(x,g=c(1,1,1,2,2,3,3))
    expect_equal(h$merge,
        rbind(
            c(-1,-2), # 1 <- g.labels
            c(-3, 1), # 1
            c(-4,-5), # 2
            c(-6,-7), # 3
            c( 2, 3),
            c( 4, 5)
        )
    )

    h2<-cutreeApriori(h)
    expect_equal(h2$merge,
        rbind(
            c(-1,-2),
            c(-3, 1)
        )
    )

    expect_equal(h2$labels,c('1','2','3'))
})

test_that("2 apriori clusters incl. 1 singleton", {
    h<-mhclust(data.frame(a=c(1,2,3), b=c(2,3,4)), g=c(1,3,3), thresh=.99)
    expect_equal(h$merge,
        rbind(
            c(-2,-3),
            c(-1, 1)
        )
    )

    h2<-cutreeApriori(h)
    expect_equal(h2$merge,
        rbind(c(-1,-2))
    )

    expect_equal(h2$labels,c('1','3'))
})

test_that("3 apriori clusters incl. 2 singletons", {
    h<-mhclust(data.frame(a=c(1,2,3,4), b=c(1,2,3,4)), g=c(1,3,3,2))
    expect_equal(h$merge,
        rbind(
            c(-2,-3),
            c(-1, 1),
            c(-4, 2)
        )
    )

    h2<-cutreeApriori(h)
    expect_equal(h2$merge,
        rbind(
            c(-2,-3),
            c(-1,1))
    )

    expect_equal(h2$labels,c('2','1','3'))
})

