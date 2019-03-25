context("orderVectors")

library(mhca)

test_that("orderVectors", {

    expect_equal(orderVectors(c(),c()),numeric(0))
    expect_equal(orderVectors(c(),1:3),1:3)
    expect_equal(orderVectors(1:3,c()),1:3)
    expect_equal(orderVectors(1:3,2:4),c(1,2,4,3,5,6))
    expect_equal(orderVectors(c(1,3,2),c(1.1,1.9,4)),c(1,4,5,2,3,6))
    expect_equal(orderVectors(c(1,4,2),c(3,5)),c(1,4,2,3,5))

})
