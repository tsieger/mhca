context("orderStrata")

library(mhca)

test_that("orderStrata", {

    expect_equal(orderStrata(c(),c()),numeric(0))
    expect_equal(orderStrata(1,1),1)
    expect_equal(orderStrata(c(1,3,5, 2,4,6),c(1,1,1, 2,2,2)),c(1,4,2,5,3,6))
    expect_equal(orderStrata(c(1,3,5, 4,6,2),c(1,1,1, 2,2,2)),c(1,2,4,3,5,6))
    expect_equal(orderStrata(c(1,2,3, 4,5,6),c(1,2,3, 1,2,3)),c(1,2,3,4,5,6))
    expect_equal(orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2)),c(1,5,2,3,4,6,7,8,9))
    expect_equal(orderStrata(c(1,1, 3,4),c(1,2, 1,2)),c(1,2,3,4))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2),c(1,2, 1,2, 1,2)),c(1,2,3,5,4,6))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2, 4,3,5),c(1,2, 1,2, 1,2, 1,2,2)),c(1,2,3,5,7,4,6,8,9))
    expect_equal(orderStrata(c(1,1,1,1),c(1,1,2,2)),c(1,2,3,4))
    expect_equal(orderStrata(c(1,1,1,1),c(1,2,1,2)),c(1,3,2,4))

})
