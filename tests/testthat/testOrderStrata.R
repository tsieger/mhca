context("orderStrata")

library(mhca)

test_that("orderStrata", {

    expect_equal(orderStrata(c(),c()), numeric(0))
    expect_equal(orderStrata(1,1), 1)
    expect_equal(orderStrata(c(1,3,5, 2,4,6),c(1,1,1, 2,2,2)), c(1,4,2,5,3,6))
    expect_equal(orderStrata(c(1,3,5, 4,6,2),c(1,1,1, 2,2,2)), c(1,2,4,3,5,6))
    expect_equal(orderStrata(c(1,2,3, 4,5,6),c(1,2,3, 1,2,3)), c(1,2,3,4,5,6))
    expect_equal(orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2)), c(1,5,2,3,4,6,7,8,9))
    expect_equal(orderStrata(c(1,1, 3,4),c(1,2, 1,2)), c(1,2,3,4))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2),c(1,2, 1,2, 1,2)), c(1,2,3,5,4,6))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2, 4,3,5),c(1,2, 1,2, 1,2, 1,2,2)), c(1,2,3,5,7,4,6,8,9))
    expect_equal(orderStrata(c(1,1,1,1),c(1,1,2,2)), c(1,2,3,4))
    expect_equal(orderStrata(c(1,1,1,1),c(1,2,1,2)), c(1,3,2,4))
})

test_that("orderStrataImplementations", {

    # compare the different implementations
    expect_equal(orderStrata(c(),c(),method='log'), orderStrata(c(),c(),method='parallel'))
    expect_equal(orderStrata(1,1,method='log'), orderStrata(1,1,method='parallel'))
    expect_equal(orderStrata(c(1,3,5, 2,4,6),c(1,1,1, 2,2,2),method='log'), orderStrata(c(1,3,5, 2,4,6),c(1,1,1, 2,2,2),method='parallel'))
    expect_equal(orderStrata(c(1,3,5, 4,6,2),c(1,1,1, 2,2,2),method='log'), orderStrata(c(1,3,5, 4,6,2),c(1,1,1, 2,2,2),method='parallel'))
    expect_equal(orderStrata(c(1,2,3, 4,5,6),c(1,2,3, 1,2,3),method='log'), orderStrata(c(1,2,3, 4,5,6),c(1,2,3, 1,2,3),method='parallel'))
    expect_equal(orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2),method='log'), orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2),method='parallel'))
    expect_equal(orderStrata(c(1,1, 3,4),c(1,2, 1,2),method='log'), orderStrata(c(1,1, 3,4),c(1,2, 1,2),method='parallel'))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2),c(1,2, 1,2, 1,2),method='log'), orderStrata(c(1,1, 3,4, 2,2),c(1,2, 1,2, 1,2),method='parallel'))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2, 4,3,5),c(1,2, 1,2, 1,2, 1,2,2),method='log'), orderStrata(c(1,1, 3,4, 2,2, 4,3,5),c(1,2, 1,2, 1,2, 1,2,2),method='parallel'))
    expect_equal(orderStrata(c(1,1,1,1),c(1,1,2,2),method='log'), orderStrata(c(1,1,1,1),c(1,1,2,2),method='parallel'))
    expect_equal(orderStrata(c(1,1,1,1),c(1,2,1,2),method='log'), orderStrata(c(1,1,1,1),c(1,2,1,2),method='parallel'))
    #
    expect_equal(orderStrata(c(),c(),method='log'), orderStrata(c(),c(),method='naive'))
    expect_equal(orderStrata(1,1,method='log'), orderStrata(1,1,method='naive'))
    expect_equal(orderStrata(c(1,3,5, 2,4,6),c(1,1,1, 2,2,2),method='log'), orderStrata(c(1,3,5, 2,4,6),c(1,1,1, 2,2,2),method='naive'))
    expect_equal(orderStrata(c(1,3,5, 4,6,2),c(1,1,1, 2,2,2),method='log'), orderStrata(c(1,3,5, 4,6,2),c(1,1,1, 2,2,2),method='naive'))
    expect_equal(orderStrata(c(1,2,3, 4,5,6),c(1,2,3, 1,2,3),method='log'), orderStrata(c(1,2,3, 4,5,6),c(1,2,3, 1,2,3),method='naive'))
    expect_equal(orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2),method='log'), orderStrata(c(1,3,2,4, 1,4,2,3,5),c(1,1,1,1, 2,2,2,2,2),method='naive'))
    expect_equal(orderStrata(c(1,1, 3,4),c(1,2, 1,2),method='log'), orderStrata(c(1,1, 3,4),c(1,2, 1,2),method='naive'))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2),c(1,2, 1,2, 1,2),method='log'), orderStrata(c(1,1, 3,4, 2,2),c(1,2, 1,2, 1,2),method='naive'))
    expect_equal(orderStrata(c(1,1, 3,4, 2,2, 4,3,5),c(1,2, 1,2, 1,2, 1,2,2),method='log'), orderStrata(c(1,1, 3,4, 2,2, 4,3,5),c(1,2, 1,2, 1,2, 1,2,2),method='naive'))
    expect_equal(orderStrata(c(1,1,1,1),c(1,1,2,2),method='log'), orderStrata(c(1,1,1,1),c(1,1,2,2),method='naive'))
    expect_equal(orderStrata(c(1,1,1,1),c(1,2,1,2),method='log'), orderStrata(c(1,1,1,1),c(1,2,1,2),method='naive'))
})

test_that("orderStrata", {

    expect_equal(orderStrata(c(2,6, 6,7),c(1,1, 2,2)), c(1,2,3,4))
    expect_equal(orderStrata(c(2,6, 6,7),c(2,2, 1,1)), c(1,3,2,4))
    expect_equal(orderStrata(c(2,5, 4, 6,7, 8),c(3,3, 1, 2,2, 4)), c(1,3,2,4,5,6))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(1,1,1,1,1, 2,2,2, 3,3,3)), c(1,6,9, 2,7,10, 8,3,4,5,11))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 2,2,2, 1,1,1)), c(1,6,9, 2,7,10, 8,3,4,5,11))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 1,1,1, 2,2,2)), c(1,6,9, 2,7,10, 8,3,4,5,11))
    expect_equal(orderStrata(c(2,4,6,8,10,12,14,16,20,18,17),c(1,2,3,1,1,2,3,2,3,1,3)), c(1,2,3,4,5,6,7,8,10,9,11))
    expect_equal(orderStrata(c(1,10,2:16),c(1,1,2:16)), c(1,3:10,2,11:17))
})

test_that("orderStrataImplementations2", {

    # compare the different implementations
    expect_equal(orderStrata(c(2,6, 6,7),c(1,1, 2,2),method='log'), orderStrata(c(2,6, 6,7),c(1,1, 2,2),method='parallel'))
    expect_equal(orderStrata(c(2,6, 6,7),c(2,2, 1,1),method='log'), orderStrata(c(2,6, 6,7),c(2,2, 1,1),method='parallel'))
    expect_equal(orderStrata(c(2,5, 4, 6,7, 8),c(3,3, 1, 2,2, 4),method='log'), orderStrata(c(2,5, 4, 6,7, 8),c(3,3, 1, 2,2, 4),method='parallel'))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(1,1,1,1,1, 2,2,2, 3,3,3),method='log'), orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(1,1,1,1,1, 2,2,2, 3,3,3),method='parallel'))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 2,2,2, 1,1,1),method='log'), orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 2,2,2, 1,1,1),method='parallel'))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 1,1,1, 2,2,2),method='log'), orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 1,1,1, 2,2,2),method='parallel'))
    expect_equal(orderStrata(c(2,4,6,8,10,12,14,16,20,18,17),c(1,2,3,1,1,2,3,2,3,1,3),method='log'), orderStrata(c(2,4,6,8,10,12,14,16,20,18,17),c(1,2,3,1,1,2,3,2,3,1,3),method='parallel'))
    expect_equal(orderStrata(c(1,10,2:16),c(1,1,2:16),method='log'), orderStrata(c(1,10,2:16),c(1,1,2:16),method='parallel'))
    #
    expect_equal(orderStrata(c(2,6, 6,7),c(1,1, 2,2),method='log'), orderStrata(c(2,6, 6,7),c(1,1, 2,2),method='naive'))
    expect_equal(orderStrata(c(2,6, 6,7),c(2,2, 1,1),method='log'), orderStrata(c(2,6, 6,7),c(2,2, 1,1),method='naive'))
    expect_equal(orderStrata(c(2,5, 4, 6,7, 8),c(3,3, 1, 2,2, 4),method='log'), orderStrata(c(2,5, 4, 6,7, 8),c(3,3, 1, 2,2, 4),method='naive'))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(1,1,1,1,1, 2,2,2, 3,3,3),method='log'), orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(1,1,1,1,1, 2,2,2, 3,3,3),method='naive'))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 2,2,2, 1,1,1),method='log'), orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 2,2,2, 1,1,1),method='naive'))
    expect_equal(orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 1,1,1, 2,2,2),method='log'), orderStrata(c(2,8,17,10,18, 4,12,16, 6,14,20),c(3,3,3,3,3, 1,1,1, 2,2,2),method='naive'))
    expect_equal(orderStrata(c(2,4,6,8,10,12,14,16,20,18,17),c(1,2,3,1,1,2,3,2,3,1,3),method='log'), orderStrata(c(2,4,6,8,10,12,14,16,20,18,17),c(1,2,3,1,1,2,3,2,3,1,3),method='naive'))
    expect_equal(orderStrata(c(1,10,2:16),c(1,1,2:16),method='log'), orderStrata(c(1,10,2:16),c(1,1,2:16),method='naive'))

})
