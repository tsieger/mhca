context("testGetGroupIdx")

library(mhca)

test_that("pathiological cases", {
    expect_error(getDistGroupIdx(NULL))
    expect_error(getDistGroupIdx(1,1:3))
})

test_that("basic", {
    expect_equal(getDistGroupIdx(6,1:3),c(1,2,6))
    expect_equal(getDistGroupIdx(6,5:6),c(15))
})
