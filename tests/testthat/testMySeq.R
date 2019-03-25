context("mySeq")

library(mhca)

test_that("mySeq", {
    expect_error(mySeq(NULL))
    expect_error(mySeq(c()))
    expect_error(mySeq(c(),c()))
    expect_error(mySeq(0,1,0))

    expect_equal(mySeq(0,1),c(0,1))
    expect_equal(mySeq(0,1,2),c(0))
    expect_equal(mySeq(0,2,2),c(0,2))
    expect_equal(mySeq(0,0,0),c(0))

    expect_length(mySeq(1,0),0)
    expect_equal(mySeq(1,0,-1),c(1,0))
})
