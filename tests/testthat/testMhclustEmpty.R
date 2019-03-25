context("mhclust: pathological cases")

library(mhca)

test_that("missing/empty/invalid input recognized", {
  expect_error(mhclust())
  expect_error(mhclust(c()))
  expect_error(mhclust(NULL))
  expect_error(mhclust(1))
  expect_error(mhclust(1:10))
})
