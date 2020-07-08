context("fixNonMonotHca")

library(mhca)

test_that("empty argument handled", {
  expect_error(fixNonMonot())
  expect_error(fixNonMonot(c()))
})

test_that("eps method works", {
  hc<-hclust(dist(1:2))
  rv<-fixNonMonotHca(hc,method='eps')
  expect_equal(rv$height,1)

  # we are lazy from now on: we simply manipulate the 'height', not the data that get rise to it

  hc$height<-c(1,2,3,2,2,6)
  rv<-fixNonMonotHca(hc,method='eps')
  expect_equal(rv$height,c(1,2,3,3.5,4,6))

  hc$height<-c(1,2,3,2,2,6)
  rv<-fixNonMonotHca(hc,method='eps',eps=.1)
  expect_equal(rv$height,c(1,2,3,3.1,3.2,6))

  hc$height<-c(1,0)
  rv<-fixNonMonotHca(hc,method='eps')
  expect_equal(rv$height,c(1,2))

  hc$height<-c(1,0)
  rv<-fixNonMonotHca(hc,method='eps',eps=.1)
  expect_equal(rv$height,c(1,1.1))

  hc$height<-c(1,2,3,2,2,6,7,2,2,10)
  rv<-fixNonMonotHca(hc,method='eps')
  expect_equal(rv$height,c(1,2,3,3.5,4,6,7,7.5,8,10))
})

test_that("halfway method works", {
  hc<-hclust(dist(1:2))
  rv<-fixNonMonotHca(hc,method='halfway')
  expect_equal(rv$height,1)

  # we are lazy from now on: we simply manipulate the 'height', not the data that get rise to it

  hc$height<-c(1,2,3,2,2,6)
  rv<-fixNonMonotHca(hc,method='halfway')
  expect_equal(rv$height,c(1,2,3,4,5,6))

  hc$height<-c(1,2,3,2,2,6,7,2,2,10)
  rv<-fixNonMonotHca(hc,method='halfway')
  expect_equal(rv$height,c(1,2,3,4,5,6,7,8,9,10))

})
