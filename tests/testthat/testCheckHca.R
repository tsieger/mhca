context("mhclust: checkHca")

library(mhca)

test_that("missing entries identified", {
    expect_equal(checkHca(c(),verb=FALSE),FALSE)
    expect_equal(checkHca(list(),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=1),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=1,height=1),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=1,height=1,order=1),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=NULL           ,height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=NULL,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=1,order=NULL,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=1,order=1:2,           dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=1,order=1:2,method='m'                ),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),TRUE)
})

test_that("invalid \'merge\' identified", {
    expect_equal(checkHca(list(merge=1:2,height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=1:2,height=1:2,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-0,-1)),height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-1)),height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-2,-2)),height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-3)),height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),FALSE)

    expect_equal(checkHca(list(merge=rbind(c(-1,-2),c(1,1)),height=1:2,order=1:3,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2),c(2,2)),height=1:2,order=1:3,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2),c(2,3)),height=1:2,order=1:3,method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2),c(2,1)),height=1:2,order=1:3,method='m',dist.method='m'),verb=FALSE),FALSE)
})

test_that("correct \'order\' passed", {
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=1,order=1:2,method='m',dist.method='m'),verb=FALSE),TRUE)
})

test_that("invalid \'order\' identified", {
    expect_equal(checkHca(list(merge=rbind(c(-1,-2)),height=1,order=0,method='m',dist.method='m'),verb=FALSE),FALSE)

    expect_equal(checkHca(list(merge=rbind(c(-1,-2),c(1,-3)),height=1:2,order=c(2,2,2),method='m',dist.method='m'),verb=FALSE),FALSE)
    expect_equal(checkHca(list(merge=rbind(c(-1,-2),c(1,-3)),height=1:2,order=c(1,1,3),method='m',dist.method='m'),verb=FALSE),FALSE)
})
