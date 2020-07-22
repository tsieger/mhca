catnl<-structure(
function # Cat with trailing newline.
##description<<
## \code{\link{cat}} with trailing newline.
(..., ##<< arguments to cat
sep=' '
) {
  cat(...,'\n',sep=sep)
  ### None.
},ex=function() {
  catnl('Hello, world!')
})
