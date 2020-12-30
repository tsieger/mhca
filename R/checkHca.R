checkHca<-structure(
function # Check HCA consistency.
##description<<
## This function checks that the results of HCA is consistent. It
## prints info on key components.
(
h, ##<< hierarchical clustering (e.g. the result of
## \code{\link[stats]{hclust}} or \code{\link[mhca]{mhclust}}).
verb = TRUE, ##<< verbosity level
dbg = 0 ##<< debug level
) {
  ok<-TRUE
  if (verb) {
    disp<-function(...) catnl(...)
  } else {
    disp<-function(...) {}
  }

  disp('\'merge\' component',ifelse(!is.null(h$merge),'present','*** MISSING ***'))
  if (is.null(h$merge)) ok<-FALSE
  disp('\'merge\' component of ',ifelse(is.matrix(h$merge),'correct','*** INCORRECT ***'),' type',sep='')
  if (!is.matrix(h$merge)) ok<-FALSE
  disp('\'height\' component',ifelse(!is.null(h$height),'present','*** MISSING ***'))
  if (is.null(h$height)) ok<-FALSE
  disp('\'order\' component',ifelse(!is.null(h$order),'present',' ***MISSING ***'))
  if (is.null(h$order)) ok<-FALSE
  #disp('\'labels\' component',ifelse(!is.null(h$labels),'present','*** MISSING ***'))
  #if (is.null(h$labels)) ok<-FALSE
  disp('\'method\' component',ifelse(!is.null(h$method),'present','*** MISSING '))
  if (is.null(h$method)) ok<-FALSE
  disp('\'dist.method\' component',ifelse(!is.null(h$dist.method),'present','*** MISSING ***'))
  if (is.null(h$dist.method)) ok<-FALSE

  # number of clusters
  n<-length(h$height)

  if (is.null(h$merge) || !is.matrix(h$merge)) h$merge<-matrix(0,1,2)
  disp('the number of rows of \'merge\' of ',nrow(h$merge),' ',ifelse(nrow(h$merge)==n,'corresponds','*** DOES NOT correspond ***'),' with the \'height\' of ',n,' entry(ies)',sep='')
  if (nrow(h$merge)!=n) ok<-FALSE
  disp('the length of \'order\' of ',length(h$order),' ',ifelse(length(h$order)==n+1,'corresponds','*** DOES NOT correspond ***'),' with the \'height\' of ',n,' entry(ies)',sep='')
  if (length(h$order)!=n+1) ok<-FALSE
  #if (length(h$labels)!=n+1) stop('the length of \'labels\'',length(h$labels),', but \'height\' has ',n,' entry(ies)',sep='')

  # check 'merge'
  obs<--h$merge[h$merge<0]
  if (dbg) {cat('sort(obs):\n');print(sort(obs))}
  disp(length(obs),' observation(s) ',ifelse(length(obs)==n+1,'correspond(s)','*** DO(ES) NOT correspond ***'),' with the \'height\' of ',n,' entry(ies)',sep='')
  if (length(obs)!=n+1) ok<-FALSE
  disp(length(obs),' unique observation ID(s) ',ifelse(length(unique(obs))==n+1,'correspond(s)','*** DO(ES) NOT correspond ***'),' with the \'height\' of ',n,' entry(ies)',sep='')
  if (length(unique(obs))!=n+1) ok<-FALSE
  if (length(obs)>0) {
      disp('min observation ID ',min(obs),' is ',ifelse(min(obs)==1,'','*** NOT *** '),'1',sep='')
      if (min(obs)!=1) ok<-FALSE
      disp('max observation ID ',max(obs),' is ',ifelse(min(obs)==1,'','*** NOT *** '),n+1,sep='')
      if (max(obs)!=n+1) ok<-FALSE
  }
  cls<-h$merge[h$merge>0]
  if (dbg) {cat('sort(cls):\n');print(sort(cls))}
  disp(length(cls),' cluster(s) ',ifelse(length(cls)==n-1,'correspond(s)','*** DO(ES) NOT correspond ***'),' with the \'height\' of ',n,' entry(ies)',sep='')
  if (length(cls)!=n-1) ok<-FALSE
  disp(length(cls),' unique cluster ID(s) ',ifelse(length(unique(cls))==n-1,'correspond(s)','*** DO(ES) NOT correspond ***'),' with the \'height\' of ',n,' entry(ies)',sep='')
  if (length(unique(cls))!=n-1) ok<-FALSE
  if (length(cls)>0) {
      disp('min cluster ID ',min(cls),' is ',ifelse(min(cls)==1,'','*** NOT *** '),'1',sep='')
      if (min(cls)!=1) ok<-FALSE
      disp('max cluster ID ',max(cls),' is ',ifelse(max(cls)==n-1,'','*** NOT *** '),n-1,sep='')
      if (max(cls)!=n-1) ok<-FALSE
  }
  mx<-apply(h$merge,1,max)
  i<-which(mx>seq(along=mx))
  if (length(i)>0) {
    disp('*** invalid clusters in \'merge\' in row ',i,': ',h$merge[i,1],', ',h$merge[i,2],' ***',sep='')
  } else {
    disp('clusters in \'merge\' OK')
  }

  # check 'order'
  disp(length(unique(h$order)),' unique numbers in \'order\' ',ifelse(length(unique(h$order))==n+1,'is','*** IS NOT ***'),' ',n+1,sep='')
  if (length(unique(h$order))!=n+1) ok<-FALSE
  if (length(h$order)>0) {
      disp('min order ',min(h$order),' is ',ifelse(min(h$order)==1,'','*** NOT *** '),' 1',sep='')
      if (min(h$order)!=1) ok<-FALSE
      disp('max order ',max(h$order),' is ',ifelse(max(h$order)==n+1,'','*** NOT *** '),' ',n+1,sep='')
      if (max(h$order)!=n+1) ok<-FALSE
  }

  return(ok)
  ### A boolean value indicating whether \code{h} is valid or not.

},ex=function() {
  hGood<-hclust(dist(1:3))
  checkHca(hGood)

  hBad<-list(merge=rbind(c(-1,-2),c(-1,2)),height=1:2,order=1:3,method='fake-method',dist.method='fake-dist.method')
  checkHca(hBad)
})
