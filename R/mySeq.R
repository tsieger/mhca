mySeq<-structure(function # mySeq
### Sequence generation resembling matlab ':' operator.
##details<<
## The difference from R `seq' is the behaviour in case when
## from=2, to=1 and by=1 - in that case R seq raises an error, while
## matlab returns an empty sequence.
## This function returns an empty sequence as well in this case.
(from,to,by=1) {

    if (!length(from)) stop('"from" must of length 1')
    if (!length(to)) stop('"to" must of length 1')
    if (!length(by)) stop('"by" must of length 1')
    if (by==0 && from!=to) stop('invalid "by"')

    if ((from<to) == (by>0)) return(seq(from,to,by))
    else if (from==to) return(from)
    else return(vector(class(from),0))
},ex=function() {
    mySeq(1,0)
})
