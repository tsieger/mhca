printWithName<-function (x,dec.digits=4) {
    if (typeof(x)=='integer') {
        print.single <- function(x) cat(sprintf('%*d',dec.digits,x))
        print.matrix.elem <- function(x) cat(sprintf('% *d',width,x))
    } else {
        print.single <- function(x) cat(sprintf('%.*f',dec.digits,x))
        # todo: hardcoded dec.digits, use eval(parse(dec.digits))
        print.matrix.elem <- function(x) cat(sprintf('% *.4f',width,x))
    }
    cat(paste(deparse(substitute(x)),'\n',sep=''))
    if (is.numeric(x)) {
        if (is.matrix(x)) {
            width<-1
            if (isTRUE(any(x>0))) width<-max(width,ceiling(log10(max(x[x>0],na.rm=TRUE))))
            if (isTRUE(any(x<0))) width<-max(width,1+ceiling(log10(max(-x[x<0],na.rm=TRUE))))
            for (r in 1:nrow(x)) {
                cat('    ')
                for (c in 1:ncol(x)) {
                    print.matrix.elem(x[r,c])
                    cat('    ')
                }
                cat('\n')
            }
        }
        else {
          cat('    ')
          print.single(x)
          cat('\n')
        }
    } else print(x)
}
