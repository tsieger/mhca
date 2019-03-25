/*
 * Mahalanobis distance-based hierarchical clustering.
 *
 * This is an optimization rewrite of the R implementation contained
 * in "../R/mhca_Rimpl.R". Nontrivial parts of the original R code are
 * preserved and prefixed with `//R:'.
 */

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

// enable debugs and progress info
#define ENABLE_DEBUGS

#ifdef ENABLE_DEBUGS
// debug print (unconditional)
# define DBGU(...) Rprintf(__VA_ARGS__)
// debug code run conditionally on the value of dbgLevel
# define DBG_CODE(dbgLevel,...) if (dbg>=dbgLevel) __VA_ARGS__
// debug print (conditional on the value of dbgLevel)
# define DBG(dbgLevel,...) DBG_CODE(dbgLevel,DBGU(__VA_ARGS__))
#else // debugs completely disabled
# define DBGU(...)
# define DBG_CODE(dbgLevel,...)
# define DBG(dbgLevel,...)
#endif

// assert
#if 01 // being on the safe side
# define ASSERT(cnd,msg) if (!(cnd)) error(msg)
#else // ignoring errors
# define ASSERT(cnd,msg)
#endif

// check for NULL arguments
#define IS_R_NULL(x) (x==R_NilValue)

// convert C 0-based indices to R 1-based indices
#define C2R(n) (n+1)
// convert R 1-based indices to C 1-based indices
#define R2C(n) (n-1)

// Memory allocation.
// By default, the R managed heap is used, which does not need to be
// freed explicitly (see R_alloc documentation, e.g. in "Writing R
// exntensions"). However, if another memory allocation gets used,
// calls to MEM_ALLOC should carefully be tied with calls to MEM_FREE,
// which is currently not the case with print* debug routines.
#define MEM_ALLOC(n,size) R_alloc(n,size)
#define MEM_FREE(p)

typedef int Num; // e.g. cluster number
typedef int LongNum; // e.g. index into distance matrix (needs to be an int in order e.g. to be passed to 'idamax')
typedef unsigned int MatrixIndex; // e.g. index into distance matrix (in our internal code)

/* Global shared string buffers to hold temporary strings.
 */
char *strBuf=NULL; // a long buffer (to hold a sequence of `O(n^2)' real numbers)
char *strBufShort=NULL; // a short buffer (to hold a sequence of `O(n)' integers)
char *strBufShort2=NULL; // a short buffer (to hold a sequence of `O(n)' integers)

#ifdef ENABLE_DEBUGS
/*
 * Get members of a cluster, a sequence of string form of IDs of observations.
 *
 * Arguments:
 *  clusterMembers - array of cluster members
 *  memberCount - number of members
 *  buf - buffer to store the string prepresentation of members
 */
const char *getMembers(const Num *clusterMembers,Num memberCount,char *buf) {
    Num i;
    Num pos=0;
    for (i=0;i<memberCount;i++) {
        sprintf(buf+pos,"%d",C2R(clusterMembers[i]));
        pos+=strlen(buf+pos);
        if (i<memberCount-1) {
            buf[pos++]=' ';
        }
    }
    buf[pos]=0;
    return buf;
}

/*
 * Print matrix of doubles.
 *
 * Arguments:
 *  name - matrix identification
 *  x - matrix
 *  rows - number of rows
 *  cols - number of columns
 */
void printDoubleMatrix(const char *name,const double *x,MatrixIndex rows,MatrixIndex cols) {
    //Rprintf("printDoubleMatrix(cols=%d)\n",cols);
    MatrixIndex i,j;
    Rprintf("%s:\n",name);
    if (rows>0 && cols>0) {
        for (i=0;i<rows;i++) {
            MatrixIndex pos=0;
            for (j=0;j<cols;j++) {
                sprintf(strBuf+pos,"% .4f%s",x[i+rows*j],j<cols-1?"  ":"");
                pos+=strlen(strBuf+pos);
            }
            strBuf[pos]=0;
            Rprintf("    %s\n",strBuf);
        }
    }
}

/*
 * Print a row of a matrix of doubles.
 *
 * Arguments:
 *  name - matrix identification
 *  x - matrix
 *  rows - number of rows
 *  cols - number of columns
 *  row - row index
 */
void printDoubleMatrixRow(const char *name,const double *x,MatrixIndex rows,MatrixIndex cols,MatrixIndex row) {
    MatrixIndex j;
    Rprintf("%s:\n",name);
    if (rows>0 && cols>0) {
        MatrixIndex pos=0;
        for (j=0;j<cols;j++) {
            sprintf(strBuf+pos,"% .4f%s",x[row+rows*j],j<cols-1?"  ":"");
            pos+=strlen(strBuf+pos);
        }
        strBuf[pos]=0;
        Rprintf("    %s\n",strBuf);
    }
}

/*
 * Print distance matrix (the lower triangle annotated with indices).
 *
 * Arguments:
 *   name - matrix identification
 *   x - distance matrix
 *   n - number of observations the distance matrix was computed for
 *   inversionFactor - the distances are `inversionFactor - x[...]'
 */
void printDistMatrix(const char *name,const double *x,Num n,double inversionFactor) {
    Num i,j,pos;

    /*    1 2 3 4
     * 2  1
     * 3  2 5
     * 4  3 6 8
     * 5  4 7 9 10
     *
     * (i,j) @ n*(i-1) - i*(i-1)/2 + j-i
     * i = cluster from, in columns
     * j = cluster to, in rows
     */
    Rprintf("%s:\n",name);
    if (n>0) {

        // print header
        pos=0;
        for (i=0;i<n-1;i++) { // i iterates over columns
            sprintf(strBuf+pos," %9d",C2R(i));
            pos+=strlen(strBuf+pos);
        }
        strBuf[pos]=0;
        Rprintf("        %s\n",strBuf);

        // print data rows
        for (j=1;j<n;j++) { // cluster to, in rows
            pos=0;
            for (i=0;i<n-1;i++) { // cluster from, in columns
                if (i<j) {
                    // note the '(long)x' below not to let the 'int*int' overflow
                    sprintf(strBuf+pos," %9.6f",inversionFactor-x[R2C(((long)n)*(C2R(i)-1) - ((long)C2R(i))*(C2R(i)-1)/2 + C2R(j)-C2R(i))]);
                } else {
                    sprintf(strBuf+pos,"          ");
                }
                pos+=strlen(strBuf+pos);
            }
            strBuf[pos]=0;
            // prepend and append row number
            Rprintf("    %3d %s    %3d \n",j+1,strBuf,j+1);
        }

        // print trailer
        pos=0;
        for (i=0;i<n-1;i++) {
            sprintf(strBuf+pos," %9d",i+1);
            pos+=strlen(strBuf+pos);
        }
        strBuf[pos]=0;
        Rprintf("        %s\n",strBuf);
    }
}

/*
 * Print matrix of integers (Num's).
 *
 * Arguments:
 *  name - matrix identification
 *  x - matrix
 *  rows - number of rows
 *  cols - number of columns
 */
void printNumMatrix(const char *name,const Num *x,MatrixIndex rows,MatrixIndex cols) {
    MatrixIndex i,j;
    Rprintf("%s:\n",name);
    if (cols>0) {
        for (i=0;i<rows;i++) {
            MatrixIndex pos=0;
            for (j=0;j<cols;j++) {
                sprintf(strBuf+pos,"%d%s",x[i+rows*j],j<cols-1?"  ":"");
                pos+=strlen(strBuf+pos);
            }
            strBuf[pos]=0;
            Rprintf("    %s\n",strBuf);
        }
    }
}

/*
 * Print matrix of long integers (LongNum's).
 *
 * Arguments:
 *  name - matrix identification
 *  x - matrix
 *  rows - number of rows
 *  cols - number of columns
 */
void printLongNumMatrix(const char *name,const LongNum *x,MatrixIndex rows,MatrixIndex cols) {
    MatrixIndex i,j;
    Rprintf("%s:\n",name);
    if (cols>0) {
        for (i=0;i<rows;i++) {
            MatrixIndex pos=0;
            for (j=0;j<cols;j++) {
                sprintf(strBuf+pos,"%d%s",x[i+rows*j],j<cols-1?"  ":"");
                pos+=strlen(strBuf+pos);
            }
            strBuf[pos]=0;
            Rprintf("    %s\n",strBuf);
        }
    }
}
#endif // ENABLE_DEBUGS

/*
 * Compute the sample covariance estimate.
 *
 * Arguments:
 *   x - data matrix of size `r x c', stored in column-major style
 *   r - number of rows of `x'
 *   c - number of columns of `x'
 *   res - where to hold the resulting `c x c' matrix
 *   covMeansTmp - temporary storage for means of x
 */
double *cov(double *x, Num r, Num c, double *res,double *covMeansTmp) {
    Num i,i1,i2,j,offset,offset1,offset2,offset3,offset4;
    double *means=covMeansTmp;
    for (offset=i=0;i<c;i++) {
        double sum=0.0;
        for (j=0;j<r;j++) {
            // sum+=x[j+i*r];
            sum+=x[j+offset];
        }
        means[i]=sum/r;
        offset+=r;
    }
    for (offset1=offset3=i1=0;i1<c;i1++) {
        for (offset2=offset4=i2=0;i2<c;i2++) {
            if (i2>=i1) {
                /* computing lower diagonal */
                double sum=0.0;
                for (j=0;j<r;j++) {
                    // sum+=(x[j+i1*r]-means[i1])*(x[j+i2*r]-means[i2]);
                    sum+=(x[j+offset1]-means[i1])*(x[j+offset2]-means[i2]);
                }
                // res[i1+i2*c]=sum/(r-1);
                res[i1+offset4]=sum/(r-1);
            } else {
                /* computing upper diagonal from lower diagonal */
                // res[i1+i2*c]=res[i2+i1*c];
                res[i1+offset4]=res[i2+offset3];
            }
            offset2+=r;
            offset4+=c;
        }
        offset1+=r;
        offset3+=c;
    }
    return res;
}

/*
 * Compute the mean Mahalanobis distance of vectors in rows of a
 * matrix. The Mahalanobis distance of a vector `Xi' of size `p x 1',
 * is defined as
 *
 *    d = sqrt( Xi' * IC * Xi )
 *
 * where IC is the inverse covariance matrix estimate.
 *
 * Arguments:
 *   X - matrix formed by vectors `Xi' in rows, stored in column-major style
 *   n - number of rows of `X'
 *   p - number of columns of `X'
 *   IC - inverse of the covariance matrix
 *   buf - buffer of working matrix of size `n x p'
 *   dbg
 */
double computeMahalanobisDistance(double *X,Num n,Num p,double *IC,double *buf,int dbg) {
    Num i,j;
    double r;

    /*
     * dsymm(const char *side, const char *uplo, const int *m,const int *n, const double *alpha,
     *       const double *a, const int *lda,const double *b, const int *ldb,const double *beta,
     *       double *c, const int *ldc);
     * in general:
     *    alpha*B * A + beta*C  ->  C, A symmetric
     * we use only:
     *    X    *    IC    ->   R
     * {n x p} * {p x p}    {n x p}
     *
     * see http://www.math.utah.edu/software/lapack/lapack-blas/dsymm.html
     */
    double alpha=1.0;
    double beta=0.0;
    F77_CALL(dsymm)("r",// A appears on the right side, i.e. B*A
                    "l",// lower part of A is to be referenced
                    &n,&p, // rowR, colR
                    &alpha,
                    IC,&p, // IC, ld = p
                    X,&n, // X, ld = n
                    &beta,
                    buf,&n); // R, ld = n
    DBG_CODE(4,printDoubleMatrix("result",buf,n,p));

    // mean Mahalanobis distance from origin to rows of X
    // R: mean(sqrt(rowSums((X%*%IC)*X)))
    for (i=0;i<n;i++) { // TODO: optimize?
        double sum=0.0;
        Num offset;
        // sum rows
        for (offset=j=0;j<p;j++) {
            // sum+=buf[i+n*j]*X[i+n*j];
            sum+=buf[i+offset]*X[i+offset];
            offset+=n;
        }
        // take the square root and store in the first column of "result"
        buf[i]=sqrt(sum);
    }
    i=1; // lda of result
    r=F77_NAME(dasum)(&n,buf,&i)/n;
    return r;
}

/*
 * Construct matrix `Xi', a matrix representing cluster `i'. The
 * resulting matrix either consists of all cluster members, or of the
 *  centroid of the cluster.
 *
 * Arguments:
 *  x - matrix of size `clusterCount x p', stored in column-major style
 *  n - number of rows of `x'
 *  p - number of columns of `x'
 *  centroid - matrix of centroids, of size `clusterCount x p',
 *      stored in column-major style
 *  cluster - index of cluster for which to compute the representants
 *  clusterMembers - members of the given cluster
 *  memberCount - number of members of the given cluster
 *  y - where to store the resulting matrix, of size `yn x p'
 *  ynPtr - pointer to hold `yn', the number of rows of the result
 *  quickMode - boolean flag determining whether to use all members of
 *      clusterFrom, or only its centroid
 * dbg
 */
Num constructRepresentantsOfCluster(double *x,Num n,Num p,double *centroid,
                                    Num cluster,Num *clusterMembers,Num memberCount,
                                    double *y,int quickMode,int dbg) {
    Num i,j,yn;

    if (quickMode) {
        Num offset;
        yn=1;
        for (offset=i=0;i<p;i++) {
            // y[i]=centroid[cluster+n*i];
            y[i]=centroid[cluster+offset];
            offset+=n;
        }
    } else {
        Num offset1,offset2;
        yn=memberCount;
        for (offset1=offset2=i=0;i<p;i++) {
            Num cumI=0;
            // offset1=yn*i;
            // offset2=n*i;
            for (j=0;j<yn;j++,cumI++) {
                y[cumI+offset1]=x[clusterMembers[j]+offset2];
            }
            offset1+=yn;
            offset2+=n;
        }
    }
    return(yn);
}

/*
 * Construct matrix `Xi - cj', where `Xi' is the matrix representing
 * cluster `i' (either all its members or its centroid), and `cj'
 * represents the centroid of the cluster `j', repeated in rows to form
 * a matrix conformable to `Xi'.
 *
 * Arguments:
 *  Xi - matrix representing cluster `i' (usually either its members,
 *      or its centroid), stored in column-major style
 *  n - number of rows of `Xi'
 *  p - number of columns of `Xi'
 *  centroid - matrix of centroids of size `clusterCount x p'
 *  clusterCount - number of clusters = number of rows of `centroid'
 *  clusterTo - index of target cluster
 *  y - where to store the resulting matrix of size `n x p'
 *  quickMode - boolean flag determining whether to use all members of
 *      clusterFrom, or only its centroid
 * dbg
 */
void constructBetweenClusterDistanceMatrixFromXi(double *Xi,Num n,Num p,
                                                 double *centroid,Num clusterCount,
                                                 Num clusterTo,double *y,
                                                 int quickMode,int dbg) {
    Num i,j,offset1,offset2;

    DBG_CODE(4,printDoubleMatrix("y",Xi,n,p));
    DBG_CODE(4,printDoubleMatrixRow("centroid[clusterTo]",centroid,clusterCount,p,clusterTo));
    //R: xc1<-xc1-matrix(centroid[c1,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE)
    for (offset1=offset2=i=0;i<p;i++) {
        // offset1=n*i;
        // offset2=clusterCount*i;
        for (j=0;j<n;j++) {
            y[j+offset1]=Xi[j+offset1]-centroid[clusterTo+offset2];
        }
        offset1+=n;
        offset2+=clusterCount;
    }
    DBG_CODE(4,printDoubleMatrix("y-centroid",y,n,p));
}

/*
 * Construct matrix `Xi - cj', where `Xi' is the matrix representing
 * cluster `i' (either all its members or its centroid), and `cj'
 * represents the centroid of the cluster `j', repeated in rows to form
 * a matrix conformable to `Xi'.
 *
 * Arguments:
 *  x - `clusterCount x p' data matrix, stored in column-major style
 *  p - number of columns of `x'
 *  n - number of rows of `x'
 *  clusterCount - number of clusters = number of rows of `centroid'
 *  centroid - matrix of centroids of size `clusterCount x p'
 *  clusterFrom - index of source cluster
 *  clusterTo - index of target cluster
 *  clusterMembers - members of the source cluster
 *  memberCount - number of members of the source cluster
 *  y - where to store resulting the matrix of size `yn x p'
 *  ynPtr - pointer to hold `yn', the number of rows of the result
 *  quickMode - boolean flag determining whether to use all members of
 *      clusterFrom, or only its centroid
 * dbg
 */
Num constructBetweenClusterDistanceMatrix(double *x,Num n,Num p,
                                          double *centroid,Num clusterCount,
                                          Num clusterFrom,Num clusterTo,
                                          Num *clusterMembers,Num memberCount,
                                          double *y,
                                          int quickMode,int dbg) {
    Num i,j,yn;
    yn=constructRepresentantsOfCluster(x,n,p,centroid,clusterFrom,
                                       clusterMembers,memberCount,
                                       y,quickMode,dbg);
    DBG_CODE(4,printDoubleMatrix("y",y,yn,p));
    DBG_CODE(4,printDoubleMatrixRow("centroid[clusterTo]",centroid,clusterCount,p,clusterTo));
    //R: xc1<-xc1-matrix(centroid[c1,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE)
    Num offset1,offset2;
    for (offset1=offset2=i=0;i<p;i++) {
        // offset1=yn*i;
        // offset2=clusterCount*i;
        for (j=0;j<yn;j++) {
            y[j+offset1]-=centroid[clusterTo+offset2];
        }
        offset1+=yn;
        offset2+=clusterCount;
    }
    DBG_CODE(4,printDoubleMatrix("y-centroid",y,yn,p));
    return yn;
}

/*
 * Mahalanobis distance-based hierarchical clustering.
 *
 * Arguments:
 *  X - data matrix of size `n x p', stored in column-major style
 *  DistX - distance matrix of X, stored as a vector holding the lower
 *       triangle of the distance matrix in the column-major style
 *  Merging - where to store the `2 x (n-1)' matrix describing the
 *       iterative merging of observations and clusters (column-major
 *       style)
 *  Height - where to store the heights of the `n'-1 clusters
 *  Thresh - the Mahalanobis threshold in the interval of (0,1) defining
 *       the minimal relative size of cluster whose distance to other
 *       clusters will be computed as a pure Mahalanobis distance.
 *  Quick - if nonzero, inter-cluster distances will be computed using
 *       centroids only, not using all the observations contained
 *       in the clusters.
 *  Normalize - if nonzero, cluster size will be ignored when computing
 *       Mahalanobis distance from the cluster. If zero, once all
 *       clusters are of at least the `Thresh' relative size, both
 *       the cluster shape and size will affect inter-cluster distance.
 *  Verb - level of verbosity, the greater the more detailed info,
 *      defaults to 0 (no info)
 * _NFull, _NLeft, _Centroid, _Members, _Invcov, _DetsSqrt,
 * _WeightFactor, _ClusterId, _Height are internal parameters used when
 *      clustering pre-clustered apriori clusters)
 */
SEXP mhclust_(SEXP X,SEXP DistX,SEXP Merging,SEXP Height,SEXP Thresh,SEXP Quick,SEXP Normalize,SEXP G,SEXP GMergingCount,SEXP Verb,
              SEXP _NFull,SEXP _NLeft,SEXP _Centroid,SEXP _Members,SEXP _Invcov,SEXP _DetsSqrt,SEXP _WeightFactor,
              SEXP _ClusterId,SEXP _ClusterSize,SEXP _MembersPoolSize) {
    int dbg;
    int quick;

    int nprot = 0;

    double *x; // input data
    Num *dimX;
    Num n,n1,n2; // number of data points
    Num nFull; // number of total data points (it equals n on primary
    // call, but it is greater than n on recursive calls)
    Num clusterCount; // number of clusters (=n-1)
    Num p; // dimensionality of x
    double *centroid; // centroids of clusters
    double *xij; // data of the merged cluster
    Num xijN;
    double *xc1,*xc2; // data of clusters

    double *covXij; // cov of xij
    Num covXijLength;
    double *invcov,*invcovMerged,*fakeInvCov; // inverse covariance matrices
    double *ic1,*ic2; // inverse covariance matrices
    double *covMeansTmp; // temporary storage for means of x during computation of cov

    double *distX,*distXTmp; // distance matrix of x
    double maxDistX; // max value of distX
    LongNum distXIdx,*distXIndices; // indices into distX
    LongNum distXLen;

    Num *merging; // output argument: cluster merging
    double *height; // output argument: cluster heights

    double thresh;
    int normalize;

    Num *g,gMergingCount;

    Num *clusterSize; // the size of clusters
    Num *clusterId; // id of cluster contained in this specific position
    double detSqrt;
    double* detsSqrt; // square root of determinants of inverse covariance matrices

    double *weightFactor,wf1; // factors weighting inverse covariance matrix against fakeInvCov

    Num **members;
    Num *membersPool;
    unsigned long membersPoolPos;
    unsigned long membersPoolSize;

    Num *otherClusters;
    Num *otherClustersTmp;
    Num c1,c2,s,i,j,k,oi;
    Num clusterSizeI,clusterSizeJ;
    Num c1n,c2n,*c1indices,*c2indices;
    Num offset,offset1,offset2;

    Num fullMahalClusterCount; // number of clusters for which full Mahalanobis distance can be computed
    int switchedToFullMahal;

    double distMaha1,distMaha2;
    double *y; // temporal storage, used for `y' for `constructBetweenClusterDistanceMatrix'
    Num yn; // number of data points in `y'
    double *distMahaBuf; // working matrix used to compute Mahalanobis distance

    int lda;
    int info;

    dbg=*INTEGER(Verb);
    DBG(2,"mhclust_ called\n");
    DBG(2,"verb: %d\n",dbg);
    /*
     * dbg = 0: no debugs
     * dbg = 1: brief progress on each step (one line)
     * dbg = 2: more detailed progress on each step (several lines) and
     *          info on input and output arguments
     * dbg = 3: even more detailed progress - many lines of values of internal variables
     * dbg = 4: even more detailed info
     * dbg = 5: debug infos
     */

    Num _nFull;
    if (!IS_R_NULL(_NFull)) {
        _nFull=*INTEGER(_NFull);
    } else {
        _nFull=0;
    }

    // sanity check: the optional '_*' arguments starting from '_nLeft' should all be present or all be NULL
    int tmp[]={
        IS_R_NULL(_NLeft),
        IS_R_NULL(_Centroid),
        IS_R_NULL(_Members),
        IS_R_NULL(_Invcov),
        IS_R_NULL(_DetsSqrt),
        IS_R_NULL(_WeightFactor),
        IS_R_NULL(_ClusterId),
        IS_R_NULL(_ClusterSize),
        IS_R_NULL(_MembersPoolSize),
        -1};
    int allOptPresent=1,allOptMissing=1;
    for (i=0;tmp[i]!=-1;i++) {
        if (tmp[i]) allOptPresent=0;
        if (!tmp[i]) allOptMissing=1;
    }
    ASSERT(allOptPresent || allOptMissing,"invalid internal parameters");
    Num _nLeft=0L;
    double *_centroid=NULL;
    double *_detsSqrt=NULL;
    double *_weightFactor=NULL;
    Num *_clusterId=NULL;
    Num *_clusterSize=NULL;
    long _membersPoolSize=0L;
    if (allOptPresent) {
        _nLeft=*INTEGER(_NLeft);
        _centroid=REAL(_Centroid);
        _detsSqrt=REAL(_DetsSqrt);
        _weightFactor=REAL(_WeightFactor);
        _clusterId=INTEGER(_ClusterId);
        _clusterSize=INTEGER(_ClusterSize);
        _membersPoolSize=*INTEGER(_MembersPoolSize);
    }

    dimX=INTEGER(GET_DIM(X));
    PROTECT(X=AS_NUMERIC(X)); nprot++;
    PROTECT(DistX=AS_NUMERIC(DistX)); nprot++;
    x=REAL(X);
    distX=REAL(DistX);
    distXLen=LENGTH(DistX);
    distXTmp=(double *)MEM_ALLOC(distXLen,sizeof(*distXTmp));
    merging=INTEGER(Merging);
    height=REAL(Height);
    thresh=*REAL(Thresh);
    DBG(2,"thresh: %g\n",thresh);
    normalize=*LOGICAL(Normalize);
    DBG(2,"normalize: %d\n",normalize);
    quick=*LOGICAL(Quick);
    DBG(2,"quick: %d\n",quick);
    if (IS_R_NULL(G)) {
        g=NULL;
    } else {
        g=INTEGER(G);
    }
    DBG(2,"g: %p\n",g);
    gMergingCount=*INTEGER(GMergingCount);
    DBG(2,"gMergingCount: %d\n",gMergingCount);

    n=dimX[0]; // number of elementary points subject to clustering (merging)
    DBG(2,"n: %d\n",n);

    // how may clusters (hot merged together yet) remain; in the
    // beginning, all points are clusters, in the end only one huge
    // cluster exists
    if (!IS_R_NULL(_NLeft)) {
        clusterCount=_nLeft;
    } else {
        clusterCount=n;
    }
    DBG(2,"clusterCount: %d\n",clusterCount);

    // check that the distance matrix can be indexed using a 'LongNum' index
    DBG(2,"n*(n-1)/2: %ld\n",((long)clusterCount)*(clusterCount-1)/2);
    DBG(2,"max representable length: %ld\n",(1l<<(sizeof(LongNum)*8-1))-1);
    if (((long)clusterCount)*(clusterCount-1)/2>(1l<<(sizeof(LongNum)*8-1))-1) {
        error("number of clusters too large");
    }
#ifdef ENABLE_DEBUGS
    // Allocate the global shared string buffer that should hold the string
    // representation of a list of cluster members (i.e. at most `n' integers)
    // as well as the representation of the distance matrix (i.e. `n*(n-1)/2'
    // real values). Thus, we reserve `c*n' space with quite generous constant
    // `c' to hold e.g. long decimal expansion of real-valued distances.
    DBG(2,"allocating %lu\n",20l*n*(n-1)/2*sizeof(char)); // note the '20l' not to let the 'int*int' overflow
    strBuf=MEM_ALLOC(20l*n*(n-1)/2,sizeof(char)); // note the '20l' not to let the 'int*int' overflow
    // two more another string buffers (shorter ones)
    strBufShort=MEM_ALLOC(20*n,sizeof(char));
    strBufShort2=MEM_ALLOC(20*n,sizeof(char));
#endif
    DBG(2,"length of distX: %ld\n",(long)distXLen);
    DBG(3,"expected distX len: %ld\n",((long)clusterCount)*(clusterCount-1)/2); // note the '(long)x' not to let the 'int*int' overflow
    if (distXLen!=((long)clusterCount)*(clusterCount-1)/2) error("invalid distXLen");
    p=dimX[1]; // number of dimensions of the feature space
    DBG(2,"p: %d\n",p);

    if (!IS_R_NULL(_NFull)) {
        nFull=_nFull;
    } else {
        nFull=n;
    }
    DBG(2,"nFull: %d\n",nFull);

    covMeansTmp=(double *)MEM_ALLOC(p,sizeof(double));
    DBG_CODE(3,printDoubleMatrix("x",x,n,p));
    covXijLength=p*p;

    DBG_CODE(3,printDoubleMatrix("distX",distX,1,distXLen));
    // convert distX to C-distX, such that we could find the maximum absolute value
    // instead of the minimum value
    lda=1;
    k=F77_CALL(idamax)(&distXLen,distX,&lda);
    ASSERT(k>0,"invalid idamax retcode");
    maxDistX=distX[k-1];
    DBG(4,"found maximum %g at %d\n",maxDistX,k-1);
    DBG_CODE(4,printDoubleMatrix("distX (orig)",distX,1,distXLen));
    for (i=0;i<distXLen;i++) distX[i]=maxDistX-distX[i];
    DBG_CODE(4,printDoubleMatrix("distX (tx'd)",distX,1,distXLen));

    // number of elementary points in each cluster
    clusterSize=(int *)MEM_ALLOC(clusterCount,sizeof(int));
    if (!IS_R_NULL(_ClusterSize)) {
        for (i=0;i<clusterCount;i++) clusterSize[i]=_clusterSize[i];
    } else {
        for (i=0;i<clusterCount;i++) clusterSize[i]=1;
    }
    members=(Num**)MEM_ALLOC(clusterCount,sizeof(Num*));
    if (!IS_R_NULL(_Members)) {
        membersPoolSize=_membersPoolSize;
        membersPool=(Num*)MEM_ALLOC(membersPoolSize,sizeof(Num));
        DBG(3,"membersPool %p: allocated %ld entries (%ld B)\n",membersPool,membersPoolSize,membersPoolSize*sizeof(Num));
        membersPoolPos=0;
        for (i=0;i<clusterCount;i++) {
            members[i]=membersPool+membersPoolPos;
            membersPoolPos+=clusterSize[i];
            ASSERT(membersPoolPos<=membersPoolSize,"membersPool exhausted");
            Num *tmp=INTEGER(VECTOR_ELT(_Members,i));
            for (j=0;j<clusterSize[i];j++) {
                members[i][j]=R2C(tmp[j]);
            }
        }
    } else {
        // To store the cluster members during the clustering, there is at
        // most the cummulative need for
        // `1 + 1 + ... + 1 + 2 + 3 + ... + n = n+(n-1)*(n+2)/2' members.
        // (The first `n' `1s' stand for individual observations, then the
        // worst scenario is a single growing cluster accumulating a single
        // observation on each step.) This is a worst-case quadratic space,
        // while the optimal case would take only `O(n*log2(n))' space.
        membersPoolSize=clusterCount+((long)clusterCount-1)*(clusterCount+2)/2;
        membersPool=(Num*)MEM_ALLOC(membersPoolSize,sizeof(Num));
        DBG(3,"membersPool %p: allocated %ld entries (%ld B)\n",membersPool,membersPoolSize,membersPoolSize*sizeof(Num));
        membersPoolPos=0;
        // create new lists of members
        for (i=0;i<clusterCount;i++) {
            members[i]=membersPool+membersPoolPos++;
            ASSERT(membersPoolPos<=membersPoolSize,"membersPool exhausted");
            *members[i]=i;
        }
    }

    // clusters being made (by merging two smaller clusters) are
    // assigned unique IDs, but reside in data structured indexed by
    // index of one of its subclusters - thus we need to map the
    // 1:clusterCount space into IDs of current clusters
    clusterId=(int *)MEM_ALLOC(clusterCount,sizeof(int));
    if (!IS_R_NULL(_ClusterId)) {
        for (i=0;i<clusterCount;i++) clusterId[i]=R2C(_clusterId[i]);
    } else {
        for (i=0;i<clusterCount;i++) clusterId[i]=i;
    }

    xij=(double *)MEM_ALLOC(n*p,sizeof(double));
    covXij=(double *)MEM_ALLOC(covXijLength,sizeof(double));
    fakeInvCov=(double *)MEM_ALLOC(p*p,sizeof(double));
    memset(fakeInvCov,0,p*p*sizeof(double));
    for (i=0;i<p;i++) {
        fakeInvCov[(p+1)*i]=1;
    }
    // if inverse of covariance matrix can't be computed, use
    // this surrogate
    //R: fakeInvCov<-diag(p)
    DBG_CODE(3,printDoubleMatrix("fakeInvCov",fakeInvCov,p,p));
    ic1=(double *)MEM_ALLOC(p*p,sizeof(double));
    ic2=(double *)MEM_ALLOC(p*p,sizeof(double));

    otherClusters=(Num *)MEM_ALLOC(clusterCount*p,sizeof(Num));
    otherClustersTmp=(Num *)MEM_ALLOC(clusterCount*p,sizeof(Num));

    distXIndices=(LongNum*)MEM_ALLOC(clusterCount-1,sizeof(LongNum));

    y=(double *)MEM_ALLOC(n*p,sizeof(double));
    distMahaBuf=(double*)MEM_ALLOC(n*p,sizeof(double));

    // inverse of covariance matrix of members of given cluster
    // (representing the shape of the clusters)
    //R: invcov<-rep(list(fakeInvCov),clusterCount)
    invcov=(double *)MEM_ALLOC(clusterCount*p*p,sizeof(double));
    invcovMerged=(double *)MEM_ALLOC(clusterCount*p*p,sizeof(double));
    if (!IS_R_NULL(_Invcov)) {
        for (i=0;i<clusterCount;i++) {
            memcpy(invcov+p*p*i,REAL(VECTOR_ELT(_Invcov,i)),p*p*sizeof(double));
        }
    } else {
        memset(invcov,0,clusterCount*p*p*sizeof(double));
        for (i=0;i<clusterCount;i++) {
            for (j=0;j<p;j++) {
                invcov[p*p*i+(p+1)*j]=1;
            }
        }
    }
    // invcov holds inverses of covariance matrices, stored one after another
    // e.g. like ic1_11,ic1_21,ic1_12,ic1_22, ic2_11,ic2_21,ic2_12,ic2_22.
    // Normalizing factor for each cluster (computed from determinant of
    // the inverse of the covariance matrix) making the N-dim volume of
    // clusters equal to 1 if `invcov[[i]]' gets divided by `detsSqrt[i]'.
    detsSqrt=(double*)MEM_ALLOC(clusterCount,sizeof(double));
    if (!IS_R_NULL(_DetsSqrt)) {
        for (i=0;i<clusterCount;i++) detsSqrt[i]=_detsSqrt[i];
    } else {
        for (i=0;i<clusterCount;i++) detsSqrt[i]=1.0;
    }

    // centroids of clusters
    //R: centroid<-x
    centroid=(double *)MEM_ALLOC(n*p,sizeof(double));
    if (!IS_R_NULL(_Centroid)) {
        memcpy(centroid,_centroid,sizeof(*_centroid)*clusterCount*p);
    } else {
        memcpy(centroid,x,sizeof(*x)*clusterCount*p);
    }
    DBG_CODE(3,printDoubleMatrix("centroid",centroid,clusterCount,p));

    xc1=(double *)MEM_ALLOC(n*p,sizeof(double));
    xc2=(double *)MEM_ALLOC(n*p,sizeof(double));

    // proportional size of each cluster
    weightFactor=(double *)MEM_ALLOC(clusterCount,sizeof(double));
    if (!IS_R_NULL(_WeightFactor)) {
        for (i=0;i<clusterCount;i++) weightFactor[i]=_weightFactor[i];
    } else {
        for (i=0;i<clusterCount;i++) weightFactor[i]=0.0;
    }
    // number of clusters whose relative size is at least
    // mahalanobis.distance.threshold

    //R: fullMahalClusterCount<-sum(clusterSize>=thresh*fullPointCount)
    fullMahalClusterCount=0;
    for (i=0;i<clusterCount;i++) {
        if (clusterSize[i]>=thresh*nFull) fullMahalClusterCount++;
    }
    DBG(3,"fullMahalClusterCount=%d\n",fullMahalClusterCount);
    // have all clusters reached the thresh relative size and have we,
    // therefore, switched into "full Mahalanobis" mode?
    switchedToFullMahal=fullMahalClusterCount==clusterCount;

    ///////////////////////////////////////////////////////////////////
    // the main loop
    ///////////////////////////////////////////////////////////////////
    // merge two closest clusters at each step `s'
    //R: for (s in mySeq(pointCount-clusterCount+1,pointCount-1)) {
    for (s=n-clusterCount;s<n-1;s++) {
        double v;
        LongNum source,target,count;
        Num distXIndicesLen;

        DBG_CODE(3, {
            DBGU("\n====================== step %d ============================\n",s);
            for (i=0;i<clusterCount;i++) {
                DBGU("members[%d]: %s\n",C2R(i),getMembers(members[i],clusterSize[i],strBuf));
            }
            for (i=0;i<clusterCount;i++) {
                DBGU("invcov[%d]: \n",i);
                printDoubleMatrix("invcov",invcov+i*p*p,p,p);
            }
            //printDoubleMatrix("detsSqrt",detsSqrt,1,clusterCount);
            for (i=0;i<clusterCount;i++) {
                DBGU("detsSqrt[%d]: %g\n",C2R(i),detsSqrt[i]);
            }
            printDoubleMatrix("centroid",centroid,clusterCount,p);
            for (i=0;i<clusterCount;i++) {
                DBGU("weightFactor[%d]: %g\n",C2R(i),weightFactor[i]);
            }
            //printDoubleMatrix("weightFactor",weightFactor,1,clusterCount);
            DBGU("clusterCount=%d\n",clusterCount);
            //printNumMatrix("clusterId",clusterId,1,clusterCount);
            for (i=0;i<clusterCount;i++) {
                DBGU("clusterId[%d]: %d\n",C2R(i),C2R(clusterId[i]));
            }
            //printNumMatrix("clusterSize",clusterSize,1,clusterCount);
            for (i=0;i<clusterCount;i++) {
                DBGU("clusterSize[%d]: %d\n",C2R(i),clusterSize[i]);
            }
            DBGU("fullMahalClusterCount=%d\n",fullMahalClusterCount);
        });

        DBG_CODE(4,{
            DBGU("maxDistX=%f\n",maxDistX);
            DBGU("distXLen=%d\n",distXLen);
            for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
            printDoubleMatrix("distX",distXTmp,1,distXLen);
        });
        DBG_CODE(3,printDistMatrix("distX",distX,clusterCount,maxDistX));

        // find two mutually nearest clusters
        //R: k<-which_min(distX,distXLen);
        // IDAMAX - return the index of the element with max abs value
        lda=1;
        k=F77_CALL(idamax)(&distXLen,distX,&lda);
        ASSERT(k>0,"invalid idamax retcode");
        v=maxDistX-distX[k-1];
        DBG(2,"found minimum %g at %d\n",v,k);
        /*
         * here 'i' and 'j' stand for 'c1' and 'c2', respectively
         * (i,j) @ k = n*(i-1) - i*(i-1)/2 + j-i
         * n*(i-1) - i*(i-1)/2 + j-i = k
         * n*i - n - i^2/2 + i/2 + 1 <= k
         * -i^2/2 + i*(n+1/2) - n - k + 1 <= 0
         * i^2 - i*2*(n+1/2) + 2*(n + k - 1) >= 0
         * D = [2*(n+1/2)]^2 - 4*2*(n+k-1) =
         *     4 * (n^2 + n + 1/4 - 2*n - 2*(k-1)) =
         *     4 * (n^2 - n + 1/4 - 2*(k-1)) =
         * i0 = {2*(n+1/2) - 2*sqrt[(n^2 - n + 1/4 - 2*(k-1))]} / 2
         * i0 = (n+1/2) - sqrt[(n+1/2)^2 - 2*(n+k-1)]
         # i = floor (i0)
         *  if j=i+1, then i=i0 (no `floor' is needed)
         *  if j>i+1, then i=floor(i0), and the result of the floor()
         * operation can't be greater than i, as the first `k*' for
         * which `floor(i0*) > i' is the `k*' corresponding to `i* =
         * i+1, j*=i*+1' (i.e. the `i*' correspoding to next following
         * column in the distance matrix), so for all `k < k*' it holds
         * that `i0 < i*' and therefore `i = floor(i0)'.
         * j = (i-1)*(i/2-n)+i+k
         */
        // we are merging clusters `c1' and `c2' into a new one
        c1=floor(clusterCount+.5-sqrt(((long)clusterCount)*clusterCount-clusterCount+.25-2l*(k-1)));
        // note the '(long)x' above not to let the 'int*int' overflow
        c2=((long)c1-1)*(c1*.5-clusterCount)+c1+k;
        // note the '(long)x' above not to let the 'int*int' overflow
        // make c1 and c2 0-based
        c1=R2C(c1);
        c2=R2C(c2);
        DBG(2,"c1=%d, c2=%d\n",C2R(c1),C2R(c2));
        ASSERT(c1>=0,"invalid c1");
        ASSERT(c2>=0,"invalid c2");
        //R: if (dbg>0) cat(sprintf('Cluster %d: depth %g, merged clusters %d and %d (%s and %s).\n',s+n,v,clusterId[c1],clusterId[c2],printMembers(members[[c1]]),printMembers(members[[c2]])))
        DBG(1,"Cluster %d: depth %g, merged clusters %d and %d ((%s) and (%s)).\n",C2R(s+n),v,C2R(clusterId[c1]),C2R(clusterId[c2]),
            getMembers(members[c1],clusterSize[c1],strBuf),
            getMembers(members[c2],clusterSize[c2],strBufShort));
        //R: merging[s,1:3]<-c(clusterId[i],clusterId[c2],v);
        merging[s]=C2R(clusterId[c1]);
        merging[s+(n-1)]=C2R(clusterId[c2]);
        height[s]=v;

        // the other clusters just need to be updated in respect to
        // distance to the newly created cluster
        //R: otherClusters1<-mySeq(1,c1-1)
        //R: otherClusters2<-mySeq(c1+1,c2-1)
        //R: otherClusters3<-mySeq(c2+1,clusterCount)
        //R: otherClusters<-c(otherClusters1,otherClusters2,otherClusters3)
        for (i=0;i<c1;i++) otherClusters[i]=i;
        for (i=c1+1;i<c2;i++) otherClusters[i-1]=i;
        for (i=c2+1;i<clusterCount;i++) otherClusters[i-2]=i;
        DBG(3,"clusterCount: %d\n",clusterCount);
        for (i=0;i<clusterCount-2;i++) otherClustersTmp[i]=C2R(otherClusters[i]);
        DBG_CODE(3,printNumMatrix("otherClusters",otherClustersTmp,1,clusterCount-2));

        // get all samples constituting the merged clusters
        //R: xij<-x[c(members[[c1]],members[[c2]]),,drop=FALSE]
        c1n=clusterSize[c1];
        c2n=clusterSize[c2];
        DBG(2,"c1n=%d, c2n=%d\n",c1n,c2n);
        ASSERT(c1n>0,"invalid c1n");
        ASSERT(c2n>0,"invalid c2n");
        xijN=c1n+c2n;
        c1indices=members[c1];
        c2indices=members[c2];
        DBG(3,"constructing xij\n");
        //R: xij<-x[c(members[[c1]],members[[c2]]),,drop=FALSE]
        DBG_CODE(3,printNumMatrix("c1indices",c1indices,1,c1n));
        for (offset1=offset2=i=0;i<p;i++) {
            // offset1=xijN*i;
            // offset2=n*i;
            for (j=0;j<c1n;j++) {
                xij[j+offset1]=x[c1indices[j]+offset2];
            }
            offset1+=xijN;
            offset2+=n;
        }
        DBG_CODE(3,printNumMatrix("c2indices",c2indices,1,c2n));
        for (offset1=offset2=i=0;i<p;i++) {
            // offset1=xijN*i;
            // offset2=n*i;
            for (j=0;j<c2n;j++) {
                xij[c1n+j+offset1]=x[c2indices[j]+offset2];
            }
            offset1+=xijN;
            offset2+=n;
        }
        DBG_CODE(2,printDoubleMatrix("xij",xij,xijN,p));

        // compute the weight factor controlling the Mahalanobis-Euclidean
        // balance (to be applied when measuring distances relatively to
        // the merged cluster)
        if (thresh > 0) {
            DBG_CODE(2,{
                DBGU("clusterSize[c1] %d\n",clusterSize[c1]);
                DBGU("clusterSize[c2] %d\n",clusterSize[c2]);
                DBGU("(clusterSize[c1] + clusterSize[c2]) / ( nFull * thresh ) %.3f\n",
                     (clusterSize[c1] + clusterSize[c2]) / ( nFull * thresh ));
                DBGU("nFull %d\n",nFull);
                DBGU("thresh %.3f\n",thresh );
            });
            wf1=fmin2(1.0, (clusterSize[c1] + clusterSize[c2]) / ( nFull * thresh ) );
        } else {
            wf1=0;
        }
        // update fullMahalClusterCount if necessary
        if (wf1==1) {
            if (weightFactor[c1] + weightFactor[c2] == 2) {
                fullMahalClusterCount = fullMahalClusterCount - 1;
            } else if (weightFactor[c1] < 1 && weightFactor[c2] < 1) {
                fullMahalClusterCount = fullMahalClusterCount + 1;
            }
        }
        weightFactor[c1]=wf1;
        DBG(2,"wf1: %g\n",wf1);
        // try to fit an ellipsoid to the merged cluster - try to
        // compute the inverse of covariance matrix
        cov(xij,xijN,p,covXij,covMeansTmp);
        DBG_CODE(2,printDoubleMatrix("covXij",covXij,p,p));
        // if the cluster consists of few members only, do not take its shape too serious:
        // round it somehow (make closer to circle) by weighting
        if (wf1<1) {
            for (i=0;i<covXijLength;i++) {
                covXij[i] = wf1 * covXij[i] + (1 - wf1) * fakeInvCov[i];
            }
        }
        DBG_CODE(2,printDoubleMatrix("updated covXij",covXij,p,p));
        // try to compute Cholesky decomposition of the covariance matrix
        // if it fails, fall back to the unit matrix
        //R: c.cholDecomp<-tryCatch(chol(covXij),error=function(e) fakeInvCov)
        F77_CALL(dpotrf)("L",// Lower diagonal is used
                         &p,covXij,&p,&info);
        DBG(4,"info: %d\n",info);
        ASSERT(info>=0,"invalid dpotrf retcode");
        DBG_CODE(4,printDoubleMatrix("L",covXij,p,p));
        if (info==0) {
            double *cholDecomp=covXij;
            //R: detSqrt<-(1/prod(diag(c.cholDecomp)))^(2/p)
            detSqrt=1.0;
            for (offset=i=0;i<p;i++) {
                // detSqrt*=cholDecomp[(p+1)*i];
                detSqrt*=cholDecomp[offset];
                offset+=p+1;
            }
            detSqrt=pow(detSqrt,-2.0/p);
            //R: invcov_merged<-chol2inv(c.cholDecomp)
            F77_CALL(dpotri)("L",// Lower diagonal is used
                             &p,covXij,&p,&info);
            DBG(4,"info: %d\n",info);
            DBG_CODE(4,printDoubleMatrix("I",covXij,p,p));
            memcpy(invcovMerged,covXij,sizeof(*covXij)*p*p);
            // invcovMerged contains only lower triangular part (plus diagonal),
            // so initialize also the upper diagonal for debugging purposes
            DBG_CODE(3,{
                for (i=0;i<p;i++) {
                    for (j=i+1;j<p;j++) {
                        invcovMerged[i+p*j]=invcovMerged[j+p*i];
                    }
                }
            });
        } else {
            detSqrt=1.0;
            memcpy(invcovMerged,fakeInvCov,sizeof(*fakeInvCov)*p*p);
        }

        DBG_CODE(3,{
            printDoubleMatrix("detSqrt",&detSqrt,1,1);
            printDoubleMatrix("invcovMerged",invcovMerged,p,p);
        });
        // compute a new center of the merged cluster
        //R: centroid[c1,]<-(clusterSize[c1]*centroid[c1,,drop=FALSE] + clusterSize[c2]*centroid[c2,,drop=FALSE])/
        //R:      (clusterSize[c1] + clusterSize[c2]) */
        DBG_CODE(2,{
            sprintf(strBuf,"centroid[c1=%d]",C2R(c1));
            printDoubleMatrixRow(strBuf,centroid,clusterCount,p,c1);
            sprintf(strBuf,"centroid[c2=%d]",C2R(c2));
            printDoubleMatrixRow(strBuf,centroid,clusterCount,p,c2);
        });
        clusterSizeI=clusterSize[c1];
        clusterSizeJ=clusterSize[c2];
        for (offset=i=0;i<p;i++) {
            // offset=clusterCount*i;
            centroid[c1+offset]=(clusterSizeI*centroid[c1+offset] + clusterSizeJ*centroid[c2+offset])/
                (clusterSizeI+clusterSizeJ);
            offset+=clusterCount;
        }
        DBG_CODE(2,{
            sprintf(strBuf,"updated centroid[c1=%d]",C2R(c1));
            printDoubleMatrixRow(strBuf,centroid,clusterCount,p,c1);
        });
        
        ///////////////////////////////////////////////////////////////
        // update distX if we haven't reached the point at which we switch
        // to the full Mahalanobis style or if we haven't clustered all
        // the samples in the apriori clusters (otherwise we recompute
        // the whole distX later on)
        ///////////////////////////////////////////////////////////////
        if (fullMahalClusterCount < clusterCount-1 || // we haven't reached the point of switch or
            switchedToFullMahal || // we have already switched or
            n-(clusterCount-1) < gMergingCount) { // we've not clustered all the samples in the apriori clusters
            // ( n-(clusterCount-1) is the number of samples clustered so far)

            double minNewDistX=0,tmp;

            //R: ic1<-invcovMerged
            //R: detSqrtIc1<-detSqrt
            memcpy(ic1,invcovMerged,sizeof(*invcovMerged)*p*p);
            DBG_CODE(3,printDoubleMatrix("ic1",ic1,p,p));
            if (normalize || fullMahalClusterCount < clusterCount-1) {
                DBG(3," normalizing (normalize %d, clusters with full Mahalanobis = %d, clusters =  %d)\n",
                    normalize,fullMahalClusterCount,clusterCount);
                for (i=0;i<p*p;i++) {
                    ic1[i]/=detSqrt;
                }
                DBG_CODE(3,printDoubleMatrix("normalized ic1",ic1,p,p));
            }

            //R: for (ii in seq(along=along=otherClusters)) {
            for (oi=0;oi<clusterCount-2;oi++) {
                int otherCluster=otherClusters[oi];
                int xc1MemberCount,xc2MemberCount;
                int ii;

                DBG(3,"oi: %d\n",C2R(oi));
                DBG(3,"otherCluster: %d\n",C2R(otherCluster));

                DBG_CODE(4,{
                    printDoubleMatrix("normalized ic1",ic1,p,p);
                    for (ii=0;ii<clusterCount;ii++) {
                        DBGU("MEMBERS[%d]: ",ii);
                        DBGU("%s\n",getMembers(members[ii],clusterSize[ii],strBuf));
                    }
                });
                // compute the distance from the newly merged cluster c1+c2 to cluster otherClusters(oi)

                //R: iRelDistXIdx<-c(otherClusters1*(clusterCount-(otherClusters1+1)/2)-clusterCount+c1,
                //R:     c1*(clusterCount-(c1+1)/2)-clusterCount+otherClusters2,
                //R:     c1*(clusterCount-(c1+1)/2)-clusterCount+otherClusters3)
                if (oi<c1) {
                    distXIdx=R2C(((long)C2R(otherCluster-1))*clusterCount-(((long)C2R(otherCluster))*(C2R(otherCluster)+1))/2+C2R(c1));
                    // note the '(long)x' above not to let the 'int*int' overflow
                } else if (oi<c2-1) {
                    distXIdx=R2C(((long)C2R(c1-1))*clusterCount-(((long)C2R(c1))*(C2R(c1)+1))/2+C2R(otherCluster));
                    // note the '(long)x' above not to let the 'int*int' overflow
                } else {
                    distXIdx=R2C(((long)C2R(c1-1))*clusterCount-(((long)C2R(c1))*(C2R(c1)+1))/2+C2R(otherCluster));
                    // note the '(long)x' above not to let the 'int*int' overflow
                }
                DBG(4,"distXIdx: %ld\n",(long)distXIdx);
                ASSERT(distXIdx>=0,"invalid distXIdx");

                // if we are still clustering samples of the apriori clusters, the distance
                // between clusters in disctinct apriori clusters arem by definition, se to Inf
                //R: if (!is.null(g) && g[members[[i]][1]]!=g[members[[otherClusters[ii]]][1]] && # clusters in distinct apriori clusters
                //R:    pointCount-(clusterCount-1) < gMergingCount) { # and we are not done with clustering the samples from prior clusters
                //R:    distX[iRelDistXIdx[ii]]<-Inf
                //R: }
                if (!IS_R_NULL(G) && g[*members[c1]]!=g[*members[otherCluster]] && // clusters in distinct apriori clusters
                    n-(clusterCount-1) < gMergingCount) { // and we are not done with clustering the samples from prior clusters
                    distX[distXIdx]=0; // set the value meaning Inf (as distX is inverted and we find the maximum absolute value,
                    // 0 denote value that can't be found if there are still some unclustered samples)
                } else {

                    if (quick) {
                        //R: xc1<-centroid[otherClusters[oi],,drop=FALSE]
                        xc1MemberCount=1;
                        for (offset=i=0;i<p;i++) {
                            // xc1[i]=centroid[otherCluster+clusterCount*i];
                            xc1[i]=centroid[otherCluster+offset];
                            offset+=clusterCount;
                        }
                    } else {
                        //R: xc1<-x[members[[otherClusters[oi]]],,drop=FALSE]
                        xc1MemberCount=clusterSize[otherCluster];
                        DBG(3,"xc1MemberCount: %d\n",xc1MemberCount);
                        ASSERT(xc1MemberCount>0,"invalid xc1MemberCount");
                        for (offset1=offset2=i=0;i<p;i++) {
                            // offset1=xc1MemberCount*i;
                            // offset2=n*i;
                            Num cumI=0;
                            Num *mmbrs=members[otherCluster];
                            for (j=0;j<xc1MemberCount;j++,cumI++) {
                                xc1[cumI+offset1]=x[mmbrs[j]+offset2];
                            }
                            offset1+=xc1MemberCount;
                            offset2+=n;
                        }
                    }
                    DBG_CODE(4,{
                        printDoubleMatrix("xc1",xc1,xc1MemberCount,p);
                        sprintf(strBuf,"centroid[c1=%d]:",C2R(c1));
                        printDoubleMatrixRow(strBuf,centroid,clusterCount,p,c1);
                    });
                    //R: xc1<-xc1-matrix(centroid[c1,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE)
                    for (offset1=offset2=i=0;i<p;i++) {
                        // offset1=xc1MemberCount*i;
                        // offset2=clusterCount*i;
                        for (j=0;j<xc1MemberCount;j++) {
                            xc1[j+offset1]-=centroid[c1+offset2];
                        }
                        offset1+=xc1MemberCount;
                        offset2+=clusterCount;
                    }
                    DBG_CODE(4,printDoubleMatrix("centered xc1",xc1,xc1MemberCount,p));
                    // distMaha1 holds a vector of squares of mahalanobis distances:
                    // mean Mahalanobis distance from c1+c2 to some other cluster
                    //R: distMaha1<-mean(sqrt(rowSums((xc1%*%ic1)*xc1)))
                    distMaha1=computeMahalanobisDistance(xc1,xc1MemberCount,p,ic1,distMahaBuf,dbg);
                    DBG(3,"distMaha1 %g\n",distMaha1);

                    // compute the distance from cluster otherClusters(oi) to the newly merged cluster c1+c2
                    if (quick) {
                        //R: xc2<-centroid[c1,,drop=FALSE]
                        xc2MemberCount=1;
                        for (offset=i=0;i<p;i++) {
                            // xc2[i]=centroid[c1+clusterCount*i];
                            xc2[i]=centroid[c1+offset];
                            offset+=clusterCount;
                        }
                    } else {
                        xc2MemberCount=xijN;
                        memcpy(xc2,xij,sizeof(*xij)*xijN*p);
                    }
                    DBG_CODE(4,{
                        printDoubleMatrix("xc2",xc2,xc2MemberCount,p);
                        sprintf(strBuf,"centroid[oi=%d]",C2R(otherCluster));
                        printDoubleMatrixRow(strBuf,centroid,clusterCount,p,otherCluster);
                    });

                    //R: xc2<-xc2-matrix(centroid[otherClusters[oi],,drop=FALSE],nrow(xc2),ncol(xc2),byrow=TRUE)
                    for (offset1=offset2=i=0;i<p;i++) {
                        // offset1=xc2MemberCount*i;
                        // offset2=clusterCount*i;
                        for (j=0;j<xc2MemberCount;j++) {
                            xc2[j+offset1]-=centroid[otherCluster+offset2];
                        }
                        offset1+=xc2MemberCount;
                        offset2+=clusterCount;
                    }
                    DBG_CODE(4,printDoubleMatrix("centered xc2",xc2,xc2MemberCount,p));

                    //R: ic2<-invcov[[otherClusters[ii]]]
                    memcpy(ic2,invcov+p*p*otherCluster,p*p*sizeof(*ic2));

                    DBG_CODE(4,printDoubleMatrix("ic2",ic2,p,p));

                    if (normalize || fullMahalClusterCount < clusterCount-1) {
                        DBG(4," normalizing (normalize %d, clusters with full Mahalanobis = %d, clusters =  %d)\n",
                            normalize,fullMahalClusterCount,clusterCount);
                        double detSqrt2=detsSqrt[otherCluster];
                        for (i=0;i<p*p;i++) {
                            ic2[i]/=detSqrt2;
                        }
                        DBG_CODE(4,printDoubleMatrix("normalized ic2",ic2,p,p));
                    }

                    // mean Mahalanobis distance from otherCluster to c1+c2
                    //R: distMaha2<-mean(sqrt(rowSums((xc2%*%ic2)*xc2)))
                    distMaha2=computeMahalanobisDistance(xc2,xc2MemberCount,p,ic2,distMahaBuf,dbg);
                    DBG(3,"distMaha2: %g\n",distMaha2);

                    // merge the clusterId(c1) <-> clusterId(c2) distances
                    //R: distX[iRelDistXIdx[ii]]<-mean(c(distMaha1,distMaha2))
                    tmp=distX[distXIdx]=maxDistX-(distMaha1+distMaha2)/2;
                    if (tmp<minNewDistX) minNewDistX=tmp;
                }

                DBG(3,"otherCluster: %d\n",C2R(otherCluster));
                //DBG("Dist from %d=(%s)\n",C2R(clusterId[otherCluster]),printMembers(members[otherCluster]));
                DBG(3,"Dist from %d=(%s) to %d=(%s %s): %g.\n",
                    C2R(clusterId[otherCluster]),getMembers(members[otherCluster],clusterSize[otherCluster],strBuf),
                    C2R(s+n),
                    getMembers(members[c1],clusterSize[c1],strBufShort),
                    getMembers(members[c2],clusterSize[c2],strBufShort2),
                    maxDistX-distX[distXIdx]);
            }
            // move distX values if some distance was below 0 (in that case
            // finding the maximum absolute value (of maxDistX-dist) would fail
            // to get the minimum distance)
            if (minNewDistX<0) {
                DBG(4,"moving distX by %f, old maxDistX %f\n",minNewDistX,maxDistX);
                DBG(4,"distX[0]: %f\n",distX[0]);
                for (i=0;i<distXLen;i++) distX[i]-=minNewDistX;
                maxDistX-=minNewDistX;
                DBG(4,"new maxDistX %f\n",maxDistX);
                DBG(4,"new distX[0] %f\n",distX[0]);
            }
        }

        // clusters clusterId[c1] and clusterId[c2] merged, remove
        // clusterId[c2]-related info, put info about the newly created
        // cluster at position occupied by clusterId[c1] previously
        //DBG(2,"  updating...\n");
        //R: members[[c1]]<-c(members[[c1]],members[[c2]])
        n1=clusterSize[c1];
        n2=clusterSize[c2];
        DBG(2,"  members: n1=%d, n2=%d\n",n1,n2);
        ASSERT(n1>0,"invalid n1");
        ASSERT(n2>0,"invalid n2");
        DBG(5,"  membersPoolPos: %lu\n",membersPoolPos);
        Num *newMembers=membersPool+membersPoolPos;
        membersPoolPos+=n1+n2;
        ASSERT(membersPoolPos<=membersPoolSize,"membersPool exhausted");
        DBG(5,"  copying members of cluster1\n");
        for (i=0;i<n1;i++) {
            newMembers[i]=members[c1][i];
        }
        DBG(5,"  copying members of cluster2\n");
        for (i=0;i<n2;i++) {
            newMembers[n1+i]=members[c2][i];
        }
        DBG(5,"  updating members\n");
        members[c1]=newMembers;
        //R: members<-members[-c2]
        memmove(members+c2,members+c2+1,sizeof(*members)*(clusterCount-1-c2));
        DBG(5,"  members updated\n");
        //R: invcov<-invcov[-c2]
        memmove(invcov+p*p*c2,invcov+p*p*(c2+1),sizeof(*invcov)*p*p*(clusterCount-1-c2));
        //R: invcov[[c1]]<-invcov_merged
        memcpy(invcov+c1*p*p,invcovMerged,sizeof(*invcov)*p*p);
        //R: detsSqrt<-detsSqrt[-c2]
        memmove(detsSqrt+c2,detsSqrt+c2+1,sizeof(*detsSqrt)*(clusterCount-1-c2));
        //R: detsSqrt[i]<-detSqrt
        detsSqrt[c1]=detSqrt;
        //R: centroid<-centroid[-c2,]
        for (i=0;i<p-1;i++) memmove(centroid+c2+i*(clusterCount-1),centroid+c2+i*clusterCount+1,sizeof(*centroid)*(clusterCount-1));
        memmove(centroid+c2+(p-1)*(clusterCount-1),centroid+c2+(p-1)*clusterCount+1,sizeof(*centroid)*(clusterCount-1-c2));
        //R: weightFactor<-weightFactor[-c2]
        memmove(weightFactor+c2,weightFactor+c2+1,sizeof(*weightFactor)*(clusterCount-1-c2));
        //R: c2RelDistXIdx<-c(otherClusters1*(clusterCount-(otherClusters1+1)/2)-clusterCount+c2,
        //R:    otherClusters2*(clusterCount-(otherClusters2+1)/2)-clusterCount+c2,
        //R:     c2*(clusterCount-(c2+1)/2)-clusterCount+otherClusters3,
        //R:     c1*(clusterCount-(c1+1)/2)-clusterCount+c2)
        //R: distX<-distX[-c2RelDistXIdx]
        /* the following formulas are simply optimized versions of the indices into the distance matrix,
         * consider that
         *   i*(2*n-(i+1))/2-n+j
         * equals
         *   n*i - (i*(i+1))/2 - n + j
         *   n*(i-1) - (i*(i-1+2))/2 + j
         *   n*(i-1) - (i*(i-1))/2 + j-i
         */
        DBG(5,"  updating distX\n");
        // note the use of '2l' below not to let the 'int*int' overflow
        for (i=0;i<c1;i++) distXIndices[i]=R2C(C2R(i)*(2l*clusterCount-(C2R(i)+1))/2-clusterCount+C2R(c2));
        distXIndices[c1]=R2C(C2R(c1)*(2l*clusterCount-(C2R(c1)+1))/2-clusterCount+C2R(c2));
        for (i=c1+1;i<c2;i++) distXIndices[i]=R2C(C2R(i)*(2l*clusterCount-(C2R(i)+1))/2-clusterCount+C2R(c2));
        for (i=c2+1;i<clusterCount;i++) distXIndices[i-1]=R2C(C2R(c2)*(2l*clusterCount-(C2R(c2)+1))/2-clusterCount+C2R(i));
        distXIndicesLen=clusterCount-1;
        DBG_CODE(4,{
            DBGU("c1=%d, c2=%d, clusterCount=%d\n",c1,c2,clusterCount);
            printLongNumMatrix("distXIndices part 1 (... c1)",distXIndices,1,c1);
            printLongNumMatrix("distXIndices part 2 (c1)",distXIndices+c1,1,1);
            printLongNumMatrix("distXIndices part 3 (c1 ... c2)",distXIndices+c1+1,1,c2-c1-1);
            printLongNumMatrix("distXIndices part 4 (c2 ...)",distXIndices+c2,1,clusterCount-c2-1);
            printLongNumMatrix("distXIndices",distXIndices,1,distXIndicesLen);
            printDistMatrix("distX pre update",distX,clusterCount,maxDistX);
        });
        /*
         * source:  0 1 2 3 4 5 6 7 8 9 A B C
         * remove:      2,    5,  7,8,    B C
         * keep:    0-1   3-4   6     9-A
         *
         * target pos=2
         *
         * remove 2: copy from pos 2+1 to target pos=2 in length of 5-2-1 = 2:
         *       0 1 2 3 4 5 6 7 8 9 A B C
         *           ^ < <
         *       0 1 3 4 4 5 6 7 8 9 A B C
         *               ^ new target pos
         *
         * remove 5: copy from pos 5+1 to target pos=4 in length of 7-5-1 = 1:
         *       0 1 3 4 4 5 6 7 8 9 A B C
         *               ^   <
         *       0 1 3 4 6 5 6 7 8 9 A B C
         *                 ^ new target pos
         *
         * remove 7,8: copy from pos 8+1 to target pos=5 in length of B-8-1 = 2:
         *       0 1 3 4 6 5 6 7 8 9 A B C
         *                 ^       < <
         *       0 1 3 4 6 9 A 7 8 9 A B C
         *                     ^ new target pos
         *
         * remove B,C: no need to copy anything
         *       0 1 3 4 6 9 A 7 8 9 A B C
         *                     ^  target pos
         */
        // length of distXIndices is clusterCount-1,
        // subtract another 1 not to stop at the last but one element of distXIndices
        target=distXIndices[0];
        i=0;
        while (i<distXIndicesLen) {
            int ii=i;
            while (ii+1<distXIndicesLen && distXIndices[ii+1]==distXIndices[ii]+1) ii++;
            // distXIndices[i] is the index of the first entry to remove,
            // distXIndices[ii] is the index of the last entry to remove,
            // distXIndices[ii+1] (if exists) is the index of the next first entry to remove
            if (distXIndices[ii]==distXLen-1) break; // no more entries to move
            source=distXIndices[ii]+1;
            if (ii+1<distXIndicesLen) {
                count=distXIndices[ii+1]-distXIndices[ii]-1;
                DBG(4," distXIndices[i=%ld]=%ld, distXIndices[ii=%ld]=%ld, distXIndices[ii+1=%ld]=%ld\n",
                    i,distXIndices[i],ii,distXIndices[ii],ii+1,distXIndices[ii+1]);
            } else {
                count=distXLen-distXIndices[ii]-1;
                DBG(4," distXIndices[i=%ld]=%ld, distXIndices[ii=%ld]=%ld\n",
                    i,distXIndices[i],ii,distXIndices[ii]);
            }
            DBG_CODE(5,{
                for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
                printDoubleMatrix("distX pre move",distXTmp,1,((long)clusterCount-1)*(clusterCount-2)/2);
            });
            DBG(4,"memmove target %ld (%p, value %g), source %ld (%p, value %g), len %ld (%ld)\n",
                target,distX+target,maxDistX-distX[target],source,distX+source,maxDistX-distX[source],count,sizeof(*distX)*count);
            memmove(distX+target,distX+source,sizeof(*distX)*count);
            DBG_CODE(5,{
                DBGU("target %ld (%p, value %g)\n",target,distX+target,maxDistX-distX[target]);
                for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
                printDoubleMatrix("distX post move",distXTmp,1,((long)clusterCount-1)*(clusterCount-2)/2);
            });
            target+=count;
            i=ii+1;
        }
        distXLen-=clusterCount-1;
        DBG_CODE(3,{
            DBGU("clusterCount: %d\n",clusterCount);
            DBGU("distXLen: %d\n",distXLen);
            for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
            printDoubleMatrix("distX",distXTmp,1,((long)clusterCount-1)*(clusterCount-2)/2);
        });
        DBG_CODE(3,printDistMatrix("distX",distX,clusterCount-1,maxDistX));
        DBG(5,"  distX updated\n");
        // update clusterCount, clusterSize, clusterId
        clusterCount--;
        clusterSize[c1]=clusterSize[c1]+clusterSize[c2];
        //R: clusterSize<-clusterSize[-c2]
        memmove(clusterSize+c2,clusterSize+c2+1,sizeof(*clusterSize)*(clusterCount-c2));
        DBG_CODE(4,printNumMatrix("clusterSize",clusterSize,1,clusterCount));
        clusterId[c1]=n+s;
        //R: clusterId<-clusterId[-c2]
        memmove(clusterId+c2,clusterId+c2+1,sizeof(*clusterId)*(clusterCount-c2));
        DBG_CODE(4,printNumMatrix("clusterId",clusterId,1,clusterCount));

        // recompute distX in case we reached the point at which we switch to the full Mahalanobis style
        // and we've also clustered all the samples in the apriori clusters
        if (fullMahalClusterCount == clusterCount && // we reached the switch point and
            !switchedToFullMahal && // we haven't swicthed yet and
            gMergingCount <= n-clusterCount) { // we have clustered all samples in the apriori clusters
            int idx,i1,i2;
            double tmp,minNewDistX=0;

            switchedToFullMahal=1;

            DBG(2,"Recomputing all distances.\n");
            idx=0;

            for (i1=0;i1<clusterCount-1;i1++) {

                Num xc1n=constructRepresentantsOfCluster(x,n,p,centroid,i1,
                                                         members[i1],clusterSize[i1],
                                                         xc1,quick,dbg);
                DBG_CODE(3,printDoubleMatrix("xc1Orig",xc1,xc1n,p));

                memcpy(ic2,invcov+i1*p*p,sizeof(*invcov)*p*p);
                DBG_CODE(3,printDoubleMatrix("ic2",ic2,p,p));
                if (normalize) {
                    DBG(3," normalizing ic2\n");
                    detSqrt=detsSqrt[i1];
                    for (i=0;i<p*p;i++) {
                        ic2[i]/=detSqrt;
                    }
                    DBG_CODE(3,printDoubleMatrix("normalized ic2",ic2,p,p));
                }

                for (i2=i1+1;i2<clusterCount;i2++) {
                    DBG(2,"recomputing dist from %d to %d\n",C2R(i1),C2R(i2));

                    // compute the distance from cluster (i1) to (i2)
                    constructBetweenClusterDistanceMatrixFromXi(xc1,xc1n,p,centroid,clusterCount,i2,y,quick,dbg);
                    DBG_CODE(3,printDoubleMatrix("xc1",y,xc1n,p));

                    memcpy(ic1,invcov+i2*p*p,sizeof(*invcov)*p*p);
                    DBG_CODE(3,printDoubleMatrix("ic1",ic1,p,p));
                    if (normalize) {
                        DBG(3," normalizing ic1\n");
                        detSqrt=detsSqrt[i2];
                        for (i=0;i<p*p;i++) {
                            ic1[i]/=detSqrt;
                        }
                        DBG_CODE(3,printDoubleMatrix("normalized ic1",ic1,p,p));
                    }
                    distMaha1=computeMahalanobisDistance(y,xc1n,p,ic1,distMahaBuf,dbg);
                    DBG(2,"distMaha1: %g\n",distMaha1);

                    // compute the distance from cluster (i2) to (i1)
                    yn=constructBetweenClusterDistanceMatrix(x,n,p,centroid,clusterCount,
                                                             i2,i1,members[i2],clusterSize[i2],y,
                                                             quick,dbg);
                    DBG_CODE(3,printDoubleMatrix("xc2",y,yn,p));
                    distMaha2=computeMahalanobisDistance(y,yn,p,ic2,distMahaBuf,dbg);
                    DBG(2,"distMaha2: %g\n",distMaha2);

                    // merge the clusterId(i1) <-> clusterId(i2) distances
                    tmp=distX[idx]=maxDistX-(distMaha1+distMaha2)/2;
                    if (tmp<minNewDistX) minNewDistX=tmp;
                    idx++;
                }
            }

            // move distX values if some distance was below 0 (in that case
            // finding the maximum absolute value (of maxDistX-dist) would fail
            // to get the minimum distance)
            if (minNewDistX<0) {
                DBG(4,"moving distX by %f, old maxDistX %f\n",minNewDistX,maxDistX);
                DBG(4,"distX[0] %f\n",distX[0]);
                for (i=0;i<distXLen;i++) distX[i]-=minNewDistX;
                maxDistX-=minNewDistX;
                DBG(4,"new maxDistX %f\n",maxDistX);
                DBG(4,"new distX[0] %f\n",distX[0]);
            }
        }

        R_CheckUserInterrupt();
    }

    DBG(2,"mhclust_ cleaning\n");

    MEM_FREE(distXTmp);
    MEM_FREE(strBuf);
    MEM_FREE(strBufShort);
    MEM_FREE(strBufShort2);
    MEM_FREE(covMeansTmp);
    MEM_FREE(clusterSize);
    MEM_FREE(membersPool);
    MEM_FREE(members);
    MEM_FREE(clusterId);
    MEM_FREE(xij);
    MEM_FREE(covXij);
    MEM_FREE(fakeInvCov);
    MEM_FREE(ic1);
    MEM_FREE(ic2);
    MEM_FREE(otherClusters);
    MEM_FREE(otherClustersTmp);
    MEM_FREE(distXIndices);
    MEM_FREE(y);
    MEM_FREE(distMahaBuf);
    MEM_FREE(invcov);
    MEM_FREE(invcovMerged);
    MEM_FREE(detsSqrt);
    MEM_FREE(centroid);
    MEM_FREE(xc1);
    MEM_FREE(xc2);
    MEM_FREE(weightFactor);

    UNPROTECT(nprot);

    DBG(2,"mhclust_ finishes\n");

    return(R_NilValue);
}
