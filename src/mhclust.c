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

#if 01
// debug print (unconditional)
# define DBGU(...) Rprintf(__VA_ARGS__)
// debug code run conditionally on the value of dbgLevel
# define DBG_CODE(dbgLevel,...) if (dbg>=dbgLevel) __VA_ARGS__
// debug print (conditional on the value of dbgLevel)
# define DBG(dbgLevel,...) DBG_CODE(dbgLevel,DBGU(__VA_ARGS__))
#else // debug completely disabled
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

// TODO:optimize if necessary
#define MEM_ALLOC_FAST MEM_ALLOC
#define MEM_FREE_FAST MEM_FREE

typedef int Num;

typedef struct _NumList {
    Num *data;
    Num n;
    Num capacity;
} NumList;
NumList allocateNumList(Num capacity) {
    NumList list;
    list.data=(Num*)MEM_ALLOC(capacity,sizeof(Num));
    list.capacity=capacity;
    list.n=0;
    return list;
}
void deallocateNumList(NumList list) {
    if (list.data!=NULL) MEM_FREE(list.data);
}

/*
 * Print members of a cluster, a sequence of observation IDs.
 #
 * TODO: the return value should be freed if the memory allocation changes
 */
const char *printMembers(NumList list) {
    char *buf=MEM_ALLOC(20*list.n,sizeof(char));
    char tmp[20];
    Num i;
    Num pos=0;
    for (i=0;i<list.n;i++) {
        sprintf(tmp,"%d",C2R(list.data[i]));
        memcpy(buf+pos,tmp,strlen(tmp));
        pos+=strlen(tmp);
        if (i<list.n-1) {
            buf[pos++]=' ';
        }
    }
	buf[pos]=0;
    return buf;
}

/*
 * Print matrix of doubles.
 *
 * TODO: the return value should be freed if the memory allocation changes
 */
void printDoubleMatrix(const char *name,double *x,Num rows,Num cols) {
    //Rprintf("printDoubleMatrix(cols=%d)\n",cols);
    char *buf=MEM_ALLOC(100*cols,sizeof(char));
    char tmp[100];
    Num i,j;
    Rprintf("%s:\n",name);
    if (rows>0 && cols>0) {
        for (i=0;i<rows;i++) {
            Num pos=0;
            for (j=0;j<cols;j++) {
                sprintf(tmp,"%.4f%s",x[i+rows*j],j<cols-1?"  ":"");
                memcpy(buf+pos,tmp,strlen(tmp));
                pos+=strlen(tmp);
            }
            buf[pos]=0;
            Rprintf("    %s\n",buf);
        }
    }
}

/*
 * Print a row of a matrix of doubles.
 *
 * TODO: the return value should be freed if the memory allocation changes
 */
void printDoubleMatrixRow(const char *name,double *x,Num rows,Num cols,Num row) {
    char *buf=MEM_ALLOC(100*cols,sizeof(char));
    char tmp[100];
    Num j;
    Rprintf("%s:\n",name);
    if (rows>0 && cols>0) {
        Num pos=0;
        for (j=0;j<cols;j++) {
            sprintf(tmp,"%.3f%s",x[row+rows*j],j<cols-1?"  ":"");
            memcpy(buf+pos,tmp,strlen(tmp));
            pos+=strlen(tmp);
        }
        buf[pos]=0;
        Rprintf("    %s\n",buf);
    }
}
/*
 * Print distance matrix (the lower triangle annotated with indices).
 *
 * TODO: the return value should be freed if the memory allocation changes
 */
void printDistMatrix(const char *name,double *x,Num n,double inversionFactor) {
    char *buf=MEM_ALLOC(100*n,sizeof(char));
    char tmp[100];
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
            sprintf(tmp," %9d",C2R(i));
            memcpy(buf+pos,tmp,strlen(tmp));
            pos+=strlen(tmp);
        }
        buf[pos]=0;
        Rprintf("        %s\n",buf);

        // print data rows
        for (j=1;j<n;j++) { // cluster to, in rows
            pos=0;
            for (i=0;i<n-1;i++) { // cluster from, in columns
                if (i<j) {
                    sprintf(tmp," %9.6f",inversionFactor-x[R2C(n*(C2R(i)-1) - C2R(i)*(C2R(i)-1)/2 + C2R(j)-C2R(i))]);
                } else {
                    sprintf(tmp,"          ");
                }
                memcpy(buf+pos,tmp,strlen(tmp));
                pos+=strlen(tmp);
            }
            buf[pos]=0;
            // prepend and append row number
            Rprintf("    %3d %s    %3d \n",j+1,buf,j+1);
        }

        // print trailer
        pos=0;
        for (i=0;i<n-1;i++) {
            sprintf(tmp," %9d",i+1);
            memcpy(buf+pos,tmp,strlen(tmp));
            pos+=strlen(tmp);
        }
        buf[pos]=0;
        Rprintf("        %s\n",buf);
    }
}

/*
 * Print matrix of integers (Num's).
 *
 * TODO: the return value should be freed if the memory allocation changes
 */
void printNumMatrix(const char *name,Num *x,Num rows,Num cols) {
    char *buf=MEM_ALLOC(20*cols,sizeof(char));
    char tmp[20];
    Num i,j;
    Rprintf("%s:\n",name);
    if (cols>0) {
        for (i=0;i<rows;i++) {
            Num pos=0;
            for (j=0;j<cols;j++) {
                sprintf(tmp,"%d%s",x[i+rows*j],j<cols-1?"  ":"");
                memcpy(buf+pos,tmp,strlen(tmp));
                pos+=strlen(tmp);
            }
            buf[pos]=0;
            Rprintf("    %s\n",buf);
        }
    }
}

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
    Num i,i1,i2,j;
    double *means=covMeansTmp;
    for (i=0;i<c;i++) {
        double sum=0.0;
        for (j=0;j<r;j++) {
            sum+=x[j+i*r];
        }
        means[i]=sum/r;
    }
    for (i1=0;i1<c;i1++) {
        for (i2=0;i2<c;i2++) {
            if (i2>=i1) {
                /* computing lower diagonal */
                double sum=0.0;
                for (j=0;j<r;j++) {
                    sum+=(x[j+i1*r]-means[i1])*(x[j+i2*r]-means[i2]);
                }
                res[i1+i2*c]=sum/(r-1);
            } else {
                /* computing upper diagonal from lower diagonal */
                res[i1+i2*c]=res[i2+i1*c];
            }
        }
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
 *   dbg
 */
double computeMahalanobisDistance(double *X,Num n,Num p,double *IC,int dbg) {
    Num i,j;
    double *result=(double*)MEM_ALLOC_FAST(n*p,sizeof(double));
    double r;

    /*
     * dsymm(const char *side, const char *uplo, const int *m,const int *n, const double *alpha,
     *       const double *a, const int *lda,const double *b, const int *ldb,const double *beta,
     *       double *c, const int *ldc);
     * in general:
     *    alpha*B * A + beta*C  ->  C
     * we use only:
     *    X    *    IC    ->   R
     * {n x p} * {p x p}    {n x p}
     * i.e.
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
                    result,&n); // R, ld = n
    DBG_CODE(4,printDoubleMatrix("result",result,n,p));

    // mean Mahalanobis distance from origin to rows of X
    // R: mean(sqrt(rowSums((X%*%IC)*X)))
    for (i=0;i<n;i++) { // TODO: optimize?
        double sum=0.0;
        // sum rows
        for (j=0;j<p;j++) {
            sum+=result[i+n*j]*X[i+n*j];
        }
        // take the square root and store in the first column of "result"
        result[i]=sqrt(sum);
    }
    i=1; // lda of result
    r=F77_NAME(dasum)(&n,result,&i)/n;
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
 *  members - members of individual clusters (see the definition of
 *      `members' below)
 *  i=clusterFrom - index of cluster for which to compute the
 *      representants
 *  y - where to store the resulting matrix, of size `yn x p'
 *  ynPtr - pointer to hold `yn', the number of rows of the result
 *  quickMode - boolean flag determining whether to use all members of
 *      clusterFrom, or only its centroid
 * dbg
 */
Num constructRepresentantsOfCluster(double *x,Num n,Num p,
                                    double *centroid,NumList *members,
                                    Num cluster,double *y,int quickMode,int dbg) {
    Num i,j,yn;

    if (quickMode) {
        yn=1;
        for (i=0;i<p;i++) {
            y[i]=centroid[cluster+n*i];
        }
    } else {
        yn=members[cluster].n;
        Num *mmbrs=members[cluster].data;
        for (i=0;i<p;i++) {
            Num offset1=yn*i;
            Num offset2=n*i;
            Num cumI=0;
            for (j=0;j<yn;j++,cumI++) {
                y[cumI+offset1]=x[mmbrs[j]+offset2];
            }
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
 *  j=clusterTo - index of target cluster
 *  y - where to store the resulting matrix of size `n x p'
 *  quickMode - boolean flag determining whether to use all members of
 *      clusterFrom, or only its centroid
 * dbg
 */
void constructBetweenClusterDistanceMatrixFromXi(double *Xi,Num n,Num p,
                                                 double *centroid,Num clusterCount,
                                                 Num clusterTo,double *y,
                                                 int quickMode,int dbg) {
    Num i,j;

    DBG_CODE(4,printDoubleMatrix("y",Xi,n,p));
    DBG_CODE(4,printDoubleMatrixRow("centroid[clusterTo]",centroid,clusterCount,p,clusterTo));
    /* xc1<-xc1-matrix(centroid[c1,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE) */
    for (i=0;i<p;i++) {
        long offset1=n*i;
        long offset2=clusterCount*i;
        for (j=0;j<n;j++) {
            y[j+offset1]=Xi[j+offset1]-centroid[clusterTo+offset2];
        }
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
 *  members - members of individual clusters (see the definition of
 *      members)
 *  i=clusterFrom - index of source cluster
 *  j=clusterTo - index of target cluster
 *  y - where to store resulting the matrix of size `yn x p'
 *  ynPtr - pointer to hold `yn', the number of rows of the result
 *  quickMode - boolean flag determining whether to use all members of
 *      clusterFrom, or only its centroid
 * dbg
 */
Num constructBetweenClusterDistanceMatrix(double *x,Num n,Num p,
                                          double *centroid,Num clusterCount,
                                          NumList *members,
                                          Num clusterFrom,Num clusterTo,
                                          double *y,
                                          int quickMode,int dbg) {
    Num i,j,yn;
    yn=constructRepresentantsOfCluster(x,n,p,centroid,members,
                                       clusterFrom,y,quickMode,dbg);
    DBG_CODE(4,printDoubleMatrix("y",y,yn,p);
    DBG_CODE(4,printDoubleMatrixRow("centroid[clusterTo]",centroid,clusterCount,p,clusterTo));
    //R: xc1<-xc1-matrix(centroid[c1,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE)
    for (i=0;i<p;i++) {
        long offset1=yn*i;
        long offset2=clusterCount*i;
        for (j=0;j<yn;j++) {
            y[j+offset1]-=centroid[clusterTo+offset2];
        }
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
 */
SEXP mhclust_(SEXP X,SEXP DistX,SEXP Merging,SEXP Height,SEXP Thresh,SEXP Quick,SEXP Normalize,SEXP Verb) {
    int dbg;
    int quick;

    int nprot = 0;

    double *x; // input data
    Num *dimX;
    Num n,n1,n2; // number of data points
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
    Num distXIdx,*distXIndices; // indices into distX
    Num distXLen;

    Num *merging; // output argument: cluster merging
    double *height; // output argument: cluster heights

    double thresh;
    int normalize;

    Num *clusterSize;// TODO: holds the same info as members[i].n - unfity?
    Num *clusterId; // id of cluster contained at this specific position
    double detSqrt;
    double* detsSqrt; // square root of determinants of inverse covariance matrices

    double *weightFactor,wf1; // factors weighting inverse covariance matrix against fakeInvCov

    NumList *members;
    Num *otherClusters;
    Num *otherClustersTmp;
    Num c1,c2,s,i,j,k,oi;
    Num cap1,cap2;
    Num clusterSizeI,clusterSizeJ;
    Num c1n,c2n,*c1indices,*c2indices;

    Num fullMahalClusterCount; // number of clusters for which full Mahalanobis distance can be computed
    int switchedToFullMahal;

    double distMaha1,distMaha2;
    double *y; // temporal storage, used for `y' for `constructBetweenClusterDistanceMatrix'
    Num yn; // number of data points in `y'

    int lda;
    int info;

    char strBuffer[100];

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

    n=dimX[0]; // number of elementary points subject to clustering (merging)
    DBG(2,"n=%d\n",n);
    if (distXLen!=n*(n-1)/2) error("invalid distXLen");
    p=dimX[1]; // number of dimensions of the feature space
    DBG(2,"p=%d\n",p);
    covMeansTmp=(double *)MEM_ALLOC(p,sizeof(double));
    DBG_CODE(3,printDoubleMatrix("x",x,n,p));
    covXijLength=p*p;

    DBG(2,"length of distX: %d\n",distXLen);
    DBG_CODE(3,printDoubleMatrix("distX",distX,1,distXLen));
    // convert distX to C-distX, such that we could find the maximum absolute value
    // instead of the minimum value
    lda=1;
    k=F77_CALL(idamax)(&distXLen,distX,&lda);
    ASSERT(k>0,"invalid idamax retcode");
    maxDistX=distX[k-1];
    DBG(4,"found maximum %g at %d\n",maxDistX,k-1));
    DBG_CODE(4,printDoubleMatrix("distX (orig)",distX,1,distXLen));
    for (i=0;i<distXLen;i++) distX[i]=maxDistX-distX[i];
    DBG_CODE(4,printDoubleMatrix("distX (tx'd)",distX,1,distXLen));

    // how may clusters (hot merged together yet) remain; in the
    // beginning, all points are clusters, in the end only one huge
    // cluster exists
    clusterCount=n;
    // number of elementary points in each cluster
    clusterSize=(int *)MEM_ALLOC(2*clusterCount-1,sizeof(int));
    for (i=0;i<clusterCount;i++) clusterSize[i]=1;
    for (i=0;i<clusterCount-1;i++) clusterSize[clusterCount+i]=0;

    // clusters being made (by merging two smaller clusters) are
    // assigned unique IDs, but reside in data structured indexed by
    // index of one of its subclusters - thus we need to map the
    // 1:n space into IDs of current clusters
    clusterId=(int *)MEM_ALLOC(n,sizeof(int));
    for (i=0;i<n;i++) clusterId[i]=i;

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

    distXIndices=(int*)MEM_ALLOC(clusterCount-1,sizeof(int));

    y=(double *)MEM_ALLOC(n*p,sizeof(double));

    // members (elementary observations) of each cluster
    members=(NumList *)MEM_ALLOC(n,sizeof(NumList));
    for (i=0;i<n;i++) {
        members[i]=allocateNumList(10);
        members[i].data[0]=i;
        members[i].n=1;
    }

    // inverse of covariance matrix of members of given cluster
    // (representing the shape of the clusters)
    //R: invcov<-rep(list(fakeInvCov),clusterCount)
    invcov=(double *)MEM_ALLOC(clusterCount*p*p,sizeof(double));
    invcovMerged=(double *)MEM_ALLOC(clusterCount*p*p,sizeof(double));
    memset(invcov,0,clusterCount*p*p*sizeof(double));
    for (i=0;i<n;i++) {
        for (j=0;j<p;j++) {
            invcov[p*p*i+(p+1)*j]=1;
        }
    }
    // invcov holds inverses of covariance matrices, stored one after another
    // e.g. like ic1_11,ic1_21,ic1_12,ic1_22, ic2_11,ic2_21,ic2_12,ic2_22.
    // Normalizing factor for each cluster (computed from determinant of
    // the inverse of the covariance matrix) making the N-dim volume of
    // clusters equal to 1 if `invcov[[i]]' gets divided by `detsSqrt[i]'.
    detsSqrt=(double*)MEM_ALLOC(clusterCount,sizeof(double));
    for (i=0;i<clusterCount;i++) detsSqrt[i]=1.0;

    // centroids of clusters
    //R: centroid<-x
    centroid=(double *)MEM_ALLOC(n*p,sizeof(double));
    memcpy(centroid,x,sizeof(*x)*n*p);
    DBG_CODE(3,printDoubleMatrix("centroid",centroid,n,p));

    xc1=(double *)MEM_ALLOC(n*p,sizeof(double));
    xc2=(double *)MEM_ALLOC(n*p,sizeof(double));

    // proportional size of each cluster
    weightFactor=(double *)MEM_ALLOC(clusterCount,sizeof(double));
    // number of clusters whose relative size is at least
    // mahalanobis.distance.threshold

    fullMahalClusterCount=0;
    // have all clusters reached the thresh relative size and have we,
    // therefore, switched into "full Mahalanobis" mode?
    switchedToFullMahal=0;

    ///////////////////////////////////////////////////////////////////
    // the main loop
    ///////////////////////////////////////////////////////////////////
    // merge two closest clusters at each step `s'
    for (s=0;s<n-1;s++) {
        double v;
        Num source,target,count;
        Num distXIndicesLen;

        DBG_CODE(3, {
            DBGU("\n====================== step %d ============================\n",s);
            for (i=0;i<clusterCount;i++) {
                DBGU("members[%d]: %s\n",C2R(i),printMembers(members[i]));
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
            //printNumMatrix("clusterSize",clusterSize,1,2*n-1);
            for (i=0;i<2*n-1;i++) {
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
        DBG(4,"found minimum %g at %d\n",v,k);

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
        c1=floor(clusterCount+.5-sqrt(clusterCount*clusterCount-clusterCount+.25-2*(k-1)));
        c2=(c1-1)*(c1*.5-clusterCount)+c1+k;
        // make c1 and c2 0-based
        c1=R2C(c1);
        c2=R2C(c2);
        DBG(2,"c1=%d, c2=%d\n",C2R(c1),C2R(c2));

        //R: if (dbg>0) cat(sprintf('Cluster %d: depth %g, merged clusters %d and %d (%s and %s).\n',s+n,v,clusterId[c1],clusterId[c2],printMembers(members[[c1]]),printMembers(members[[c2]])))
        DBG(1,"Cluster %d: depth %g, merged clusters %d and %d ((%s) and (%s)).\n",C2R(s+n),v,C2R(clusterId[c1]),C2R(clusterId[c2]),printMembers(members[c1]),printMembers(members[c2]));
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
        DBG(3,"clusterCount %d\n",clusterCount);
        for (i=0;i<clusterCount-2;i++) otherClustersTmp[i]=C2R(otherClusters[i]);
        DBG_CODE(3,printNumMatrix("otherClusters",otherClustersTmp,1,clusterCount-2));

        // get all samples constituting the merged clusters
        //R: xij<-x[c(members[[c1]],members[[c2]]),,drop=FALSE]
        c1n=members[c1].n;
        c2n=members[c2].n;
        xijN=c1n+c2n;
        c1indices=members[c1].data;
        c2indices=members[c2].data;
        DBG(2,"c1n=%d, c2n=%d\n",c1n,c2n);
        DBG(3,"constructing xij\n");
        //R: xij<-x[c(members[[c1]],members[[c2]]),,drop=FALSE]
        DBG_CODE(3,printNumMatrix("c1indices",c1indices,1,c1n));
        for (i=0;i<p;i++) {
            int offset1=xijN*i;
            int offset2=n*i;
            for (j=0;j<c1n;j++) {
                xij[j+offset1]=x[c1indices[j]+offset2];
            }
        }
        DBG_CODE(3) printNumMatrix("c2indices",c2indices,1,c2n);
        for (i=0;i<p;i++) {
            int offset1=xijN*i;
            int offset2=n*i;
            for (j=0;j<c2n;j++) {
                xij[c1n+j+offset1]=x[c2indices[j]+offset2];
            }
        }
        DBG_CODE(2) printDoubleMatrix("xij",xij,xijN,p);

        // compute the weight factor controlling the Mahalanobis-Euclidean
        // balance (to be applied when measuring distances relatively to
        // the merged cluster)
        if (thresh > 0) {
            DBG_CODE(2,{
                DBGU("clusterSize[clusterId[c1]] %d\n",clusterSize[clusterId[c1]]);
                DBGU("clusterSize[clusterId[c2]] %d\n",clusterSize[clusterId[c2]]);
                DBGU("(clusterSize[clusterId[c1]] + clusterSize[clusterId[c2]]) / ( n * thresh ) %.3f\n",
                     (clusterSize[clusterId[c1]] + clusterSize[clusterId[c2]]) / ( n * thresh ));
                DBGU("n %d\n",n);
                DBGU("thresh %.3f\n",thresh );
            });
            wf1=fmin2(1.0, (clusterSize[clusterId[c1]] + clusterSize[clusterId[c2]]) / ( n * thresh ) );
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
        DBG_CODE(4,printDoubleMatrix("L",covXij,p,p));
        ASSERT(info>=0,"invalid dpotrf retcode");
        if (info==0) {
            double *cholDecomp=covXij;
            //R: detSqrt<-(1/prod(diag(c.cholDecomp)))^(2/p)
            detSqrt=1.0;
            for (i=0;i<p;i++) {
                detSqrt*=cholDecomp[(p+1)*i];
            }
            detSqrt=pow(detSqrt,-2.0/p);
            //R: invcov_merged<-chol2inv(c.cholDecomp)
            F77_CALL(dpotri)("L",// Lower diagonal is used
                             &p,covXij,&p,&info);
            DBG(4,"info: %d\n",info);
            DBG_CODE(4,printDoubleMatrix("I",covXij,p,p));
            memcpy(invcovMerged,covXij,sizeof(*covXij)*p*p);
            // invcovMerged contains only lower triangular part(plus diagonal),
            // so initialize also the upper diagonal for debugging purposes
            if (dbg) {
                for (i=0;i<p;i++) {
                    for (j=i+1;j<p;j++) {
                        invcovMerged[i+p*j]=invcovMerged[j+p*i];
                    }
                }
            }
        } else {
            detSqrt=1.0;
            memcpy(invcovMerged,fakeInvCov,sizeof(*fakeInvCov)*p*p);
        }

        DBG_CODE(3,{
            printDoubleMatrix("detSqrt",&detSqrt,1,1);
            printDoubleMatrix("invcovMerged",invcovMerged,p,p);
        });
        // compute a new center of the merged cluster
        //R: centroid[c1,]<-(clusterSize[clusterId[c1]]*centroid[c1,,drop=FALSE] + clusterSize[clusterId[c2]]*centroid[c2,,drop=FALSE])/
        //R:      (clusterSize[clusterId[c1]] + clusterSize[clusterId[c2]]) */
        DBG_CODE(2,{
            sprintf(strBuffer,"centroid[c1=%d]",C2R(c1));
            printDoubleMatrixRow(strBuffer,centroid,clusterCount,p,c1);
            sprintf(strBuffer,"centroid[c2=%d]",C2R(c2));
            printDoubleMatrixRow(strBuffer,centroid,clusterCount,p,c2);
        });
        clusterSizeI=clusterSize[clusterId[c1]];
        clusterSizeJ=clusterSize[clusterId[c2]];
        for (i=0;i<p;i++) {
            long offset=clusterCount*i;
            centroid[c1+offset]=(clusterSizeI*centroid[c1+offset] + clusterSizeJ*centroid[c2+offset])/
                (clusterSizeI+clusterSizeJ);
        }
        DBG_CODE(2,{
            sprintf(strBuffer,"updated centroid[c1=%d]",C2R(c1));
            printDoubleMatrixRow(strBuffer,centroid,clusterCount,p,c1);
        });
        
        ///////////////////////////////////////////////////////////////
        // update distX
        ///////////////////////////////////////////////////////////////
        if (!(fullMahalClusterCount == clusterCount-1 && !switchedToFullMahal)) {
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

                DBG(2,"oi %d\n",C2R(oi));
                DBG(2,"otherCluster %d\n",C2R(otherCluster));

                for (ii=0;ii<clusterCount;ii++) {
                    DBG(3,"MEMBERS[%d]: ",ii);
                    DBG(3,"%s\n",printMembers(members[ii]));
                }
                // compute the distance from the newly merged cluster c1+c2 to cluster otherClusters(oi)

                if (quick) {
                    //R: xc1<-centroid[otherClusters[oi],,drop=FALSE]
                    xc1MemberCount=1;
                    for (i=0;i<p;i++) {
                        xc1[i]=centroid[otherCluster+clusterCount*i];
                    }
                } else {
                    //R: xc1<-x[members[[otherClusters[oi]]],,drop=FALSE]
                    xc1MemberCount=members[otherCluster].n;
                    for (i=0;i<p;i++) {
                        Num offset1=xc1MemberCount*i;
                        Num offset2=n*i;
                        Num cumI=0;
                        Num *mmbrs=members[otherCluster].data;
                        for (j=0;j<xc1MemberCount;j++,cumI++) {
                            xc1[cumI+offset1]=x[mmbrs[j]+offset2];
                        }
                    }
                }
                DBG_CODE(3,printDoubleMatrix("xc1",xc1,xc1MemberCount,p));
                DBG_CODE(3,{
                    sprintf(strBuffer,"centroid[c1=%d]",C2R(c1));
                    printDoubleMatrixRow(strBuffer,centroid,clusterCount,p,c1);
                });
                //R: xc1<-xc1-matrix(centroid[c1,,drop=FALSE],nrow(xc1),ncol(xc1),byrow=TRUE)
                for (i=0;i<p;i++) {
                    long offset1=xc1MemberCount*i;
                    long offset2=clusterCount*i;
                    for (j=0;j<xc1MemberCount;j++) {
                        xc1[j+offset1]-=centroid[c1+offset2];
                    }
                }
                DBG_CODE(3,printDoubleMatrix("centered xc1",xc1,xc1MemberCount,p));
                // distMaha1 holds a vector of squares of mahalanobis distances:
                // mean Mahalanobis distance from c1+c2 to some other cluster
                //R: distMaha1<-mean(sqrt(rowSums((xc1%*%ic1)*xc1)))
                distMaha1=computeMahalanobisDistance(xc1,xc1MemberCount,p,ic1,dbg);
                DBG(2,"distMaha1 %g\n",distMaha1);

                // compute the distance from cluster otherClusters(oi) to the newly merged cluster c1+c2
                if (quick) {
                    //R: xc2<-centroid[c1,,drop=FALSE]
                    xc2MemberCount=1;
                    for (i=0;i<p;i++) {
                        xc2[i]=centroid[c1+clusterCount*i];
                    }
                } else {
                    xc2MemberCount=xijN;
                    memcpy(xc2,xij,sizeof(*xij)*xijN*p);
                }
                DBG_CODE(3,printDoubleMatrix("xc2",xc2,xc2MemberCount,p));
                DBG_CODE(3,{
                    sprintf(strBuffer,"centroid[oi=%d]",C2R(otherCluster));
                    printDoubleMatrixRow(strBuffer,centroid,clusterCount,p,otherCluster);
                });

                //R: xc2<-xc2-matrix(centroid[otherClusters[oi],,drop=FALSE],nrow(xc2),ncol(xc2),byrow=TRUE)
                for (i=0;i<p;i++) {
                    long offset1=xc2MemberCount*i;
                    long offset2=clusterCount*i;
                    for (j=0;j<xc2MemberCount;j++) {
                        xc2[j+offset1]-=centroid[otherCluster+offset2];
                    }
                }
                DBG_CODE(3,printDoubleMatrix("centered xc2",xc2,xc2MemberCount,p));

                //R: ic2<-invcov[[otherClusters[ii]]]
                for (i=0;i<p*p;i++) {
                    ic2[i]=invcov[p*p*otherCluster+i];
                }

                DBG_CODE(2,printDoubleMatrix("ic2",ic2,p,p));

                if (normalize || fullMahalClusterCount < clusterCount-1) {
                    DBG(3," normalizing (normalize %d, clusters with full Mahalanobis = %d, clusters =  %d)\n",
                        normalize,fullMahalClusterCount,clusterCount);
                    for (i=0;i<p*p;i++) {
                        ic2[i]/=detsSqrt[otherCluster];
                    }
                    DBG_CODE(2,printDoubleMatrix("normalized ic2",ic2,p,p));
                }

                // mean Mahalanobis distance from otherCluster to c1+c2
                //R: distMaha2<-mean(sqrt(rowSums((xc2%*%ic2)*xc2)))
                distMaha2=computeMahalanobisDistance(xc2,xc2MemberCount,p,ic2,dbg);
                DBG(2,"distMaha2 %g\n",distMaha2);

                // merge the clusterId(c1) <-> clusterId(c2) distances
                //R: iRelDistXIdx<-c(otherClusters1*(clusterCount-(otherClusters1+1)/2)-clusterCount+c1,
                //R:     c1*(clusterCount-(c1+1)/2)-clusterCount+otherClusters2,
                //R:     c1*(clusterCount-(c1+1)/2)-clusterCount+otherClusters3)
                //R: distX[iRelDistXIdx[ii]]<-mean(c(distMaha1,distMaha2))
                if (oi<c1) {
                    distXIdx=R2C(C2R(otherCluster-1)*clusterCount-(C2R(otherCluster)*(C2R(otherCluster)+1))/2+C2R(c1));
                } else if (oi<c2-1) {
                    distXIdx=R2C(C2R(c1-1)*clusterCount-(C2R(c1)*(C2R(c1)+1))/2+C2R(otherCluster));
                } else {
                    distXIdx=R2C(C2R(c1-1)*clusterCount-(C2R(c1)*(C2R(c1)+1))/2+C2R(otherCluster));
                }
                DBG(4,"distXIdx %d\n",distXIdx);
                tmp=distX[distXIdx]=maxDistX-(distMaha1+distMaha2)/2;
                if (tmp<minNewDistX) minNewDistX=tmp;

                DBG(3,"otherCluster %d\n",C2R(otherCluster));
                //DBG("Dist from %d=(%s)\n",C2R(clusterId[otherCluster]),printMembers(members[otherCluster]));
                DBG(2,"Dist from %d=(%s) to %d=(%s %s): %g.\n",
                    C2R(clusterId[otherCluster]),printMembers(members[otherCluster]),
                    C2R(s+n),printMembers(members[c1]),printMembers(members[c2]),
                    maxDistX-distX[distXIdx]);
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

         // clusters clusterId[c1] and clusterId[c2] merged, remove
         // clusterId[c2]-related info, put info about the newly created
         // cluster at position occupied by clusterId[c1] previously
        //DBG(2,"  updating...\n");
        //R: members[[c1]]<-c(members[[c1]],members[[c2]])
        n1=members[c1].n;
        n2=members[c2].n;
        cap1=members[c1].capacity;
        cap2=members[c2].capacity;
        DBG(2,"  members: n1=%d, n2=%d\n",n1,n2);
        if (cap1>=n1+n2 || cap2>=n1+n2) {
            // reuse of one of members[c1] members[c2]
            if (cap1>=cap2) {
                // copy members[c2] at the end of members[c1]
                for (i=0;i<n2;i++) {
                    members[c1].data[n1+i]=members[c2].data[i];
                }
                members[c1].n+=n2;
                MEM_FREE(members[c2].data);
                members[c2].data=NULL;
            } else {
                // copy members[c1] at the end of members[c2]
                for (i=0;i<n1;i++) {
                    members[c2].data[n2+i]=members[c1].data[i];
                }
                members[c2].n+=n1;
                MEM_FREE(members[c1].data);
                members[c1]=members[c2];
                members[c2].data=NULL;
            }
        } else {
            // allocate a new buffer for members[c1]
            Num *newData=(Num *)MEM_ALLOC(2*(n1+n2),sizeof(Num));
            for (i=0;i<n1;i++) {
                newData[i]=members[c1].data[i];
            }
            for (i=0;i<n2;i++) {
                newData[n1+i]=members[c2].data[i];
            }
            MEM_FREE(members[c1].data);
            MEM_FREE(members[c2].data);
            members[c1].n=n1+n2;
            members[c1].data=newData;
            members[c2].data=NULL;
        }
        //R: members<-members[-c2]
        //R: for (i=c2;i<clusterCount-1;i++) members[i]=members[i+1];
        memmove(members+c2,members+c2+1,sizeof(*members)*(clusterCount-1-c2));
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
        for (i=0;i<c1;i++) distXIndices[i]=R2C(C2R(i)*(2*clusterCount-(C2R(i)+1))/2-clusterCount+C2R(c2));
        distXIndices[c1]=R2C(C2R(c1)*(2*clusterCount-(C2R(c1)+1))/2-clusterCount+C2R(c2));
        for (i=c1+1;i<c2;i++) distXIndices[i]=R2C(C2R(i)*(2*clusterCount-(C2R(i)+1))/2-clusterCount+C2R(c2));
        for (i=c2+1;i<clusterCount;i++) distXIndices[i-1]=R2C(C2R(c2)*(2*clusterCount-(C2R(c2)+1))/2-clusterCount+C2R(i));
        distXIndicesLen=clusterCount-1;
        DBG_CODE(4,{
            DBGU("c1=%d, c2=%d, clusterCount=%d\n",c1,c2,clusterCount);
            printNumMatrix("distXIndices part 1 (... c1)",distXIndices,1,c1);
            printNumMatrix("distXIndices part 2 (c1)",distXIndices+c1,1,1);
            printNumMatrix("distXIndices part 3 (c1 ... c2)",distXIndices+c1+1,1,c2-c1-1);
            printNumMatrix("distXIndices part 4 (c2 ...)",distXIndices+c2,1,clusterCount-c2-1);
            printNumMatrix("distXIndices",distXIndices,1,distXIndicesLen);
            DBG_CODE(4,printDistMatrix("distX pre update",distX,clusterCount,maxDistX));
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
            } else {
                count=distXLen-distXIndices[ii]-1;
            }
            if (ii+1<distXIndicesLen) {
                DBG(4," distXIndices[i=%d]=%d, distXIndices[ii=%d]=%d, distXIndices[ii+1=%d]=%d\n",
                    i,distXIndices[i],ii,distXIndices[ii],ii+1,distXIndices[ii+1]);
            } else {
                DBG(4," distXIndices[i=%d]=%d, distXIndices[ii=%d]=%d\n",i,distXIndices[i],ii,distXIndices[ii]);
            }
            DBG_CODE(5,{
                for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
                printDoubleMatrix("distX pre move",distXTmp,1,(clusterCount-1)*(clusterCount-2)/2);
            });
            DBG(4,"memmove target %d (%p, value %g), source %d (%p, value %g), len %d (%d)\n",
                target,distX+target,maxDistX-distX[target],source,distX+source,maxDistX-distX[source],count,sizeof(*distX)*count);
            memmove(distX+target,distX+source,sizeof(*distX)*count);
            DBG(5,"target %d (%p, value %g)\n",target,distX+target,maxDistX-distX[target]);
            DBG_CODE(5,{
                for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
                printDoubleMatrix("distX post move",distXTmp,1,(clusterCount-1)*(clusterCount-2)/2);
            });
            target+=count;
            i=ii+1;
        }
        distXLen-=clusterCount-1;
        DBG_CODE(3) {
            DBGU("clusterCount %d\n",clusterCount);
            DBGU("distXLen %d\n",distXLen);
            for (i=0;i<distXLen;i++) distXTmp[i]=maxDistX-distX[i];
            printDoubleMatrix("distX",distXTmp,1,(clusterCount-1)*(clusterCount-2)/2);
        }
        DBG_CODE(3,printDistMatrix("distX",distX,clusterCount-1,maxDistX));

        // update clusterCount, clusterSize, clusterId
        clusterCount--;
        clusterSize[n+s]=clusterSize[clusterId[c1]]+clusterSize[clusterId[c2]];
        clusterId[c1]=n+s;
        //R: clusterId<-clusterId[-c2]
        memmove(clusterId+c2,clusterId+c2+1,sizeof(*clusterId)*(clusterCount-c2));
        DBG_CODE(4,printNumMatrix("clusterId",clusterId,1,clusterCount));

        if (fullMahalClusterCount == clusterCount && !switchedToFullMahal) {
            int idx,i1,i2;
            double tmp,minNewDistX=0;

            switchedToFullMahal=1;

            DBG(2,"Recomputing all distances.\n");
            idx=0;

            for (i1=0;i1<clusterCount-1;i1++) {

                Num xc1n=constructRepresentantsOfCluster(x,n,p,centroid,members,
                                                         i1,xc1,quick,dbg);
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
                    distMaha1=computeMahalanobisDistance(y,xc1n,p,ic1,dbg);
                    DBG(2,"distMaha1 %g\n",distMaha1);

                    // compute the distance from cluster (i2) to (i1)
                    yn=constructBetweenClusterDistanceMatrix(x,n,p,centroid,clusterCount,
                                                             members,i2,i1,y,
                                                             quick,dbg);
                    DBG_CODE(3,printDoubleMatrix("xc2",y,yn,p));
                    distMaha2=computeMahalanobisDistance(y,yn,p,ic2,dbg);
                    DBG(2,"distMaha2 %g\n",distMaha2);

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
    MEM_FREE(covMeansTmp);
    MEM_FREE(clusterSize);
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
    for (i=0;i<n;i++) {
        deallocateNumList(members[i]);
    }
    MEM_FREE(members);
    MEM_FREE(invcov);
    MEM_FREE(invcovMerged);
    MEM_FREE(detsSqrt);
    MEM_FREE(centroid);
    MEM_FREE(xc1);
    MEM_FREE(xc2);
    MEM_FREE(weightFactor);

	UNPROTECT(nprot);

	DBG(2,"mhclust_ finishes\n");

	//return(NEW_LIST(0));
    return(R_NilValue);
}
