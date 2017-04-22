#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#ifndef SCALAR
   #define SCALAR double
#endif

typedef SCALAR *vector, **matrix;

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))
/*
void nrerror(char error_text[]);
*/

int *Alloc_ivector(long nl, long nh);
double *Alloc_dvector(long nl, long nh);
double **Alloc_dmatrix(long nrl, long nrh, long ncl, long nch);
void Free_ivector(int *v, long nl, long nh);
void Free_dvector(double *v, long nl, long nh);
void Free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

#endif /* _NR_UTILS_H_ */
