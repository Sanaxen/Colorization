#ifndef _INC_NR_LINEAR
#define _INC_NR_LINEAR

#define EPS 1.0e-14
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//		2.7 Sparse Linear Systems																	//
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Sparse-matrix ストレージデータ構造
typedef struct sparse_type {
	int n;
	double thresh;
	double *sa;
	unsigned long *ija;
	int k_index_max;
	int nonzero_num;
	int cur_row;	/* 1 <= cur_row <= n */
	int cur_col;	/* 1 <= cur_col <= n */
	int k_index;
} sparse_t;

typedef struct sparse_nzero_data_type{
	int	index[2];
	double data;
	char flg;
}sparse_nzero_data_t;


int sparse_init(int n, int nmax, double thresh, sparse_t* sparse);
int sparse_term(sparse_t* sparse);

int sparse_next_where( sparse_t* sparse, int *row, int *col);
int sparse_set( double a, sparse_t* sparse);
int SLS_dsprsin3(int nonzeronum, sparse_nzero_data_t *a, int n, double thresh, sparse_t* sparse);
int SLS_dsprsin2(double *a, int n, double thresh, sparse_t* sparse);

int SLS_dsprsin(double **a, int n, double thresh, unsigned long nmax, double sa[], unsigned long ija[]);

int SLS_dsprsax(const double sa[], const unsigned long ija[], const double x[], double b[], const unsigned long n);

int SLS_dsprstx(const double sa[], const unsigned long ija[], const double x[], double b[], const unsigned long n);

int SLS_linbcg(const unsigned long n, double b[], double x[], const int itol, const double tol, const int itmax, int *iter, double *err, const double sa[], const unsigned long ija[]);

float SLS_pythag(float a, float b);

double SLS_snrm(const unsigned long n, const double sx[], const int itol);

inline double SLS_snrm1(const unsigned long n, const double sx[])	//SLS_snrm(... itol==1)
{
	double ans = 0.0;
#pragma omp parallel for reduction(+:ans)
	for (long i = 1; i <= n; i++) ans = ans + sx[i] * sx[i];
	return ans;
}


inline int SLS_atimes(const unsigned long n, double x[],  double r[], const int itrnsp, const double sa[], const unsigned long ija[])
{
	//These are double versions of sprsax and sprstx.
	if (itrnsp) return SLS_dsprstx(sa, ija, x, r, n);
	else return SLS_dsprsax(sa, ija, x, r, n);
}

inline int SLS_asolve(const unsigned long n, const double b[], double x[], const int itrnsp, const double sa[], const unsigned long ija[])
{
	//unsigned long i;
#pragma omp parallel for
	for (long i = 1; i <= n; i++) x[i] = (sa[i] != 0.0 ? b[i] / sa[i] : b[i]);
	//The matrix
	//e
	//A is the diagonal part of A, stored in the rst n elements of sa. Since the
	//transpose matrix has the same diagonal, the ag itrnsp is not used.
	return 0;
}

#endif

