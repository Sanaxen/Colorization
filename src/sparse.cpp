//Linear Algebra functions from Numerical Recipes
//
#ifdef _MSC_VER
#pragma warning(disable: 4244 4305)
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "NR_utility.h"
#include "sparse.h"




// Sparse-matrix ストレージデータ構造の初期化
int sparse_init(int n, int nmax, double thresh, sparse_t* sparse)
{
	int i,j;
	unsigned long k;
	double	*saP;
	unsigned long *ijaP;

	memset(sparse, '\0', sizeof(sparse_t));
	sparse->n = n;
	sparse->thresh = thresh;
	sparse->k_index_max = nmax;

	saP = (double*)malloc( sparse->k_index_max*sizeof(double));
	if ( saP == NULL ){
		return -99;
	}
	memset(saP, '\0', sparse->k_index_max*sizeof(double));

	ijaP = (unsigned long*)malloc( sparse->k_index_max*sizeof(unsigned long));
	if ( ijaP == NULL ){
		return -99;
	}
	memset(ijaP, '\0', sparse->k_index_max*sizeof(unsigned long));

	sparse->cur_col = 1;
	sparse->cur_row = 1;
	sparse->sa = saP;
	sparse->ija = ijaP;
	sparse->ija[1] = n+2;
	sparse->k_index = n+1;
	return 0;
}

int sparse_term(sparse_t* sparse)
{
	if ( sparse->ija != NULL ) free (sparse->ija);
	if ( sparse->sa != NULL ) free (sparse->sa);
	sparse->ija = NULL;
	sparse->sa = NULL;
	return 0;
}




//次にセットすべき配列要素のインデックス要求（０オリジン・インデックス)
int sparse_next_where( sparse_t* sparse, int *row, int *col)
{
	*col = sparse->cur_col-1;
	*row = sparse->cur_row-1;
	if ( sparse->cur_row >= sparse->n && sparse->cur_col >= sparse->n+1){
		return -1;
	}
	return 0;
}

//要求された場所の配列要素の設定
int sparse_set( double a, sparse_t* sparse)
{
	double	*saP, *stmp;
	unsigned long *ijaP, *itmp;


	if ( sparse->k_index+2 >= sparse->k_index_max ){
		stmp = sparse->sa;
		itmp = sparse->ija;

		sparse->k_index_max += 1000;
		saP = (double*)malloc( sparse->k_index_max*sizeof(double));
		if ( saP == NULL ){
			return -99;
		}
		memset(saP, '\0', sparse->k_index_max*sizeof(double));

		ijaP = (unsigned long*)malloc( sparse->k_index_max*sizeof(unsigned long));
		if ( ijaP == NULL ){
			return -99;
		}
		memset(ijaP, '\0', sparse->k_index_max*sizeof(unsigned long));

		memcpy(saP, sparse->sa, (sparse->k_index+1)*sizeof(double));
		memcpy(ijaP, sparse->ija, (sparse->k_index+1)*sizeof(unsigned long));
		sparse->sa = saP;
		sparse->ija = ijaP;
		free( stmp );
		free( itmp );
	}

	// 対角成分
	if ( sparse->cur_col == sparse->cur_row ){
		if ( sparse->cur_row == sparse->n ){
			sparse->ija[sparse->cur_row+1] = sparse->k_index+1;
		}
		sparse->sa[sparse->cur_col] = a;
		sparse->cur_col++;
		return 0;
	}

	//非ゼロ要素の場合
	if ( fabs(a) > sparse->thresh ){
		sparse->nonzero_num++;
		sparse->k_index++;
		sparse->sa[sparse->k_index]=a;
		sparse->ija[sparse->k_index] = sparse->cur_col;
	}
	sparse->cur_col++;
	//列の終端なら次の行の先頭を要求する
	if ( sparse->cur_col == sparse->n+1 ){
		sparse->cur_col = 1;
		sparse->ija[sparse->cur_row+1] = sparse->k_index+1;
		printf("@@@%d\n", sparse->k_index+1);
		sparse->cur_row++;
		return 0;
	}
	return 0;
}

int SLS_dsprsin2(double *a, int n, double thresh, sparse_t* sparse)
{
	int i,j;
	int	nmax;
	unsigned long k;
	int	stat;

	k=n+1;
	for (i=1;i<=n;i++) {			//Loop over rows.
		for (j=1;j<=n;j++) {		//Loop over columns.
			if (fabs(a[n*(i-1)+(j-1)]) >= thresh && i != j) {
				k++;
			}
		}
	}
	nmax = k;
	stat = sparse_init(n, k+2, thresh, sparse);
	if ( stat != 0 ){
		return stat;
	}

	for (j=1;j<=n;j++) sparse->sa[j]=a[n*(j-1)+(j-1)];	//Store diagonal elements.
	sparse->ija[1]=n+2;				//Index to 1st row o-diagonal element, if any.
	k=n+1;
	for (i=1;i<=n;i++) {			//Loop over rows.
		for (j=1;j<=n;j++) {		//Loop over columns.
			if (fabs(a[n*(i-1)+(j-1)]) >= thresh && i != j) {
				if (++k > nmax) return -1;//nrerror("sprsin: nmax too small");
				sparse->sa[k]=a[n*(i-1)+(j-1)];	//Store o-diagonal elements and their columns.
				sparse->ija[k]=j;
				sparse->k_index++;
			}
		}
		sparse->ija[i+1]=k+1;			//As each row is completed, store index to
		//next.
	}
	return 0;
}

int SLS_dsprsin3(int nonzeronum, sparse_nzero_data_t *z, int n, double thresh, sparse_t* sparse)
{
	int i,j, ii;
	int	nmax;
	unsigned long k;
	int	stat;
	int	flg;


	stat = sparse_init(n, nonzeronum+2, thresh, sparse);
	if ( stat != 0 ){
		return stat;
	}

	for (j=1;j<=n;j++){
		for ( ii = 0; ii < nonzeronum; ii++ )
		{
			if (z[ii].index[0] == z[ii].index[1] && j - 1 == z[ii].index[0])
			{
				sparse->sa[j]=z[ii].data;	//Store diagonal elements.
				z[ii].flg = '1';
				break;
			}
		}
	}

	sparse->ija[1]=n+2;				//Index to 1st row o-diagonal element, if any.
	k=n+1;
	for (i=1;i<=n;i++) {			//Loop over rows.
		for (j=1;j<=n;j++) {		//Loop over columns.
			if ( i == j ) continue;
			for ( ii = 0; ii < nonzeronum; ii++ )
			{
				if (z[ii].flg) continue;

				if ( i-1 == z[ii].index[0] && j-1 == z[ii].index[1] ){
					k++;
					sparse->sa[k]=z[ii].data;	//Store diagonal elements.
					sparse->ija[k]=j;
					sparse->k_index++;
					z[ii].flg = '1';
					break;
				}
			}
		}
		sparse->ija[i+1]=k+1;			//As each row is completed, store index to
		//next.
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//		2.7 Sparse Linear Systems																	//
//////////////////////////////////////////////////////////////////////////////////////////////////////
int SLS_dsprsin(double **a, int n, double thresh, unsigned long nmax, double sa[], unsigned long ija[])
//Converts a square matrix a[1..n][1..n] into row-indexed sparse storage mode. Only ele-
//ments of a with magnitude . thresh are retained. Output is in two linear arrays with dimen-
//sion nmax (an input parameter): sa[1..] contains array values, indexed by ija[1..]. The
//number of elements lled of sa and ija on output are both ija[ija[1]-1]-1 (see text).
{
	int i,j;
	unsigned long k;
	for (j=1;j<=n;j++) sa[j]=a[j][j];	//Store diagonal elements.
	ija[1]=n+2;				//Index to 1st row o-diagonal element, if any.
	k=n+1;
	for (i=1;i<=n;i++) {			//Loop over rows.
		for (j=1;j<=n;j++) {		//Loop over columns.
			if (fabs(a[i][j]) >= thresh && i != j) {
				if (++k > nmax) return -1;//nrerror("sprsin: nmax too small");
				sa[k]=a[i][j];	//Store o-diagonal elements and their columns.
				ija[k]=j;
			}
		}
		ija[i+1]=k+1;			//As each row is completed, store index to
		//next.
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int SLS_dsprsax(const double sa[], const unsigned long ija[], const double x[], double b[], const unsigned long n)
//Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x[1..n], giving
//a vector b[1..n].
{
	//unsigned long i,k;
	//if (ija[1] != n+2) return -1;//nrerror("sprsax: mismatched vector and matrix");
#pragma omp parallel for
	for (long i = 1; i <= n; i++) {
		b[i]=sa[i]*x[i];		//Start with diagonal term.
		for (long k=ija[i];k<=ija[i+1]-1;k++)//Loop over o-diagonal terms.
		b[i] += sa[k]*x[ija[k]];
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int SLS_dsprstx(const double sa[], const unsigned long ija[], const double x[], double b[], const unsigned long n)
//Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a vector
//x[1..n], giving a vector b[1..n].
{
	//unsigned long i,j,k;
	//if (ija[1] != n+2) return -1; //nrerror("mismatched vector and matrix in sprstx");
#pragma omp parallel for
	for (long i = 1; i <= n; i++) b[i] = sa[i] * x[i];		//Start with diagonal terms.

	//ここは並列化出来ない。
	for (long i = 1; i <= n; i++) {				//Loop over o-diagonal terms.
		for (long k = ija[i]; k <= ija[i + 1] - 1; k++) {
			long j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
int SLS_linbcg(const unsigned long n, double b[], double x[], const int itol, const double tol, const int itmax, int *iter, double *err, const double sa[], const unsigned long ija[])
//Solves A . x = b for x[1..n], given b[1..n], by the iterative biconjugate gradient method.
//On input x[1..n] should be set to an initial guess of the solution (or all zeros); itol is 1,2,3,
//or 4, specifying which convergence test is applied (see text); itmax is the maximum number
//of allowed iterations; and tol is the desired convergence tolerance. On output, x[1..n] is
//reset to the improved solution, iter is the number of iterations actually taken, and err is the
//estimated error. The matrix A is referenced only through the user-supplied routines atimes,
//which computes the product of either A or its transpose on a vector; and asolve, which solves
//e
//A . x = b or
//e
//A
//T
//. x = b for some preconditioner matrix
//e
//A (possibly the trivial diagonal part of A).
{
	unsigned long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;	//Double precision is a good idea in this routine.
	
	p=Alloc_dvector(1,n);
	pp=Alloc_dvector(1,n);
	r=Alloc_dvector(1,n);
	rr=Alloc_dvector(1,n);
	z=Alloc_dvector(1,n);
	zz=Alloc_dvector(1,n);

#pragma omp parallel for
	for (long j = 0; j < n + 1; j++){
			p[j] = 0.0;
			pp[j] = 0.0;
			r[j] = 0.0;
			rr[j] = 0.0;
			z[j] = 0.0;
			zz[j] = 0.0;
	}
    bkden = 1.0;

	//Calculate initial residual.
	*iter=0;
	SLS_atimes(n,x,r,0, sa, ija);	//Input to atimes is x[1..n], output is r[1..n];
	//the nal 0 indicates that the matrix (not its
	//transpose) is to be used.
#pragma omp parallel for
	for (long j = 1; j <= n; j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	//SLS_atimes(n,r,rr,0);  //Uncomment this line to get the \minimum resid-
	//ual" variant of the algorithm.
	znrm = 1.0;
	if (itol == 1) {
		bnrm=SLS_snrm1(n,b);
		SLS_asolve(n,r,z,0,sa, ija);	//Input to asolve is r[1..n], output is z[1..n];
		//the nal 0 indicates that the matrix
		//e
		//A (not
		//its transpose) is to be used.
	}
	else if (itol == 2) {
		SLS_asolve(n,b,z,0,sa, ija);
		bnrm=SLS_snrm(n,z,itol);
		SLS_asolve(n,r,z,0,sa, ija);
	}
	else if (itol == 3 || itol == 4) {
		SLS_asolve(n,b,z,0,sa, ija);
		bnrm=SLS_snrm(n,z,itol);
		SLS_asolve(n,r,z,0,sa, ija);
		znrm=SLS_snrm(n,z,itol);
	} else return -1; //nrerror("illegal itol in linbcg");

	while (*iter <= itmax) {	//Main loop.
		++(*iter);
		zm1nrm=znrm;
		SLS_asolve(n,rr,zz,1,sa, ija);	//Final 1 indicates use of transpose matrix
		//e
		//A
		//T
		//.
		bknum = 0.0;
#pragma omp parallel for reduction(+:bknum)
		for (long j = 1; j <= n; j++) bknum = bknum + z[j] * rr[j];
		//Calculate coe.cient bk and direction vectors p and pp.
		if (*iter == 1) {
#pragma omp parallel for
			for (long j = 1; j <= n; j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
#pragma omp parallel for
			for (long j = 1; j <= n; j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;		//Calculate coe.cient ak, new iterate x, and new
		//residuals r and rr.
		SLS_atimes(n,p,z,0,sa, ija);

		akden = 0.0;
#pragma omp parallel for reduction(+:akden)
		for (long j = 1; j <= n; j++) akden = akden + z[j] * pp[j];

		ak = 1.0*bknum / akden;
		//if (fabs(ak) > 0.00001) ak = 0.00001;
		SLS_atimes(n,pp,zz,1,sa, ija);
#pragma omp parallel for
		for (long j = 1; j <= n; j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];

			//if (x[j] > 1.0)
			//{
			//	r[j] /= fabs(x[j]);
			//	rr[j] /= fabs(x[j]);
			//	x[j] = 0.99;
			//}
		}
		SLS_asolve(n,r,z,0,sa, ija); //Solve
		//e
		//A . z = r and check stopping criterion.
		if (itol == 1){
			znrm=1.0;
			*err=SLS_snrm1(n,r)/bnrm;
		}
		else if (itol == 2){
			znrm=1.0;
			*err=SLS_snrm(n,z,itol)/bnrm;
		}else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=SLS_snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*SLS_snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;		//Error may not be accurate, so loop again.
				continue;
			}
			xnrm=SLS_snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;		//Error may not be accurate, so loop again.
				continue;
			}
		}
		printf("iter=%4d err=%12.6f\r\r\r\r\r\r\r\r\r\r\r\r\r",*iter,*err);
		if (*err <= tol) break;
	}
	Free_dvector(p,1,n);
	Free_dvector(pp,1,n);
	Free_dvector(r,1,n);
	Free_dvector(rr,1,n);
	Free_dvector(z,1,n);
	Free_dvector(zz,1,n);
	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
double SLS_snrm(const unsigned long n, const double sx[], const int itol)
//Compute one of two norms for a vector sx[1..n], as signaled by itol. Used by linbcg.
{
	//unsigned long i,isamax;
	
	if (itol <= 3) {
		double ans = 0.0;
#pragma omp parallel for reduction(+:ans)
		for (long i = 1; i <= n; i++) ans = ans + sx[i] * sx[i];	//Vector magnitude norm.
		return sqrt(ans);
	} else {
		unsigned long isamax;
		isamax = 1;
		for (long i=1;i<=n;i++) {		//Largest component norm.
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
