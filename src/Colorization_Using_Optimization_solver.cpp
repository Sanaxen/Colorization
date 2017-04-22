#include "Colorization_Using_Optimization_solver.h"
#include "SparseMatrix.h"
#include "sparse.h"

void Reg(const BitMap& g, BitMap& c)
{
	c.Copy(g);
	//for (int y = 0; y < g.H(); y++)
	//{
	//	for (int x = 0; x < g.W(); x++)
	//	{
	//		c.cell(x, y).r /= 255.0;
	//		c.cell(x, y).g /= 255.0;
	//		c.cell(x, y).b /= 255.0;
	//	}
	//}
}

void X255(const BitMap& g, BitMap& c)
{
	c.Copy(g);
	//for (int y = 0; y < g.H(); y++)
	//{
	//	for (int x = 0; x < g.W(); x++)
	//	{
	//		c.cell(x, y).r *= 255.0;
	//		c.cell(x, y).g *= 255.0;
	//		c.cell(x, y).b *= 255.0;
	//	}
	//}
}

void Abs(const BitMap& g, BitMap& c)
{
	c.Copy(g);
	for (int y = 0; y < g.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < g.W(); x++)
		{
			c.cell(x, y).r = fabs(c.cell(x, y).r);
			c.cell(x, y).g = fabs(c.cell(x, y).g);
			c.cell(x, y).b = fabs(c.cell(x, y).b);
		}
	}
}

void Sub(const BitMap& a, BitMap& b, BitMap& c)
{

	c.Copy(a);
	for (int y = 0; y < a.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < a.W(); x++)
		{
			c.cell(x, y).r -= b.cell(x, y).r;
			c.cell(x, y).g -= b.cell(x, y).g;
			c.cell(x, y).b -= b.cell(x, y).b;
		}
	}
}

void SumRGB(const BitMap& a, BitMap& c)
{

	c.Copy(a);
	for (int y = 0; y < a.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < a.W(); x++)
		{
			c.cell(x, y).r = a.cell(x, y).r + a.cell(x, y).g + a.cell(x, y).b;
			c.cell(x, y).g = c.cell(x, y).r;
			c.cell(x, y).b = c.cell(x, y).r;
		}
	}
}

void Bool(const BitMap& a, const float eps, BitMap& c)
{

	c.Copy(a);
	for (int y = 0; y < a.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < a.W(); x++)
		{
			c.cell(x, y).r = (a.cell(x, y).r < eps) ? 0 : 255.0;
			c.cell(x, y).g = (a.cell(x, y).g < eps) ? 0 : 255.0;
			c.cell(x, y).b = (a.cell(x, y).b < eps) ? 0 : 255.0;
		}
	}
}

void YIQ(const BitMap& a, BitMap& ntcs)
{
	/*
	Y = 0.299R + 0.587G + 0.114B
	I = 0.596R - 0.274G - 0.322B
	Q = 0.211R - 0.523G + 0.312B
	*/
	ntcs.Copy(a);
	for (int y = 0; y < a.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < a.W(); x++)
		{
			ntcs.cell(x, y).r = a.cell(x, y).r*0.299 + a.cell(x, y).g*0.587 + a.cell(x, y).b*0.114;
			ntcs.cell(x, y).g = a.cell(x, y).r*0.596 - a.cell(x, y).g*0.274 - a.cell(x, y).b*0.322;
			ntcs.cell(x, y).b = a.cell(x, y).r*0.211 - a.cell(x, y).g*0.523 + a.cell(x, y).b*0.312;
		}
	}

}

void YIQ2RGB(const BitMap& a, BitMap& ntcs)
{
	ntcs.Copy(a);
	for (int y = 0; y < a.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < a.W(); x++)
		{
			ntcs.cell(x, y).r = a.cell(x, y).r*1.000 + a.cell(x, y).g*0.956 + a.cell(x, y).b*0.621;
			ntcs.cell(x, y).g = a.cell(x, y).r*1.000 - a.cell(x, y).g*0.272 - a.cell(x, y).b*0.647;
			ntcs.cell(x, y).b = a.cell(x, y).r*1.000 - a.cell(x, y).g*1.106 + a.cell(x, y).b*1.703;
		}
	}

}


void getColorExactSolver(const BitMap& colorIm, const BitMap& ntscIm_, const int windsize, BitMap& out)
{
	int n = ntscIm_.W();	//col
	int m = ntscIm_.H();	//row

	unsigned long int imgSize = n*m;

	BitMap ntscIm;
	ntscIm.Copy(ntscIm_);
	//for (int y = 0; y < imgSize; y++)
	//{
	//	ntscIm.cell(y).r /= 255.0;
	//	ntscIm.cell(y).g /= 255.0;
	//	ntscIm.cell(y).b /= 255.0;
	//}


	// nI(:, : , 1) = ntscIm(:, : , 1);
	BitMap nI;
	nI.Create(n, m);
	for (int y = 0; y < imgSize; y++)
	{
		nI.cell(y).r = ntscIm.cell(y).r;
	}

	const int wd = 1;

	std::vector<int> col_inds;
	std::vector<int> row_inds;
	std::vector<float> vals;

	col_inds.resize(imgSize*pow(2 * wd + 1, 2) + 1);
	row_inds.resize(imgSize*pow(2 * wd + 1, 2) + 1);
	vals.resize(imgSize*pow(2 * wd + 1, 2) + 1);

	int thread_num = 1;
#ifdef _OPENMP
	thread_num = omp_get_max_threads();
#endif
	std::vector<int>* row = new std::vector<int>[thread_num];
	std::vector<int>* col = new std::vector<int>[thread_num];
	std::vector<float>* var = new std::vector<float>[thread_num];

	for (int k = 0; k < thread_num; k++)
	{
		row[k].resize(imgSize);
		col[k].resize(imgSize);
		var[k].resize(imgSize);
	}


	printf("make Matrix elements start\n");
#define LOG_001	-4.60517018598809136804f

	int consts_len = 0;
	for (int i = 0; i < m; i++)
	{
		if (i % 100 == 0) fprintf(stderr, "%d/%d -->%.0f%% \r\r\r\r\r\r\r\r\r\r", i + 1, m, 100.0*(double)(i + 1) / (double)m);
#pragma omp parallel for
		for (int j = 0; j < n; j++)
		{
			if (fabs(colorIm.cell_(i, j).r) < 0.00000001f)
			{
				int thread_id = 0;
#ifdef _OPENMP
				thread_id = omp_get_thread_num();
#endif
				float min_variance = 1.0;
				float sigma;

				int tlen = 0;

				const int iis = max(0, i - wd);
				const int iie = min(i + wd, m - 1);
				const int jjs = max(0, j - wd);
				const int jje = min(j + wd, n - 1);

				float sum_sqr = 0.0, sum = 0.0;
				for (int ii = iis; ii <= iie; ii++)
				{
					for (int jj = jjs; jj <= jje; jj++)
					{
						float val = ntscIm.cell_(ii, jj).r;
						sum += val;
						sum_sqr += val*val;

						if (ii != i || jj != j)
						{
							row[thread_id][tlen] = i*n + j;
							col[thread_id][tlen] = ii*n + jj;
							var[thread_id][tlen] = val - ntscIm.cell_(i, j).r;
							if (0)
							{
								val = ntscIm.cell_(ii, jj).g - ntscIm.cell_(i, j).g;
								var[thread_id][tlen] += val*val;

								val = ntscIm.cell_(ii, jj).b - ntscIm.cell_(i, j).b;
								var[thread_id][tlen] += val*val;
							}
							var[thread_id][tlen] *= var[thread_id][tlen];
							if (var[thread_id][tlen] < min_variance) min_variance = var[thread_id][tlen];
							++tlen;
						}
					}
				}
				sigma = (sum_sqr - (sum*sum) / (float)(tlen + 1)) / (float)tlen;
				if (sigma < 0.000002f)
				{
					sigma = 0.000002f;
				}
				else if (sigma < (-min_variance / LOG_001))
				{
					sigma = -min_variance / LOG_001;
				}

				sum = 0;
				for (int ii = 0; ii < tlen; ii++)
				{
					var[thread_id][ii] = exp(-var[thread_id][ii] / sigma);
					sum += var[thread_id][ii];
				}

				if (sum == 0) continue;

#pragma omp critical
				{
					for (int ii = 0; ii < tlen; ii++)
					{
						row_inds[consts_len] = row[thread_id][ii];
						col_inds[consts_len] = col[thread_id][ii];
						vals[consts_len] = -var[thread_id][ii] / sum;
						++consts_len;
					}
				}
			}
#pragma omp critical
				{
					row_inds[consts_len] = col_inds[consts_len] = i*n + j;
					vals[consts_len] = 1.0;
					consts_len++;
				}
		}
	}
	delete[] row;
	delete[] col;
	delete[] var;

	printf("consts_len %d imgSize %d\n", consts_len, imgSize);
	printf("make Matrix elements end\n");

	double* b = new double[imgSize + 1];
	double* xx = new double[imgSize + 1];
	int	itol;
	int itmax;
	int iter;
	double tol;
	double err;
	itol = 1;
	tol = 1.0e-4;
	itmax = imgSize;
	err = 0.0;

	const float nonzeroValue = 0.00000001f;

	int stat = 0;
	sparse_t sparse;


	//Indexed Storage of Sparse Matrices
	printf("Indexed Storage of Sparse Matrix start.\n"); fflush(stdout);
	if (10){
		printf("\tSparse Matrix elements set start.\n"); fflush(stdout);

		int N = imgSize;
		int nmax = imgSize*(2 * wd + 1)*(2 * wd + 1) + 2;
		stat = sparse_init(N, nmax, nonzeroValue, &sparse);


		std::vector<int> col_min;
		std::vector<int> col_max;
		col_min.resize(N+1);
		col_max.resize(N+1);
		for (int i = 0; i <= N; i++)
		{
			col_min[i] = N+10;
			col_max[i] = -1;
		}
		std::vector<bool> nonezeroV;
		nonezeroV.resize(consts_len);

		SparseMatrix sp(N);
		for (int i = 0; i < consts_len; i++)
		{
			if (i % 100 == 0) fprintf(stderr, "set %d/%d -->%.0f%% \r\r\r\r\r\r\r\r\r\r\r\r", i + 1, consts_len, 100.0*(double)(i + 1) / (double)consts_len);
			const int ii = row_inds[i];
			const int jj = col_inds[i];

			sp.set(i + 1, ii + 1, jj + 1);

			if (col_min[ii + 1] > jj + 1) col_min[ii + 1] = jj + 1;
			if (col_max[ii + 1] < jj + 1) col_max[ii + 1] = jj + 1;
			nonezeroV[i] = false;
			if (fabs(vals[i]) >= nonzeroValue) nonezeroV[i] = true;
		}

		const int check_debug = 1;

		//値が正しくマッピングされているかチェック
		if (check_debug)
		{
#pragma omp parallel for
			for (int i = 0; i < consts_len; i++)
			{
				if (i % 100 == 0) fprintf(stderr, "chk %d/%d -->%.0f%% \r\r\r\r\r\r\r\r\r\r\r\r", i + 1, consts_len, 100.0*(double)(i + 1) / (double)consts_len);
				const int ii = row_inds[i];
				const int jj = col_inds[i];

				int k = sp.get(ii + 1, jj + 1) - 1;
				if (k != i)
				{
					printf("Matrix element error.(%d,%d)\n", ii, jj);
					fflush(stdout);

					exit(0);
				}
			}
		}
		printf("\tSparse Matrix elements set end.\n"); fflush(stdout);
		//printf("OK\n");
		//fflush(stdout);

#pragma omp parallel for
		for (int j = 1; j <= N; j++)
		{
			const int kk = sp.get(j, j);
			if (kk == 0)
			{
				printf("diagonal elements=0!!\n");		fflush(stdout);

				//sparse.sa[j] = 1;
				continue;
			}
			sparse.sa[j] = vals[kk - 1];	//Store diagonal elements.
		}

		sparse.ija[1] = N + 2;				//Index to 1st row o-diagonal element, if any.
		unsigned long k = N + 1;
		for (int i = 1; i <= N; i++) {			//Loop over rows.
			if (i % 100 == 0) fprintf(stderr, "Storage %d/%d -->%.0f%% \r\r\r\r\r\r\r\r\r\r", i, N, 100.0*(double)i / (double)N);
#pragma omp parallel for
			for (int j = col_min[i]; j <= col_max[i]; j++) {		//Loop over columns.

				if (i != j)
				{
					const int kk = sp.get(i, j);
					if (kk <= 0)
					{
						//printf("zero element!!\n");
						//exit(0);
						//zero element!!
						continue;
					}

					if (nonezeroV[kk-1])
					{
						//printf("none zero element\n");
#pragma omp critical
						{
							if (++k > nmax)
							{
								printf("sprsin: nmax too small!!\n");
								exit(0);
							}
							sparse.sa[k] = vals[kk - 1];	//Store o-diagonal elements and their columns.
							sparse.ija[k] = j;
						}
					}
				}
			}
			sparse.ija[i + 1] = k + 1;			//As each row is completed, store index to
			//next.
		}
	}
	printf("Indexed Storage of Sparse Matrix end.\n\n"); fflush(stdout);

	//find all the indices of the nonzeros in "colorlm"
	//memset(b, '\0', sizeof(double)*(imgSize + 1));

	int id = 0;
	for (int j = 0; j < n; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			b[i*n + j + 1] = 0.00001;	//このゴミは結構必要なパラメータ
			if (fabs(colorIm.cell_(i, j).r) > 0.9)
			{
				b[i*n + j + 1] = ntscIm.cell_(i, j).g;
			}
		}
	}

	printf("Conjugate Gradient Method for a Sparse System\n");
	//Sparse Linear Systems
	memset(xx, '\0', sizeof(double)*(imgSize + 1));
	stat = SLS_linbcg(imgSize, b, xx, itol, tol, itmax, &iter, &err, sparse.sa, sparse.ija);
	printf("\nstat %d err %f\n", stat, err);

	for (int j = 0; j < n; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			nI.cell_(i, j).g = xx[i*n + j + 1];
		}
	}


	//find all the indices of the nonzeros in "colorlm"
	//memset(b, '\0', sizeof(double)*(imgSize + 1));
	for (int j = 0; j < n; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			b[i*n + j + 1] = 0.00001;	//このゴミは結構必要なパラメータ
			if (fabs(colorIm.cell_(i, j).r) > 0.9)
			{
				b[i*n + j + 1] = ntscIm.cell_(i, j).b;
			}
		}
	}



	printf("Conjugate Gradient Method for a Sparse System\n");
	//Sparse Linear Systems(Conjugate Gradient Method for a Sparse System)
	memset(xx, '\0', sizeof(double)*(imgSize + 1));
	stat = SLS_linbcg(imgSize, b, xx, itol, tol, itmax, &iter, &err, sparse.sa, sparse.ija);
	printf("\nstat %d err %f\n", stat, err);
	for (int j = 0; j < n; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < m; i++)
		{
			nI.cell_(i, j).b = xx[i*n + j + 1];
		}
	}




	//for (int y = 0; y < imgSize; y++)
	//{
	//	nI.cell(y).r *= 255.0;
	//	nI.cell(y).g *= 255.0;
	//	nI.cell(y).b *= 255.0;
	//	//if (nI.cell(y).r > 0 ) printf("%f %f %f\n", nI.cell(y).r, nI.cell(y).g, nI.cell(y).b);
	//}


	delete[] xx;
	delete[] b;
	row_inds.clear();
	col_inds.clear();
	vals.clear();

	//nI.Write("xxx.bmp");

	BitMap output;
	YIQ2RGB(nI, output);
	//output.Write("xxxx.bmp");

	sparse_term(&sparse);

	out.Copy(output);
	//A = sparse(row_inds, col_inds, vals, consts_len, imgSize);
	//b = zeros(size(A, 1), 1);

	/*
	for t=2:3
	curIm=ntscIm(:,:,t);
	b(lblInds)=curIm(lblInds);
	new_vals=A\b;
	nI(:,:,t)=reshape(new_vals,n,m,1);
	end
	*/


}

