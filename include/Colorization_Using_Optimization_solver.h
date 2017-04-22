#ifndef _Colorization_Using_Optimization_solver_H
#define _Colorization_Using_Optimization_solver_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bitmap.h"

void Reg(const BitMap& g, BitMap& c);

void X255(const BitMap& g, BitMap& c);
void Abs(const BitMap& g, BitMap& c);
void Sub(const BitMap& a, BitMap& b, BitMap& c);
void SumRGB(const BitMap& a, BitMap& c);
void Bool(const BitMap& a, const float eps, BitMap& c);
void YIQ(const BitMap& a, BitMap& ntcs);
void YIQ2RGB(const BitMap& a, BitMap& ntcs);

inline double log2(double x)
{
	return log(x) / log(2);
}

inline int pow2(int n)
{
	if (n < 0)
		return 0;
	if (n == 1)
		return 1;

	return (int)pow(2.0, (floor(log2(n - 1)) + 1.0));
}

#ifdef _used_Anat_Levin_VERSION
void getvolcolor(ImageArray *pvar0, ImageArray *pvar1, double pvar6, double pvar7, ImageArray *n0);
#endif

void getColorExactSolver(const BitMap& colorIm, const BitMap& ntscIm_, const int windsize, BitMap& out);

#endif