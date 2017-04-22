#include "sparse.h"
#include "Colorization_Using_Optimization_solver.h"

#include <time.h>

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		fprintf(stderr, "%s GrayImageFile ImageMarkedFile [-wd windsize -m ImageMaskFile -o output.bmp -d]\n", argv[0]);
		fprintf(stderr, "-d -> debug\n");
		return -1;
	}

	clock_t start = clock();
	BitMap gI;
	BitMap g0;
	BitMap cI;
	BitMap c0;
	BitMap gI_cI;
	BitMap abs_gI_cI;
	BitMap sum_abs_gI_cI;
	BitMap colorIm;
	BitMap colorIm2;
	BitMap sgI;
	BitMap scI;
	BitMap ntscIm;

	int windsize = 1;
	double eps = 0.01;

	int debug = 0;
	int test = 0;
#ifdef _used_Anat_Levin_VERSION
	test = 1;
#endif

	char out[512];
	out[0] = '\0';

	char maskFile[256];
	maskFile[0] = '\0';

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-m") == 0)
		{
			strcpy(maskFile, argv[i + 1]);
		}
		if (strcmp(argv[i], "-o") == 0)
		{
			strcpy(out, argv[i + 1]);
		}
		if (strcmp(argv[i], "-wd") == 0)
		{
			windsize = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-test") == 0)
		{
			test = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-d") == 0)
		{
			debug = 1;
		}
	}

	char drive[256];
	char dir[256];
	char name[256];
	char ext[256];
	_splitpath(argv[1], drive, dir, name, ext);


	if ( out[0] == '\0')
	{
		sprintf(out, "%s%s%s_Out%s", drive, dir, name, ".bmp");
	}


	g0.Read(argv[1]);
	if (debug ) g0.Write("debug_inputGray.bmp");

	if (g0.W() < 10 || g0.H() < 10)
	{
		printf("ERROR Image size is too small.\n10x10 or more is required at a minimum\n");
		return -1;
	}

	if ( maskFile[0] == '\0' )
	{
		c0.Read(argv[2]);

		//BitMap tmp;
		//tmp.Copy(c0);

		////白黒画像とマスク付白黒画像の白黒部分で差が生じていると問題が生じる＝＞差があるとそこもマスク（拘束）とみなされてしまう。
		//for (int y = 0; y < c0.H(); y++)
		//{
		//	for (int x = 0; x < c0.W(); x++)
		//	{
		//		float r = c0.cell(x, y).r - g0.cell(x, y).r;
		//		float g = c0.cell(x, y).g - g0.cell(x, y).g;
		//		float b = c0.cell(x, y).b - g0.cell(x, y).b;
		//		if (fabs(r) < 100 && fabs(g) < 100 && fabs(b) < 100)
		//		{
		//			tmp.cell(x, y).r = g0.cell(x, y).r;
		//			tmp.cell(x, y).g = g0.cell(x, y).g;
		//			tmp.cell(x, y).b = g0.cell(x, y).b;
		//		}
		//	}
		//}
		//tmp.Write("debug_marker2.bmp");
		//exit(0);

	}else
	{
		c0.Read(maskFile);

		for (int y = 0; y < c0.H(); y++)
		{
#pragma omp parallel for
			for (int x = 0; x < c0.W(); x++)
			{
				if ( c0.cell(x,y).r < 1.0 &&  c0.cell(x,y).g < 1.0 &&  c0.cell(x,y).b < 1.0 )
				{
					c0.cell(x, y).r = g0.cell(x, y).r;
					c0.cell(x, y).g = g0.cell(x, y).g;
					c0.cell(x, y).b = g0.cell(x, y).b;
				}
			}
		}
	}

	if (windsize > g0.W() / 10 || windsize > g0.H() / 10)
	{
		windsize = 1;
	}
	if (debug) c0.Write("debug_inputMarked.bmp");
	BitMap tmp;


	//gI=double(g0)/255;
	Reg(g0, gI);

	//cI=double(c0)/255;
	Reg(c0, cI);


	//gI-cI 
	Sub(gI, cI, gI_cI);

	if (debug)
	{
		X255(gI_cI, tmp);
		tmp.Write("debug_grayImg-markedImag.bmp");
	}
	//abs(gI-cI)
	Abs(gI_cI, abs_gI_cI);
	if (debug)
	{
		X255(abs_gI_cI, tmp);
		tmp.Write("debug_abs(grayImg-markedImag).bmp");
	}

	//sum(abs(gI-cI),3)
	SumRGB(abs_gI_cI, sum_abs_gI_cI);
	

	if (debug)
	{
		X255(sum_abs_gI_cI, tmp);
		tmp.Write("debug_Restraint_condition.bmp");
	}
	Bool(sum_abs_gI_cI, eps, colorIm);
	colorIm2.Copy(colorIm);


	if (debug)
	{
		X255(colorIm2, tmp);
		tmp.Write("debug_Restraint_conditionFlag.bmp");

		for (int y = 0; y < g0.H(); y++)
		{
#pragma omp parallel for
			for (int x = 0; x < g0.W(); x++)
			{
				if (colorIm2.cell(x, y).r > 0.99)
				{
					tmp.cell(x, y).r = c0.cell(x, y).r;
					tmp.cell(x, y).g = c0.cell(x, y).g;
					tmp.cell(x, y).b = c0.cell(x, y).b;
				}
				else
				{
					tmp.cell(x, y).r = 0;
					tmp.cell(x, y).g = 0;
					tmp.cell(x, y).b = 0;
				}
			}
		}
		tmp.Write("debug_Restraint_conditionColorFlag.bmp");

	}

	//input image to YIQ map NTSC image.
	YIQ(gI, sgI);
	if (debug)
	{
		X255(sgI, tmp);
		tmp.Write("debug4.bmp");
	}

	YIQ(cI, scI);
	if (debug)
	{
		X255(scI, tmp);
		tmp.Write("debug5.bmp");
	}

	ntscIm.Copy(g0);
	for (int y = 0; y < g0.H(); y++)
	{
#pragma omp parallel for
		for (int x = 0; x < g0.W(); x++)
		{
			ntscIm.cell(x, y).r = sgI.cell(x, y).r;
			ntscIm.cell(x, y).g = scI.cell(x, y).g;
			ntscIm.cell(x, y).b = scI.cell(x, y).b;
		}
	}
	if (debug)
	{
		X255(ntscIm, tmp);
		tmp.Write("debug6.bmp");
	}

#if 0
	int max_d = floor(log(min(ntscIm.H(), ntscIm.W())) / log(2) - 2.0);
	int iu = floor(ntscIm.H() / (pow(2.0, (max_d - 1))))*(pow(2.0,(max_d - 1)));
	int ju = floor(ntscIm.W() / (pow(2.0, (max_d - 1))))*(pow(2.0, (max_d - 1)));

	printf("%d %d %d\n", max_d, iu, ju);

	BitMap colorIm3;
	BitMap ntscIm2;

	colorIm3.Create(ju, iu);
	ntscIm2.Create(ju, iu);

	for (int y = 0; y < iu; y++)
	{
		for (int x = 0; x < ju; x++)
		{
			colorIm3.cell(x, y).r = colorIm2.cell(x, y).r;
			colorIm3.cell(x, y).g = colorIm2.cell(x, y).g;
			colorIm3.cell(x, y).b = colorIm2.cell(x, y).b;

			ntscIm2.cell(x, y).r = ntscIm.cell(x, y).r;
			ntscIm2.cell(x, y).g = ntscIm.cell(x, y).g;
			ntscIm2.cell(x, y).b = ntscIm.cell(x, y).b;
		}
	}
#else
	int iu = pow2(ntscIm.H());
	int ju = pow2(ntscIm.W());

	//printf("%d %d\n", iu, ju);

	BitMap colorIm3;
	BitMap ntscIm2;

	colorIm3.Create(ju, iu);
	ntscIm2.Create(ju, iu);

	for (int y = 0; y < iu; y++)
	{
#pragma omp parallel for
		for (int x = 0; x < ju; x++)
		{
			if (x >= colorIm2.W() || y >= colorIm2.H())
			{
				continue;
			}
			colorIm3.cell(x, y).r = colorIm2.cell(x, y).r;
			colorIm3.cell(x, y).g = colorIm2.cell(x, y).g;
			colorIm3.cell(x, y).b = colorIm2.cell(x, y).b;

			ntscIm2.cell(x, y).r = ntscIm.cell(x, y).r;
			ntscIm2.cell(x, y).g = ntscIm.cell(x, y).g;
			ntscIm2.cell(x, y).b = ntscIm.cell(x, y).b;
		}
	}

#endif

	BitMap output_fin;

#ifdef _used_Anat_Levin_VERSION
	int Multi_Grid_solver = 1;
#else
	int Multi_Grid_solver = 0;
#endif

#ifdef _used_Anat_Levin_VERSION
	if ( Multi_Grid_solver)
	{
		ImageArray colorIm3p(colorIm3.GetImage());
		ImageArray ntscIm2p(ntscIm2.GetImage());
		ImageArray n0(ntscIm2.GetImage());

		//Multi Grid solver
		getvolcolor(&colorIm3p, &ntscIm2p, 5, 1, &n0);


		BitMap output;
		BitMap output2;

		output.Create(ju, iu);

		n0.ToImage(output.GetImage());
		if (debug)
		{
			X255(output, tmp);
			tmp.Write("debug7.bmp");
		}
		YIQ2RGB(output, output2);
		
		if (debug)
		{
			X255(output2, tmp);
			tmp.Write("debug8.bmp");
		}

		output_fin.Create(colorIm2.W(), colorIm2.H());

		for (int y = 0; y < colorIm2.H(); y++)
		{
#pragma omp parallel for
			for (int x = 0; x < colorIm2.W(); x++)
			{
				output_fin.cell(x, y).r = output2.cell(x, y).r;
				output_fin.cell(x, y).g = output2.cell(x, y).g;
				output_fin.cell(x, y).b = output2.cell(x, y).b;
			}
		}

		X255(output_fin, tmp);
		if ( debug ) tmp.Write("debug9.bmp");

		tmp.Write(out);
	}
#endif

	if ( !test )
	{
		getColorExactSolver(colorIm2, ntscIm, windsize, output_fin);
		if (debug) output_fin.Write("debug10.bmp");
		output_fin.Write(out);
	}

	printf("Time %.2f[sec]\n", (double)(clock() - start) / (double)CLOCKS_PER_SEC);
	return 0;
}