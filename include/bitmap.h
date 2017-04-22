#ifndef __BITMAP_H_INCLUDED__
#define __BITMAP_H_INCLUDED__

#define FILEHEADERSIZE 14
#define INFOHEADERSIZE 40
#define HEADERSIZE (FILEHEADERSIZE+INFOHEADERSIZE)

#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include<math.h>

#include <omp.h>

//#define BMP_USE_ACC

//#define BMP_USE_GPU

#include <amp.h>
#include <amp_math.h>
using namespace concurrency;
#ifdef BMP_USE_GPU
#define BMP_ACC_RESTRICTION restrict(amp)
#define ACC_TYPE "GPU"
#define FAST_MATH	fast_math::
#else if
#define BMP_ACC_RESTRICTION restrict(cpu)
#define ACC_TYPE "CPU"
#define FAST_MATH
#endif


//#define OMP_SCHEDULE_BMP schedule(guided)
#define OMP_SCHEDULE_BMP schedule(static)
#define LINELENGMAX	(4096*10)

typedef struct Rgb_{
	float b;
	float g;
	float r;
	~ Rgb_(){}
	Rgb_(){}
	inline Rgb_(int x, int y, int z)
	{
		r = x;
		g = y;
		b = z;
	}
	inline Rgb_(const int* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
	}
	inline Rgb_(const unsigned char* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
	}
	inline Rgb_(const unsigned char x)
	{
		r = x;
		g = x;
		b = x;
	}
} Rgb;

typedef struct iRgb_{
	float b;
	float g;
	float r;
	~ iRgb_(){}
	iRgb_(){}
	inline iRgb_(int x, int y, int z)
	{
		r = x;
		g = y;
		b = z;
	}
	inline iRgb_(const int* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
	}
	inline iRgb_(const unsigned int* x)
	{
		r = x[0];
		g = x[1];
		b = x[2];
	}
} iRgb;


typedef struct{
	unsigned int height;
	unsigned int width;
	Rgb *data;
}Image;

class ImageArray
{
public:
	unsigned int height;
	unsigned int width;
	double* data;

	ImageArray(Image* img)
	{
		height = img->height;
		width = img->width;
		data = new double[3 * height*width];

		int i = 0;
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				data[i] = img->data[y*width + x].r;
				i++;
			}
		}
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				data[i] = img->data[y*width + x].g;
				i++;
			}
		}
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				data[i] = img->data[y*width + x].b;
				i++;
			}
		}
	}
	~ImageArray()
	{
		delete[] data;
	}

	ImageArray(ImageArray* src)
	{
		height = src->height;
		width = src->width;
		data = new double[3 * height*width];

		memcpy(data, src->data, 3 * height*width*sizeof(double));
	}

	void Copy(ImageArray* src)
	{
		height = src->height;
		width = src->width;
		memcpy(data, src->data, 3 * height*width*sizeof(double));
	}

	void ToImage(Image* img)
	{
		int i = 0;
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				img->data[y*width + x].r = data[i];
				i++;
			}
		}
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				img->data[y*width + x].g = data[i];
				i++;
			}
		}
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				img->data[y*width + x].b = data[i];
				i++;
			}
		}
	}
};


class BitMap
{
	Image* data;
	Image *Read_Bmp(char *filename);
	void Free_Image(Image *img);
	int Write_Bmp(char *filename, Image *img);
	Image *Create_Image(int width, int height);
	int WriteText(char* filename, Image* img);
	Image* Read_Text(char *filename);
	Image* Read_Csv(char* filename, double min, double max);
	int Write_Csv(char* filename, Image* img, int rgb);
	int Write_Csv(char* filename, Image* img, int rgb, double min, double max);

public:

	inline BitMap()
	{
		data = NULL;
	}

	inline ~BitMap()
	{
		Clear();
	}

	inline void Clear()
	{
		if ( data ) Free_Image(data);
		data = NULL;
	}

	void Create(int width, int height)
	{
		data = Create_Image(width, height);
	}

	inline Image* GetImage() const
	{ 
		return data;
	}
	inline int W() const
	{
		return data->width;
	}
	inline int H() const
	{
		return data->height;
	}
	inline void Copy(const BitMap& bmp)
	{
		Clear();
		Create(bmp.data->width, bmp.data->height);
		memcpy(data->data, bmp.data->data, sizeof(Rgb)*bmp.data->width*bmp.data->height);
	}

	void Write(char *filename)
	{
		Write_Bmp(filename, data);
	}
	void Read(char *filename)
	{
		data = Read_Bmp( filename);
	}

	void WriteText(char *filename)
	{
		WriteText( filename, data);
	}

	void ReadText(char *filename)
	{
		Read_Text(filename);
	}

	inline Rgb& cell(const int i, const int j) const
	{
		//return *(data->data+((data->height-j-1)*data->width + i));
		return *(data->data + ((j)*data->width + i));
	}
	inline Rgb& cell_(const int i, const int j) const
	{
		//return *(data->data + ((data->height - i - 1)*data->width + j));
		return *(data->data + (i*data->width + j));
	}
	inline Rgb& cell(const int i) const
	{
		return *(data->data + i);
	}


	void ReadCsv(char* filename, double min, double max)
	{
		if ( data ) Free_Image(data);
		data = Read_Csv(filename, min, max);
	}
	void ReadCsv(double* value, int w, int h, double min, double max)
	{
		if ( data ) Free_Image(data);
		data = Read_Csv(value, w, h, min, max);
	}

	void WriteCsv(char* filename, int rgb)
	{
		Write_Csv( filename, data, rgb);
	}
	void WriteCsv(char* filename, int rgb, double min, double max)
	{
		Write_Csv( filename, data, rgb, min, max);
	}

	void Reverse();

	void Offset(int size);
	void convolve_smooth(int* mask, double mat[3][3]);
	
	void convolve_smooth(int* mask)
	{
		double tmp[3],conv[3][3]={{1.0, 1.0, 1.0},
		{1.0, 1.0, 1.0},
		{1.0, 1.0, 1.0}};
		convolve_smooth(mask, conv);
	}
	
	Image* Read_Csv(double* value, int w, int h, double min, double max);


	static int colortableNum;
	static unsigned char colorTbl[1024][3];
	void ColorTable();
	void ColorTable(BitMap& colormap);
	void ColorTable(int startIndex, int endIndex,unsigned char start[3], unsigned char end[3]);
	void ColorLevel( double min, double max, double* z, double zmask, unsigned char* maskcolor=NULL, int* top=NULL, double* elv=NULL);

	inline void decompose3(BitMap &r, BitMap &g, BitMap &b)
	{

		r.Create(W(), H());
		g.Create(W(), H());
		b.Create(W(), H());

		for (int y = 0; y < H(); y++)
		{
			for (int x = 0; x <W(); x++)
			{
				r.cell(x, y).r = cell(x, y).r;
				g.cell(x, y).g = cell(x, y).g;
				b.cell(x, y).b = cell(x, y).b;
			}
		}
	}
	void compose3(const BitMap &r, const BitMap &g, const BitMap &b)
	{
		for (int y = 0; y < H(); y++)
		{
#pragma omp parallel for
			for (int x = 0; x <W(); x++)
			{
				float pr = r.cell(x, y).r;
				float pg = g.cell(x, y).g;
				float pb = b.cell(x, y).b;

				// クリッピング
				pr = min(max(pr, 0.f), 255.f);
				pg = min(max(pg, 0.f), 255.f);
				pb = min(max(pb, 0.f), 255.f);
				cell(x, y) = Rgb(pr, pg, pb);
			}
		}
	}

	void resize(const float scalex, const float scaley, BitMap& dst)
	{
		int w = W()*scalex;
		int h = H()*scaley;
		dst.Create(w, h);

		printf("%d %d\n", dst.W(), dst.H());
		fflush(stdout);

		for (int y = 0; y < dst.H(); y++)
		{
#pragma omp parallel for
			for (int x = 0; x <dst.W(); x++)
			{
				dst.cell(x, y) = Rgb(cell((int)floor(x / scalex), (int)floor(y / scaley)));
			}
		}
		printf("%d %d\n", dst.W(), dst.H());
		fflush(stdout);
	}
};


#endif /*__BITMAP_H_INCLUDED__*/
