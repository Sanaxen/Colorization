#include"bitmap.h"

#define COLORTBL_MAX	1024
unsigned char BitMap::colorTbl[COLORTBL_MAX][3];
int BitMap::colortableNum;

#define STB_IMAGE_IMPLEMENTATION
#include "stb\stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb\stb_image_write.h"

#include <amp.h>
#include <amp_math.h>
using namespace concurrency;

Image *BitMap::Read_Bmp(char *filename)
{
	int i, j;
	int real_width;
	unsigned int width, height;
	unsigned int color;
	FILE *fp;
	unsigned char header_buf[HEADERSIZE];
	unsigned char *bmp_line_data;
	Image *img;

#if 0
	if((fp = fopen(filename, "rb")) == NULL){
		fprintf(stderr, "Error: %s could not read.\n", filename);
		return NULL;
	}

	fread(header_buf, sizeof(unsigned char), HEADERSIZE, fp);
	if(strncmp((char*)header_buf, "BM", 2)){
		fprintf(stderr, "Error: %s is not Bitmap file.\n", filename);
		return NULL;
	}

	memcpy(&width, header_buf + 18, sizeof(width));
	memcpy(&height, header_buf + 22, sizeof(height));
	memcpy(&color, header_buf + 28, sizeof(unsigned int));

	if(color != 24){
		fprintf(stderr, "Error: %s is not 24bit color image\n", filename);
		return NULL;
	}

	real_width = width*3 + width%4;

	if((bmp_line_data = (unsigned char *)malloc(sizeof(unsigned char)*real_width)) == NULL){
		fprintf(stderr, "Error: Allocation error.\n");
		return NULL;
	}

	if((img = Create_Image(width, height)) == NULL){
		free(bmp_line_data);
		fclose(fp);
		return NULL;
	}

	for(i=0; i<height; i++){
		fread(bmp_line_data, 1, real_width, fp);
		for(j=0; j<width; j++){
			img->data[(height-i-1)*width + j].b = bmp_line_data[j*3];
			img->data[(height-i-1)*width + j].g = bmp_line_data[j*3 + 1];
			img->data[(height-i-1)*width + j].r = bmp_line_data[j*3 + 2];
		}
	}

	free(bmp_line_data);

	fclose(fp);
#else

	unsigned char *data = 0;
	int x, y;
	int nbit;
	data = stbi_load(filename, &x, &y, &nbit, 0);
	if (data == NULL)
	{
		printf("image file[%s] read error.\n", filename);
		return NULL;
	}
	printf("height %d   width %d \n", y, x);

	height = y;
	width = x;
	if ((img = Create_Image(width, height)) == NULL)
	{
		return NULL;
	}

	for (i = 0; i<height; i++){
		for (j = 0; j<width; j++){
			if (nbit == 1)	//8bit
			{
				int pos = (i*width + j);
				img->data[pos].r = data[pos];
				img->data[pos].g = data[pos];
				img->data[pos].b = data[pos];
				//img->data[pos].alp = 255;
			}
			if (nbit == 2)	//16bit
			{
				int pos = (i*width + j);
				img->data[pos].r = data[pos * 2 + 0];
				img->data[pos].g = data[pos * 2 + 0];
				img->data[pos].b = data[pos * 2 + 0];
				//img->data[pos].alp = data[pos * 2 + 1];
			}
			if (nbit == 3)	//24
			{
				int pos = (i*width + j);
				img->data[pos].r = data[pos * 3 + 0];
				img->data[pos].g = data[pos * 3 + 1];
				img->data[pos].b = data[pos * 3 + 2];
				//img->data[pos].alp = 255;
			}
			if (nbit == 4)	//32
			{
				int pos = (i*width + j);
				img->data[pos].r = data[pos * 4 + 0];
				img->data[pos].g = data[pos * 4 + 1];
				img->data[pos].b = data[pos * 4 + 2];
				//img->data[pos].alp = data[pos * 4 + 3];
			}
		}
	}
	stbi_image_free(data);

#endif
	return img;
}


void BitMap::Offset( int size )
{
	if ( size == 0 ) return;

	Image* d = Create_Image( data->width+size*2, data->height+size*2);

	if ( size > 0 )
	{
		int ii = 0;
		int jj = 0;
		for(int i=size; i< d->height-size; i++, ii++){
			jj = 0;
			for(int j=size; j<d->width-size; j++, jj++){
				d->data[i*d->width + j].b = data->data[ii*data->width + jj].b;
				d->data[i*d->width + j].g = data->data[ii*data->width + jj].g;
				d->data[i*d->width + j].r = data->data[ii*data->width + jj].r;
			}
		}
	}else
	{
		size = -size;
		int ii = size;
		int jj = 0;
		for(int i=0; i< d->height; i++, ii++){
			jj = 0;
			for(int j=0; j<d->width; j++, jj++){
				d->data[i*d->width + j].b = data->data[ii*data->width + jj+size].b;
				d->data[i*d->width + j].g = data->data[ii*data->width + jj+size].g;
				d->data[i*d->width + j].r = data->data[ii*data->width + jj+size].r;
			}
		}
	}
	free(data->data);
	data = d;
}


int BitMap::Write_Bmp(char *filename, Image *img)
{
	int i, j;
	FILE *fp;
	int real_width;
	unsigned char **bmp_line_data;
	unsigned char header_buf[HEADERSIZE];
	unsigned int file_size;
	unsigned int offset_to_data;
	unsigned long info_header_size;
	unsigned int planes;
	unsigned int color;
	unsigned long compress;
	unsigned long data_size;
	long xppm;
	long yppm;

#if 10

	/// 画像の書き出し.
	char drive[256];
	char dir[256];
	char name[256];
	char ext[256];
	_splitpath(filename, drive, dir, name, ext);


	const int w = W();
	const int h = H();
	std::vector<unsigned char> image_buffer(w * h * 3);

	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			int pos = (y*w+x) * 3;
			unsigned char r;
			unsigned char g;
			unsigned char b;
			
			if (cell(x, y).r < 0) r = 0;
			else if (cell(x, y).r > 255) r = 255;
			else r = cell(x, y).r;

			if (cell(x, y).g < 0) g = 0;
			else if (cell(x, y).g > 255) g = 255;
			else g = cell(x, y).g;

			if (cell(x, y).b < 0) b = 0;
			else if (cell(x, y).b > 255) b = 255;
			else b = cell(x, y).b;

			image_buffer[pos + 0] = r;
			image_buffer[pos + 1] = g;
			image_buffer[pos + 2] = b;
		}
	}
	if ( stricmp(ext, ".png") == 0 )
	{
		stbi_write_png(filename, w, h, STBI_rgb, &(*image_buffer.begin()), 0);
	}else
	if ( stricmp(ext, ".bmp") == 0 )
	{
		stbi_write_bmp(filename, w, h, STBI_rgb, &(*image_buffer.begin()));
	}
	else
	if ( stricmp(ext, ".tga") == 0 )
	{
		stbi_write_tga(filename, w, h, STBI_rgb, &(*image_buffer.begin()));
	}else
	{
		stbi_write_bmp((std::string(filename) + std::string(".bmp")).c_str(), w, h, STBI_rgb, &(*image_buffer.begin()));
	}

	return 0;
#else
	if((fp = fopen(filename, "wb")) == NULL){
		fprintf(stderr, "Error: %s could not open.\n", filename);
		return 1;
	}

	real_width = img->width*3 + img->width%4;

	file_size = img->height * real_width + HEADERSIZE;
	offset_to_data = HEADERSIZE;
	info_header_size = INFOHEADERSIZE;
	planes = 1;
	color = 24;
	compress = 0;
	data_size = img->height * real_width;
	xppm = 1;
	yppm = 1;
	
	header_buf[0] = 'B';
	header_buf[1] = 'M';
	memcpy(header_buf + 2, &file_size, sizeof(file_size));
	header_buf[6] = 0;
	header_buf[7] = 0;
	header_buf[8] = 0;
	header_buf[9] = 0;
	memcpy(header_buf + 10, &offset_to_data, sizeof(file_size));
	header_buf[11] = 0;
	header_buf[12] = 0;
	header_buf[13] = 0;

	memcpy(header_buf + 14, &info_header_size, sizeof(info_header_size));
	header_buf[15] = 0;
	header_buf[16] = 0;
	header_buf[17] = 0;
	memcpy(header_buf + 18, &img->width, sizeof(img->width));
	memcpy(header_buf + 22, &img->height, sizeof(img->height));
	memcpy(header_buf + 26, &planes, sizeof(planes));
	memcpy(header_buf + 28, &color, sizeof(color));
	memcpy(header_buf + 30, &compress, sizeof(compress));
	memcpy(header_buf + 34, &data_size, sizeof(data_size));
	memcpy(header_buf + 38, &xppm, sizeof(xppm));
	memcpy(header_buf + 42, &yppm, sizeof(yppm));
	header_buf[46] = 0;
	header_buf[47] = 0;
	header_buf[48] = 0;
	header_buf[49] = 0;
	header_buf[50] = 0;
	header_buf[51] = 0;
	header_buf[52] = 0;
	header_buf[53] = 0;

	fwrite(header_buf, sizeof(unsigned char), HEADERSIZE, fp);
	
	const int hsz = img->height;
	const int wsz = img->width;

	if((bmp_line_data = (unsigned char **)malloc(sizeof(unsigned char*)*img->height)) == NULL){
		fprintf(stderr, "Error: Allocation error.\n");
		fclose(fp);
		return 1;
	}
	for ( int i = 0; i < hsz; i++ )
	{
		if((bmp_line_data[i] = (unsigned char *)malloc(sizeof(unsigned char)*real_width)) == NULL){
			fprintf(stderr, "Error: Allocation error.\n");
			fclose(fp);
			return 1;
		}
	}

#pragma omp parallel for OMP_SCHEDULE_BMP
	for(int i=0; i< hsz; i++){
		for(int j=0; j< wsz; j++){

			int b = (int)(img->data[(img->height - i - 1)*img->width + j].b  + 0.5);
			int g = (int)(img->data[(img->height - i - 1)*img->width + j].g  + 0.5);
			int r = (int)(img->data[(img->height - i - 1)*img->width + j].r  + 0.5);
			if (r < 0) r = 0;
			if (g < 0) g = 0;
			if (b < 0) b = 0;
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;

			bmp_line_data[i][j * 3] = b;
			bmp_line_data[i][j * 3 + 1] = g;
			bmp_line_data[i][j * 3 + 2] = r;
		}
		for(int j=img->width*3; j<real_width; j++){
			bmp_line_data[i][j] = 0;
		}
	}

	for(i=0; i< hsz; i++){
		fwrite(bmp_line_data[i], sizeof(unsigned char), real_width, fp);
		free(bmp_line_data[i]);
	}
	free(bmp_line_data);

	fclose(fp);
#endif

	return 0;
}

Image *BitMap::Create_Image(int width, int height)
{
	Image *img;

	if((img = (Image *)malloc(sizeof(Image))) == NULL){
		fprintf(stderr, "Allocation error\n");
		return NULL;
	}

	if((img->data = (Rgb*)malloc(sizeof(Rgb)*width*height)) == NULL){
		fprintf(stderr, "Allocation error\n");
		free(img);
		return NULL;
	}

	memset(img->data, '\0', sizeof(Rgb)*width*height);
	img->width = width;
	img->height = height;

	return img;
}

void BitMap::Free_Image(Image *img)
{
	if ( img->data ) free(img->data);
	if ( img ) free(img);
}


int BitMap::WriteText(char* filename, Image* img)
{
	int i, j;
	FILE *fp;
	int real_width;
	unsigned char *bmp_line_data;

	if((fp = fopen(filename, "w")) == NULL){
		fprintf(stderr, "Error: %s could not open.\n", filename);
		return 1;
	}

	real_width = img->width*3 + img->width%4;

	
	fprintf(fp, "%d %d\n", img->width, img->height);
	if((bmp_line_data = (unsigned char *)malloc(sizeof(unsigned char)*real_width)) == NULL){
		fprintf(stderr, "Error: Allocation error.\n");
		fclose(fp);
		return 1;
	}

	for(i=0; i<img->height; i++){
		for(j=0; j<img->width; j++){
			bmp_line_data[j*3]		=	img->data[(img->height - i - 1)*img->width + j].b;
			bmp_line_data[j*3 + 1]	=	img->data[(img->height - i - 1)*img->width + j].g;
			bmp_line_data[j*3 + 2]	=	img->data[(img->height - i - 1)*img->width + j].r;
		}
		for(j=img->width*3; j<real_width; j++){
			bmp_line_data[j] = 0;
		}

		for ( int k = 0; k < real_width; k++ )
		{
			fprintf(fp, "%d\n", bmp_line_data[k]);
		}
	}

	free(bmp_line_data);

	fclose(fp);

	return 0;
}

Image* BitMap::Read_Text(char *filename)
{
	int i, j;
	int real_width;
	unsigned int width, height;
	unsigned int color;
	FILE *fp;
	unsigned char header_buf[HEADERSIZE];
	unsigned char *bmp_line_data;
	Image *img;

	if((fp = fopen(filename, "r")) == NULL){
		fprintf(stderr, "Error: %s could not read.\n", filename);
		return NULL;
	}

	char* buf = new char[LINELENGMAX];
	fgets(buf, LINELENGMAX, fp);

	sscanf(buf, "%d %d", &width, &height);
	color = 24;

	real_width = width*3 + width%4;

	if((bmp_line_data = (unsigned char *)malloc(sizeof(unsigned char)*real_width)) == NULL){
		fprintf(stderr, "Error: Allocation error.\n");
		delete [] buf;
		return NULL;
	}

	if((img = Create_Image(width, height)) == NULL){
		free(bmp_line_data);
		fclose(fp);
		delete [] buf;
		return NULL;
	}

	for(i=0; i<height; i++){
		for ( int k = 0; k < real_width; k++ )
		{
			int dmy;
			fgets(buf, LINELENGMAX, fp);
			sscanf(buf, "%d",  &dmy);
			bmp_line_data[k] = (unsigned char)dmy;
		}
		for(j=0; j<width; j++){
			img->data[(height-i-1)*width + j].b = bmp_line_data[j*3];
			img->data[(height-i-1)*width + j].g = bmp_line_data[j*3 + 1];
			img->data[(height-i-1)*width + j].r = bmp_line_data[j*3 + 2];
		}
	}

	free(bmp_line_data);

	fclose(fp);
	delete [] buf;

	return img;
}




Image* BitMap::Read_Csv(char* filename, double min, double max)
{
	int i, j;
	int real_width;
	unsigned int width, height;
	unsigned int color;
	FILE *fp;
	unsigned char header_buf[HEADERSIZE];
	unsigned char *bmp_line_data;
	Image *img;

	if((fp = fopen(filename, "r")) == NULL){
		fprintf(stderr, "Error: %s could not read.\n", filename);
		return NULL;
	}

	char* buf = new char[LINELENGMAX];
	fgets(buf, LINELENGMAX, fp);
	color = 24;

	char* p = buf;
	width = 0;

	do{
		p = strchr(p, ',');
		if ( p )
		{
			width++;
			p++;
		}
	}while( p );
	width++;

	height = 1;
	while( fgets(buf, LINELENGMAX, fp) ) height++;


	fclose(fp);
	fp = fopen(filename, "r");

	real_width = width*3 + width%4;

	if((bmp_line_data = (unsigned char *)malloc(sizeof(unsigned char)*real_width)) == NULL){
		fprintf(stderr, "Error: Allocation error.\n");
		delete [] buf;
		return NULL;
	}
	memset(bmp_line_data, '\0', sizeof(unsigned char)*real_width);

	if((img = Create_Image(width, height)) == NULL){
		free(bmp_line_data);
		fclose(fp);
		delete [] buf;
		return NULL;
	}

	double w = 255.0/(max - min);
	for(i=0; i<height; i++){
		fgets(buf, LINELENGMAX, fp);
		p = buf;
		for ( int k = 0; k < width; k++ )
		{
			double dmy;
			sscanf(p, "%lf",  &dmy);
			dmy = (dmy - min)*w;
			if ( dmy < 0.0 ) dmy = 0.0;
			else if ( dmy > 255.0 ) dmy = 255.0;

			unsigned char* b = bmp_line_data + 3 * k;
			*(b++) = (unsigned char)dmy;
			*(b++) = (unsigned char)dmy;
			*(b++) = (unsigned char)dmy;
			//bmp_line_data[3*k] = (unsigned char)dmy;
			//bmp_line_data[3*k+1] = (unsigned char)dmy;
			//bmp_line_data[3*k+2] = (unsigned char)dmy;
			p = strchr(p, ',');
			if ( p ) p++;
		}
		//for(j=img->width*3; j<real_width; j++){
		//	bmp_line_data[j] = 0;
		//}

#pragma omp parallel for OMP_SCHEDULE_BMP
		for(int j=0; j<width; j++){
			img->data[(height-i-1)*width + j].b = bmp_line_data[j*3];
			img->data[(height-i-1)*width + j].g = bmp_line_data[j*3 + 1];
			img->data[(height-i-1)*width + j].r = bmp_line_data[j*3 + 2];
		}
	}

	free(bmp_line_data);

	fclose(fp);
	delete [] buf;


	return img;
}

Image* BitMap::Read_Csv(double* value, int w, int h, double min, double max)
{
	unsigned int width, height;
	unsigned int color;
	Image *img;

	color = 24;

	width = w;
	height = h;

	if((img = Create_Image(width, height)) == NULL){
		return NULL;
	}

	double ww = 255.0/(max - min);
#pragma omp parallel for OMP_SCHEDULE_BMP
	for(int i=0; i<height; i++)
	{
		for ( int k = 0; k < width; k++ )
		{
			double dmy = value[i*w+k];
			dmy = (dmy - min)*ww;
			if ( dmy < 0.0 ) dmy = 0.0;
			else if ( dmy > 255.0 ) dmy = 255.0;

			img->data[w*i+k] = Rgb_((unsigned char)dmy);
		}
	}
	return img;
}


int BitMap::Write_Csv(char* filename, Image* img, int rgb)
{
	int i, j;
	FILE *fp;
	int real_width;

	if((fp = fopen(filename, "w")) == NULL){
		fprintf(stderr, "Error: %s could not open.\n", filename);
		return 1;
	}

	real_width = img->width*3 + img->width%4;

	
	for(i=0; i<img->height; i++){
		switch(rgb)
		{
		case 0:
			for(j=0; j<img->width-1; j++){
				fprintf(fp, "%.1f ,", (float)img->data[(img->height - i - 1)*img->width + j].r);
			}
			fprintf(fp, "%.1f\n", (float)img->data[(img->height - i - 1)*img->width + j].r);
		break;
		case 1:
			for(j=0; j<img->width-1; j++){
				fprintf(fp, "%.1f ,", (float)img->data[(img->height - i - 1)*img->width + j].g);
			}
			fprintf(fp, "%.1f\n", (float)img->data[(img->height - i - 1)*img->width + j].g);
		break;
		case 2:
			for(j=0; j<img->width-1; j++){
				fprintf(fp, "%.1f ,", (float)img->data[(img->height - i - 1)*img->width + j].b);
			}
			fprintf(fp, "%.1f\n", (float)img->data[(img->height - i - 1)*img->width + j].b);
		break;
		}

	}
	fclose(fp);

	return 0;
}

int BitMap::Write_Csv(char* filename, Image* img, int rgb, double min, double max)
{
	int i, j;
	FILE *fp;
	int real_width;

	if((fp = fopen(filename, "w")) == NULL){
		fprintf(stderr, "Error: %s could not open.\n", filename);
		return 1;
	}

	real_width = img->width*3 + img->width%4;
	
	int minpix = 255;
	int maxpix = 0;
	int minpix2 = 255;
	int maxpix2 = 0;
	for (int i = 0; i < H(); i++)
	{
		for (int j = 0; j < W(); j++)
		{
			if ( minpix > cell(i,j).r ) minpix = cell(i,j).r;
			if ( maxpix < cell(i,j).r ) maxpix = cell(i,j).r;
			if ( minpix2 > cell(i,j).b ) minpix2 = cell(i,j).b;
			if ( maxpix2 < cell(i,j).b ) maxpix2 = cell(i,j).b;
		}
	}
	printf("pixcel(red) %d ～ %d\n", minpix, maxpix);
	printf("pixcel(blue) %d ～ %d\n", minpix2, maxpix2);
	if ( minpix == maxpix ) return 1;


	printf("0.0メートル付近は%dになります\n", (int)(-min*(float)(maxpix-minpix)/(max -min)+(float)minpix));

	float t;
	for(i=0; i<img->height; i++){
		switch(rgb)
		{
		case 0:
			for(j=0; j<img->width-1; j++){

				t = (float)(img->data[(img->height - i - 1)*img->width + j].r - minpix)/(float)(maxpix-minpix);
				fprintf(fp, "%.3f ,", min*(1.0-t) + max*t);
			}
			t = (float)(img->data[(img->height - i - 1)*img->width + j].r - minpix)/(float)(maxpix-minpix);
			fprintf(fp, "%.3f\n", min*(1.0-t) + max*t);
		break;
		case 1:
			for(j=0; j<img->width-1; j++){

				t = (float)(img->data[(img->height - i - 1)*img->width + j].g - minpix)/(float)(maxpix-minpix);
				fprintf(fp, "%.3f ,", min*(1.0-t) + max*t);
			}
			t = (float)(img->data[(img->height - i - 1)*img->width + j].g - minpix)/(float)(maxpix-minpix);
			fprintf(fp, "%.3f\n", min*(1.0-t) + max*t);
		break;
		case 2:
			for(j=0; j<img->width-1; j++){

				t = (float)(img->data[(img->height - i - 1)*img->width + j].b - minpix)/(float)(maxpix-minpix);
				fprintf(fp, "%.3f ,", min*(1.0-t) + max*t);
			}
			t = (float)(img->data[(img->height - i - 1)*img->width + j].b - minpix)/(float)(maxpix-minpix);
			fprintf(fp, "%.3f\n", min*(1.0-t) + max*t);
		break;
		case 3:
			for(j=0; j<img->width-1; j++){

				if ( img->data[(img->height - i - 1)*img->width + j].b == 0 )
				{
					t = (float)(img->data[(img->height - i - 1)*img->width + j].r)/255.0f;
					fprintf(fp, "%.3f ,", max*t);
				}else
				{
					t = (float)(img->data[(img->height - i - 1)*img->width + j].b)/255.0f;
					fprintf(fp, "%.3f ,", min*t);
				}
			}
			if ( img->data[(img->height - i - 1)*img->width + j].b == 0 )
			{
				t = (float)(img->data[(img->height - i - 1)*img->width + j].r)/255.0f;
				fprintf(fp, "%.3f\n", max*t);
			}else
			{
				t = (float)(img->data[(img->height - i - 1)*img->width + j].b)/255.0f;
				fprintf(fp, "%.3f\n", min*t);
			}
		break;
		case 4:
			for(j=0; j<img->width-1; j++){

				if ( img->data[(img->height - i - 1)*img->width + j].b == 0 )
				{
					t = (float)(img->data[(img->height - i - 1)*img->width + j].r)/(float)(maxpix);;
					fprintf(fp, "%.3f ,", max*t);
				}else
				{
					t = (float)(img->data[(img->height - i - 1)*img->width + j].b)/(float)(maxpix2);
					fprintf(fp, "%.3f ,", min*t);
				}
			}
			if ( img->data[(img->height - i - 1)*img->width + j].b == 0 )
			{
				t = (float)(img->data[(img->height - i - 1)*img->width + j].r)/(float)(maxpix);
				fprintf(fp, "%.3f\n", max*t);
			}else
			{
				t = (float)(img->data[(img->height - i - 1)*img->width + j].b)/(float)(maxpix2);
				fprintf(fp, "%.3f\n", min*t);
			}
		break;
		case 5:
			float alp = 2.0f;
			char* alp_env = getenv("BMP2CSV_CONV_COEF");
			if ( alp_env )
			{
				alp = atof(alp_env);
			}
			for(j=0; j<img->width-1; j++){

				if ( img->data[(img->height - i - 1)*img->width + j].b == 0 )
				{
					t = (float)(img->data[(img->height - i - 1)*img->width + j].r)/(float)(maxpix);;
					t = t*exp(alp*t*t)/exp(alp);
					t *= (float)(img->data[(img->height - i - 1)*img->width + j].r)/255.0f;

					fprintf(fp, "%.3f ,", max*t);
				}else
				{
					t = (float)(img->data[(img->height - i - 1)*img->width + j].b)/255.0f;;
					fprintf(fp, "%.3f ,", min*t);
				}
			}
			if ( img->data[(img->height - i - 1)*img->width + j].b == 0 )
			{
				t = (float)(img->data[(img->height - i - 1)*img->width + j].r)/(float)(maxpix);
				t =t*exp(alp*t*t)/exp(alp);
				t *= (float)(img->data[(img->height - i - 1)*img->width + j].r)/255.0f;
				fprintf(fp, "%.3f\n", max*t);
			}else
			{
				t = (float)(img->data[(img->height - i - 1)*img->width + j].b)/255.0f;;
				fprintf(fp, "%.3f\n", min*t);
			}
		break;
		}
	}
	fclose(fp);

	return 0;
}

void BitMap::Reverse()
{
#pragma omp parallel for OMP_SCHEDULE_BMP
	for(int i=0; i<data->height; i++)
	{
		for(int j=0; j<data->width; j++)
		{
			Rgb& rgb = cell(i,j);
			rgb.r = 255 - rgb.r;
			rgb.b = rgb.r;
			rgb.g = rgb.r;
		}
	}
}

void BitMap::convolve_smooth(int *mask, double conv[3][3])
{
	double tmp[3];


	Image* img = Create_Image(data->width, data->height);

	/*畳み込む分布の宣言。ここを変えるとスムージンングでなくなる
	この場合は３x３画素の平均が真ん中に入るようにつくってある。*/

	for(int j=1;j<data->height-1;j++)
	{
		for(int i=1;i<data->width-1;i++)
		{
			if ( mask )
			{
				if ( mask[(data->height-1-j)*data->width+i] == 0 )
				{
					img->data[j*data->width+i].r = data->data[j*data->width+i].r;
					img->data[j*data->width+i].g = data->data[j*data->width+i].g;
					img->data[j*data->width+i].b = data->data[j*data->width+i].b;
					continue;
				}
			}

			/*まずは、９つの値をよび、それぞれに畳み込む分布を掛け
			そのすべての値を加算した結果を tmpに保存する*/

    
			tmp[0] =
			data->data[(j-1)*data->width+(i-1)].r*conv[0][0]+
			data->data[(j-1)*data->width+i    ].r*conv[0][1]+
			data->data[(j-1)*data->width+(i+1)].r*conv[0][2]+

			data->data[j  *data->width+(i-1)].r*conv[1][0]+
			data->data[j  *data->width+i    ].r*conv[1][1]+
			data->data[j  *data->width+(i+1)].r*conv[1][2]+

			data->data[(j+1)*data->width+(i-1)].r*conv[2][0]+
			data->data[(j+1)*data->width+i    ].r*conv[2][1]+
			data->data[(j+1)*data->width+(i+1)].r*conv[2][2];

			tmp[1] =
			data->data[(j-1)*data->width+(i-1)].g*conv[0][0]+
			data->data[(j-1)*data->width+i    ].g*conv[0][1]+
			data->data[(j-1)*data->width+(i+1)].g*conv[0][2]+

			data->data[j  *data->width+(i-1)].g*conv[1][0]+
			data->data[j  *data->width+i    ].g*conv[1][1]+
			data->data[j  *data->width+(i+1)].g*conv[1][2]+

			data->data[(j+1)*data->width+(i-1)].g*conv[2][0]+
			data->data[(j+1)*data->width+i    ].g*conv[2][1]+
			data->data[(j+1)*data->width+(i+1)].g*conv[2][2];

			tmp[2] =
			data->data[(j-1)*data->width+(i-1)].g*conv[0][0]+
			data->data[(j-1)*data->width+i    ].g*conv[0][1]+
			data->data[(j-1)*data->width+(i+1)].g*conv[0][2]+

			data->data[j  *data->width+(i-1)].b*conv[1][0]+
			data->data[j  *data->width+i    ].b*conv[1][1]+
			data->data[j  *data->width+(i+1)].b*conv[1][2]+

			data->data[(j+1)*data->width+(i-1)].b*conv[2][0]+
			data->data[(j+1)*data->width+i    ].b*conv[2][1]+
			data->data[(j+1)*data->width+(i+1)].b*conv[2][2];

			///*９画素の平均なので,９で割る*/
			//tmp[0] /= 9.0;
			//tmp[1] /= 9.0;
			//tmp[2] /= 9.0;

			if ( tmp[0] < 0 ) tmp[0] = 0;
			if ( tmp[0] > 255 ) tmp[0] = 255;
			if ( tmp[1] < 0 ) tmp[1] = 0;
			if ( tmp[1] > 255 ) tmp[1] = 255;
			if ( tmp[2] < 0 ) tmp[2] = 0;
			if ( tmp[2] > 255 ) tmp[2] = 255;

			/*割った値を代入*/
			img->data[j*data->width+i].r = (unsigned char)tmp[0];
			img->data[j*data->width+i].g = (unsigned char)tmp[1];
			img->data[j*data->width+i].b = (unsigned char)tmp[2];
		}
	}

	free( data->data);
	data->data = img->data;

	img->data = NULL;
	Free_Image( img );
}

void BitMap::ColorTable(BitMap& colormap)
{
	printf("%d %d\n", colormap.W(), colormap.H());

	int n = colormap.H();
	if ( n >= 1024 ) n = 1024;
	int w = colormap.W()/2;

	colortableNum = n-1;
	printf("colortableNum %d\n", colortableNum);
	for ( int i = 0; i < n; i++ )
	{
		colorTbl[i][0] = colormap.cell(i, w).r;
		colorTbl[i][1] = colormap.cell(i, w).g;
		colorTbl[i][2] = colormap.cell(i, w).b;
		//printf("(%d) %d %d %d\n", i, colorTbl[i][0], colorTbl[i][1], colorTbl[i][2]);
	}
}

void BitMap::ColorTable()
{
	colortableNum = 255;
	// 青→緑 
	for (int i = 0; i < 64; i++) 
	{ 
		int green = i * 4;
		colorTbl[i][0] = 0;
		colorTbl[i][1] = green;
		colorTbl[i][2] = 255 - green;
	} 
	// 緑→黄 
	for (int i = 0; i < 64; i++) 
	{ 
		int red = i * 4;
		colorTbl[i+64][0] = red;
		colorTbl[i+64][1] = 255;
		colorTbl[i+64][2] = 0;
	} 
	// 黄→赤 
	for (int i = 0; i < 128; i++) 
	{ 
		int green = 255 - i * 2; 
		colorTbl[i+128][0] = 255;
		colorTbl[i+128][1] = green;
		colorTbl[i+128][2] = 0;
	}
}
void BitMap::ColorTable(int startIndex, int endIndex, unsigned char start[3], unsigned char end[3])
{
	colortableNum = 255;
	double dr = ((double)(end[0] - start[0])/(double)(endIndex-startIndex));
	double dg = ((double)(end[1] - start[1])/(double)(endIndex-startIndex));
	double db = ((double)(end[2] - start[2])/(double)(endIndex-startIndex));

	for ( int  i = startIndex; i <= endIndex; i++ )
	{
		int ir = (int)floor( start[0] + i*dr );
		int ig = (int)floor( start[1] + i*dg );
		int ib = (int)floor( start[2] + i*db );

		colorTbl[i][0] = ir;
		colorTbl[i][1] = ig;
		colorTbl[i][2] = ib;
	}
}

void BitMap::ColorLevel( double min, double max, double* z, double zmask, unsigned char* maskcolor, int* top, double* elv)
{

#ifdef BMP_USE_ACC
	const float b = 40.f;
	const float w = colortableNum/(max - min);
	
	const float zero = (0.0f-min)/(max-min);

	const int y = data->height;
	const int x = data->width;
	Rgb* p = data->data;

	const int colortableNum_ = colortableNum;
	const float zmask_f = zmask;
	const float min_f = min;
	const float max_f = max;

	int maskcolor_flg = 0;
	if ( maskcolor )
	{
		maskcolor_flg = 1;
	}


	float* z_f = new float[x*y];
#pragma omp parallel for OMP_SCHEDULE_BMP
	for ( int i = 0; i < y; i++ )
	{
		for ( int j = 0; j < x; j++ )
		{
			z_f[i*x+j] = z[(y-i-1)*x + j];
		}
	}

	// GPU側r CPU側rr array_viewでrがデストラクトされる時に r->rrへの転送があるため{}内で括る
	int* rr= new int[x*y];
	{
		extent<2> e(y,x);
		array_view<int, 2> r(e, rr);
		array_view<float, 2> z_f_(e, z_f);


		parallel_for_each(e, [=](index<2> idx) BMP_ACC_RESTRICTION
		{
			const int i = idx[0];
			const int j = idx[1];
			float zz = z_f_[idx];

			r[idx] = -2;
			if ( FAST_MATH fabs(zz) < zmask_f )
			{
				if ( maskcolor_flg )
				{
					r[idx] = -1;
				}
			}else
			{
				float dmy  = (zz-min_f)*w;

				if ( dmy <= zero )
				{
					dmy = ((zz - min_f)/(0.0f - min_f))*b;
				}else
				{
					dmy = b+((zz - 0.0f)/(max_f - 0.0f))*(colortableNum_ - b-1.0f);
				}

				if ( dmy < 0.0f ) dmy = 0.0f;
				if ( dmy > (float)colortableNum_ ) dmy = (float)colortableNum_;

				int colorindex = (int)dmy;
				if ( colorindex < 0 ) colorindex = 0;
				else if ( colorindex >= COLORTBL_MAX) colorindex = COLORTBL_MAX-1;
				r[idx] = colorindex;
			}
		});
	}
	delete [] z_f;

#pragma omp parallel for OMP_SCHEDULE_BMP
	for ( int i = 0; i <x*y; i++ )
	{
		if ( rr[i] == -2 ) continue;
		if ( maskcolor_flg ) p[i] = Rgb(maskcolor);
		if ( rr[i] >= 0 )
		{
			p[i] = Rgb(colorTbl[rr[i]]);
		}
	}
	delete [] rr;

#else
	const double b = 40;
	const double w = colortableNum/(max - min);
	
	const double zero = (0.0-min)/(max-min);

	const int y = data->height;
	const int x = data->width;
	Rgb* p = data->data;

	const int colortableNum_ = colortableNum;
	const unsigned char* colorTbl_ = &(colorTbl[0][0]);

#pragma omp parallel for OMP_SCHEDULE_BMP
	for ( int i = 0; i < y; i++ )
	{
		for ( int j = 0; j < x; j++ )
		{
			double zz = z[(y-i-1)*x + j];
			//if ( top && top[(data->height-i-1)*data->width + j] == 0 && elv )
			//{
			//	zz = zz - elv[(data->height-i-1)*data->width + j];
			//}
			if ( fabs(zz) < zmask )
			{
				if ( maskcolor )
				{
					p[i*x + j] = Rgb(maskcolor);
				}
				continue;
			}

			double dmy  = (zz-min)*w;

			if ( dmy <= zero )
			{
				dmy = ((zz - min)/(0.0 - min))*b;
			}else
			{
				dmy = b+((zz - 0.0)/(max - 0.0))*(colortableNum_ - b-1);
			}

			if ( dmy < 0.0 ) dmy = 0.0;
			if ( dmy > colortableNum_ ) dmy = colortableNum_;

			int colorindex = (int)dmy;
			if ( colorindex < 0 ) colorindex = 0;
			else if ( colorindex >= COLORTBL_MAX) colorindex = COLORTBL_MAX-1;
			p[i*x + j] = Rgb( ((unsigned char(*)[3])colorTbl_)[colorindex]);
		}
	}
#endif
}