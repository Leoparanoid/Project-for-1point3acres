/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
    unsigned char		pixel[3];
    unsigned char     dummy;
    char  		foo[8];
    int   		i, j;
    FILE			*fd;

    if (reset) {          /* open and load texture file */
        fd = fopen ("texture", "rb");
        if (fd == NULL) {
            fprintf (stderr, "texture file not found\n");
            exit(-1);
        }
        fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
        image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
        if (image == NULL) {
            fprintf (stderr, "malloc for texture image failed\n");
            exit(-1);
        }

        for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
            fread(pixel, sizeof(pixel), 1, fd);
            image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
            image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
            image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
	    }

        reset = 0;          /* init is done */
	    fclose(fd);
    }
/* bounds-test u,v to make sure nothing will overflow image array bounds */
    if (u > 1) u = 1;
    else if (u < 0) u = 0;
    if (v > 1) v = 1;
    else if (v < 0) v = 0;

/* determine texture cell corner values and perform bilinear interpolation */
	float px = u * (xs - 1);
	float py = v * (ys - 1);
	float s = px - floor(px);
	float t = py - floor(py);

/* set color to interpolated GzColor value and return */
	for (int i = 0; i < 3; i++)
		color[i] = (1 - s) * (1 - t) * image[xs * (int)floor(py) + (int)floor(px)][i] +
				   (1 - s) * t * image[xs * (int)ceil(py) + (int)floor(px)][i] +
				   s * (1 - t) * image[xs * (int)floor(py) + (int)ceil(px)][i] +
				   s * t * image[xs * (int)ceil(py) + (int)ceil(px)][i] ;

	return GZ_SUCCESS;
}

class complexNumber{
	public:
		float r;
		float i;
		complexNumber() {};
		~complexNumber() {};
};


/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	int N = 100;
	complexNumber x;	complexNumber c;

	c.r = -0.7;    c.i = 0.4;
	x.r = 2 * u - 0.9;    x.i = 2 *v - 0.9;

	int i;
	float clr;
    for(i = 0; i < N; i++) 
	{
		if ((x.r * x.r + x.i * x.i) > 5.0)
			break;
		float tempr = x.r * x.r - x.i * x.i + c.r;
		float tempi = x.r * x.i + x.i * x.r + c.i;
        x.r = tempr;
        x.i = tempi;
    }

	if (i == N)
		for(int j = 0; j < 3; j++)
			color[j] = sqrt(x.r * x.r + x.i * x.i) / 2;
	else
		for(int j = 0; j < 3; j++)
			color[j] = ((float)i) / N * 10;

    return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

