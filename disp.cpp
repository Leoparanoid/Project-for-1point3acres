/*   CS580 HW1 display functions to be completed   */

#include   "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"

#define RGBbytes	3

//function to check Res
int BoundCheck(int res, int max);

//function to check GzIntensity
GzIntensity BoundCheck(GzIntensity rgb);

//function to shift GzIntensity by 4
char Shift4PPM(GzIntensity rgb);

int GzNewFrameBuffer(char** framebuffer, int width, int height){

/* HW1.1 create a framebuffer for MS Windows display:
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- pass back pointer 
 */

	* framebuffer = new char[RGBbytes * width * height];
	if((* framebuffer) == NULL) return GZ_FAILURE;

	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes){

/* HW1.2 create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/

	* display = new GzDisplay;
	if ((*display) == NULL) return GZ_FAILURE;
	//GzDisplay includes unsigned short xres, yres, GzPixel fbuf
	(* display)->xres = BoundCheck(xRes, MAXXRES);
	(* display)->yres = BoundCheck(yRes, MAXYRES);
	(* display)->fbuf = new GzPixel[(*display)->xres * (*display)->yres];
	if ((* display)->fbuf == NULL) return GZ_FAILURE;
	
	return GZ_SUCCESS;
}

int GzFreeDisplay(GzDisplay	*display){

/* HW1.3 clean up, free memory */

	if (display == NULL) return GZ_FAILURE;
	//Clean fbuf of display then display itself
	delete[] display->fbuf;
	delete display;
	return GZ_SUCCESS;
}

int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes){

/* HW1.4 pass back values for a display */
	if (display == NULL) return GZ_FAILURE;
	* xRes = display->xres;
	* yRes = display->yres;
	return GZ_SUCCESS;
}

int GzInitDisplay(GzDisplay	*display){

/* HW1.5 set everything to some default values - start a new frame */
	if (display == NULL || display->fbuf == NULL) return GZ_FAILURE;

	for (int j = 0; display->yres > j; j++) {
		for (int i = 0; display->xres > i; i++) {
			display->fbuf[ARRAY(i, j)].red = 2000;
			display->fbuf[ARRAY(i, j)].green = 2000;
			display->fbuf[ARRAY(i, j)].blue = 2000;
			display->fbuf[ARRAY(i, j)].alpha = 1;
			display->fbuf[ARRAY(i, j)].z = MAXINT;
		}
	}
	return GZ_SUCCESS;
}

int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z){

/* HW1.6 write pixel values into the display */
	if (display == NULL || display->fbuf == NULL) return GZ_FAILURE;

	//Bound Check:: Criticl Point: i & j cant be as same as display->xres & display->yres
	if (display->xres > i && display->yres > j && i >= 0 && j >= 0) {
		display->fbuf[ARRAY(i, j)].red = BoundCheck(r);
		display->fbuf[ARRAY(i, j)].green = BoundCheck(g);
		display->fbuf[ARRAY(i, j)].blue = BoundCheck(b);
		display->fbuf[ARRAY(i, j)].alpha = a;
		display->fbuf[ARRAY(i, j)].z = z;
		return GZ_SUCCESS;
	}
	else
		return GZ_FAILURE;
}

int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z){

/* HW1.7 pass back a pixel value to the display */
	if (display == NULL || display->fbuf == NULL) return GZ_FAILURE;
	if (display->xres > i && display->yres > j && i >= 0 && j >= 0) {
		* r = display->fbuf[ARRAY(i, j)].red;
		* g = display->fbuf[ARRAY(i, j)].green;
		* b = display->fbuf[ARRAY(i, j)].blue;
		* a = display->fbuf[ARRAY(i, j)].alpha;
		* z = display->fbuf[ARRAY(i, j)].z;
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2File(FILE* outfile, GzDisplay *display){

/* HW1.8 write pixels to ppm file -- "P6 %d %d 255\r" */

	if (display == NULL || outfile == NULL || display->fbuf == NULL) return GZ_FAILURE;

	//print the header of ppm file
	fprintf(outfile, "P6 %d %d 255\r", display->xres, display->yres);

	//scan by raster order
	for (int j = 0; display->yres > j; j++) {
		for (int i = 0; display->xres > i; i++) {
			//write in the value after shifting
			fprintf(outfile, "%c%c%c", Shift4PPM(display->fbuf[ARRAY(i, j)].red), Shift4PPM(display->fbuf[ARRAY(i, j)].green), Shift4PPM(display->fbuf[ARRAY(i, j)].blue));
		}
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display){

/* HW1.9 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	if (display == NULL || framebuffer == NULL) return GZ_FAILURE;

	//scan by raster order
	for (int j = 0; display->yres > j; j++) {
		for (int i = 0; display->xres > i; i++) {
			//write in the value after shifting
			framebuffer[ARRAY(i, j) * 3] = Shift4PPM(display->fbuf[ARRAY(i, j)].blue);
			framebuffer[ARRAY(i, j) * 3 + 1] = Shift4PPM(display->fbuf[ARRAY(i, j)].green);
			framebuffer[ARRAY(i, j) * 3 + 2] = Shift4PPM(display->fbuf[ARRAY(i, j)].red);
		}
	}
	return GZ_SUCCESS;
}


int BoundCheck(int res, int max) {

	if (res > max) return max;
	else if (res < 0) return 0;
	else return res;
}

GzIntensity BoundCheck(GzIntensity rgb) {

	//MaxIntensity and MinIntensity define in gz.h
	if (rgb > MaxIntensity) return MaxIntensity;
	else if (rgb < MinIntensity) return MinIntensity;
	else return rgb;
}

char Shift4PPM(GzIntensity rgb) {

	return rgb >> 4;
}

GzDepth GzGetFbufZ(GzDisplay *display, int i, int j) {
	if (display->xres > i && display->yres > j && i >= 0 && j >= 0)
		return display->fbuf[ARRAY(i, j)].z;
}