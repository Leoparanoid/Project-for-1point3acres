#include "disp.h" /* include your own disp.h file (e.g. hw1)*/

/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#ifndef GZRENDER
#define GZRENDER
typedef struct {			/* define a renderer */
  GzDisplay		*display;
  GzCamera		camera;
  short		    matlevel;	        /* top of stack - current xform */
  GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
  GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
  GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
  GzColor		flatcolor;          /* color state for flat shaded triangles */
  int			interp_mode;
  int			numlights;
  GzLight		lights[MAX_LIGHTS];
  GzLight		ambientlight;
  GzColor		Ka, Kd, Ks;
  float		    spec;		/* specular power */
  GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
  float shift_offsetX;
  float shift_offsetY;
  float weight;
}  GzRender;
#endif

// Function declaration
// HW2
int GzNewRender(GzRender **render, GzDisplay *display);
int GzFreeRender(GzRender *render);
int GzBeginRender(GzRender	*render);
int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList);
int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList);

// HW3
int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix	matrix);
int GzPopMatrix(GzRender *render);

// HW5
int GzFreeTexture();

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);

int CreateXsp(GzRender *render);
int CreateXpi(GzRender *render);
int CreateXiw(GzRender *render);

//Vector Type
#ifndef  VERTEX
#define VERTEX
class Vertex {
public:
	float x, y, z;
	float nx, ny, nz;
	float tx, ty;
	Vertex() {};
	Vertex(GzCoord & v) {
		x = v[0];
		y = v[1];
		z = v[2];
	}
};
#endif // ! VECTOR

#ifndef  EDGE
#define EDGE
class Edge {
public:
	Vertex current;
	float slopex, slopez;
	float sloper, slopeg, slopeb;
	float slopetx, slopety;

	Edge() {};
	void InitEdge(Vertex vbegin, Vertex vend) {
		slopex = (vend.x - vbegin.x) / (vend.y - vbegin.y);
		slopez = (vend.z - vbegin.z) / (vend.y - vbegin.y);
		sloper = (vend.nx - vbegin.nx) / (vend.y - vbegin.y);
		slopeg = (vend.ny - vbegin.ny) / (vend.y - vbegin.y);
		slopeb = (vend.nz - vbegin.nz) / (vend.y - vbegin.y);
		slopetx = (vend.tx - vbegin.tx) / (vend.y - vbegin.y);
		slopety = (vend.ty - vbegin.ty) / (vend.y - vbegin.y);
	}
};
#endif // ! EDGE