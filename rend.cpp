/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#include <algorithm>
#include <vector>

#define PI 3.1415926535;
using namespace std;

// HW2
void WorldStoScreenS(vector<Vertex> & vert, GzMatrix matrix);
int ScreenCheck(GzDisplay *display, vector<Vertex> & v);
void GzMatrixMul(GzMatrix matrix1, GzMatrix matrix2, GzMatrix & matrix_out);

void InitVertex(GzCoord* coordPointer, vector<Vertex> & vert);
void InitVertex_t(GzTextureIndex* coordPointer_t, vector<Vertex> & verts);
bool SortVertexY(Vertex & v1, Vertex & v2);
void SetEdge(vector<Vertex> & vert, vector<Edge> & edges);
int TriShape(vector<Vertex> & vert, vector<Edge> & edges);
void RasterTri(GzRender *render, vector<Vertex> & vert, vector<Edge>& edges, int trishape);
void InitSpanLine(Edge & span, vector<Edge> & edges, int i, int j);
void MoveEdgeY(Edge & edge);
void MoveSpanX(Edge & span);
void RastLine(GzRender *render, Edge & span, vector<Edge> & edges, int i, int j);

float Dot(Vertex& v1, Vertex& v2);
Vertex Cross(Vertex& v1, Vertex& v2);
Vertex Normal(Vertex v1);
Vertex Normal_N(Vertex v1);
void Normal_C(GzCoord & c1);

Vertex VertexMerge(float a, Vertex & v1, float b, Vertex & v2);
Vertex VertexMerge(float a, GzCoord & c1, float b, GzCoord & c2);

int Color(GzRender *render, const Vertex v, GzColor color);

void v_perspective(vector<Vertex> &v);
void uv_perspective(Vertex &v, GzTextureIndex tex);

void shiftVerts(GzRender *render, vector<Vertex>& v);

/* NOT part of API - just for general assistance */

short ctoi(float color)		/* convert float color to GzIntensity short */
{
	return(short)((int)(color * ((1 << 12) - 1)));
}

int GzNewRender(GzRender **render, GzDisplay *display)
{
	/*
	- malloc a renderer struct
	- span interpolator needs pointer to display for pixel writes
	*/
	* render = new GzRender();
	(*render)->display = display;

	if ((*render) == NULL) return GZ_FAILURE;

	(*render)->tex_fun = NULL;
	(*render)->matlevel = 0;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			(*render)->Xsp[i][j] = 0;
			(*render)->camera.Xiw[i][j] = 0;
			(*render)->camera.Xpi[i][j] = 0;
		}

	(*render)->interp_mode = GZ_RGB_COLOR;
	GzColor Ka = DEFAULT_AMBIENT;
	GzColor Kd = DEFAULT_DIFFUSE;
	GzColor Ks = DEFAULT_SPECULAR;

	(*render)->camera.FOV = DEFAULT_FOV;

	for (int i = 0; i < 3; i++) {
		(*render)->camera.lookat[i] = 0;
		(*render)->camera.worldup[i] = 0;
		(*render)->Ka[i] = Ka[i];
		(*render)->Kd[i] = Kd[i];
		(*render)->Ks[i] = Ks[i];
	}
	(*render)->camera.worldup[1] = 1;
	(*render)->camera.position[0] = DEFAULT_IM_X;
	(*render)->camera.position[1] = DEFAULT_IM_Y;
	(*render)->camera.position[2] = DEFAULT_IM_Z;

	return GZ_SUCCESS;
}

int GzFreeRender(GzRender *render)
{
	/*
	-free all renderer resources
	*/
	if (render == NULL) return GZ_FAILURE;
	delete render;

	return GZ_SUCCESS;
}

int GzBeginRender(GzRender	*render)
{
	/*
	- set up for start of each frame - init frame buffer
	*/
	if (render == NULL) return GZ_FAILURE;
	GzInitDisplay(render->display);

	CreateXsp(render);	CreateXpi(render);	CreateXiw(render);
	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render, render->camera.Xpi);
	GzPushMatrix(render, render->camera.Xiw);

	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, GzPointer *valueList) /* void** valuelist */
{
	/*
	- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	- later set shaders, interpolaters, texture maps, and lights
	*/
	if (render == NULL) return GZ_FAILURE;
	for (int i = 0; i < numAttributes; i++) {
		switch (nameList[i]){
			case GZ_RGB_COLOR: {
				GzColor* colorPointer = (GzColor*)(valueList[i]);
				for (int j = 0; j < 3; j++)
					render->flatcolor[j] = (*colorPointer)[j];
				break;
			}
			case GZ_INTERPOLATE: {
				render->interp_mode = *((int*)(valueList[i]));
				//render->interp_mode = (int)(valueList[i]);
				//render->interp_mode = GZ_NORMALS;
				//render->interp_mode = GZ_COLOR;
				break;
			}
			case GZ_DIRECTIONAL_LIGHT: {
				for (int j = 0; j < 3; j++) {
					render->lights[render->numlights].color[j] = ((GzLight*)(valueList[i]))->color[j];
					render->lights[render->numlights].direction[j] = ((GzLight*)(valueList[i]))->direction[j];
				}
				Normal_C(render->lights[render->numlights].direction);
				render->numlights++;
				break;
			}
			case GZ_AMBIENT_LIGHT: {
				for (int j = 0; j < 3; j++)		
					render->ambientlight.color[j] = ((GzLight*)(valueList[i]))->color[j];
				break;
			}
			case GZ_AMBIENT_COEFFICIENT: {
				for (int j = 0; j < 3; j++)
					render->Ka[j] = (*(GzColor*)(valueList[i]))[j];
				break;
			}
			case GZ_DIFFUSE_COEFFICIENT: {
				for (int j = 0; j < 3; j++)
					render->Kd[j] = (*(GzColor*)(valueList[i]))[j];
				break;
			}
			case GZ_SPECULAR_COEFFICIENT: {
				for (int j = 0; j < 3; j++)
					render->Ks[j] = (*(GzColor*)(valueList[i]))[j];
				break;
			}
			case GZ_DISTRIBUTION_COEFFICIENT: {
				render->spec = (*(float*)(valueList[i]));
				break;
			}
			case GZ_TEXTURE_MAP: {
				render->tex_fun = (GzTexture)(valueList[i]);
				break;
			}
			case GZ_AASHIFTX: {
				render->shift_offsetX = (*(float*)(valueList[i]));
				break;
			}
			case GZ_AASHIFTY: {
				render->shift_offsetY = (*(float*)(valueList[i]));
				break;
			}
		}
	}
	return GZ_SUCCESS;
}

int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
	/*
	- pass in a triangle description with tokens and values corresponding to
	GZ_NULL_TOKEN:		do nothing - no values
	GZ_POSITION:		3 vert positions
	- Invoke the scan converter and return an error code
	*/
	//valueListTriangle[0] = (GzPointer)vertexList; 
	int trishape;
	GzCoord* coordPointer = (GzCoord *)(valueList[0]);
	GzCoord* coordPointer_n = (GzCoord *)(valueList[1]);
	GzTextureIndex* coordPointer_t = (GzTextureIndex*)(valueList[2]);

	vector<Vertex> verts(3);
	vector<Vertex> vnormals(3);
	vector<Vertex> vtextures(3);
	vector<Edge> edges(3);

	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_POSITION) {
			InitVertex(coordPointer, verts);
			WorldStoScreenS(verts, render->Ximage[render->matlevel - 1]);
			shiftVerts(render, verts);
		}

		else if (nameList[i] == GZ_NORMAL) {
			InitVertex(coordPointer_n, vnormals);
			WorldStoScreenS(vnormals, render->Xnorm[render->matlevel - 1]);		
			for (int j = 0; j < 3; j++) {
				Vertex temp = Normal(vnormals[j]);
				verts[j].nx = temp.x;
				verts[j].ny = temp.y;
				verts[j].nz = temp.z;
			}
		}
		else if (nameList[i] == GZ_TEXTURE_INDEX) {
			InitVertex_t(coordPointer_t, vtextures);
			for (int j = 0; j < 3; j++) {
				verts[j].tx = vtextures[j].tx;
				verts[j].ty = vtextures[j].ty;
			}
		}
	}
	//check if in the screen
	if (ScreenCheck(render->display, verts) == 0) {
		//sort verts
		GzColor color;
		sort(verts.begin(), verts.end(), SortVertexY);
		if (render->interp_mode == GZ_COLOR) {
			for (int j = 0; j < 3; j++) {
				Color(render, verts[j], color);
				verts[j].nx = color[0];
				verts[j].ny = color[1];
				verts[j].nz = color[2];
			}
		}
		v_perspective(verts);
		SetEdge(verts, edges);
		RasterTri(render, verts, edges, TriShape(verts, edges));
	}
	return GZ_SUCCESS;
}

void WorldStoScreenS(vector<Vertex>& vert, GzMatrix matrix) {
	GzMatrix temp;
	GzMatrix out;

	for (int i = 0; i < 3; i++) {
		temp[0][i] = vert[i].x;
		temp[1][i] = vert[i].y;
		temp[2][i] = vert[i].z;
		temp[3][i] = 1;
	}

	GzMatrixMul(matrix, temp, out);

	for (int i = 0; i < 3; i++) {
		vert[i].x = out[0][i] / out[3][i];
		vert[i].y = out[1][i] / out[3][i];
		vert[i].z = out[2][i] / out[3][i];
	}
}

void InitVertex(GzCoord* coordPointer, vector<Vertex>& verts) {
	verts.clear();
	for (int i = 0; i < 3; i++)
		verts.push_back(Vertex(*(coordPointer + i)));
}

void InitVertex_t(GzTextureIndex* coordPointer_t, vector<Vertex> & verts) {
	verts.clear();
	Vertex tmp;
	for (int i = 0; i < 3; i++) {
		tmp.tx = (*(coordPointer_t + i))[0];
		tmp.ty = (*(coordPointer_t + i))[1];
		verts.push_back(tmp);
	}
}

bool SortVertexY(Vertex& v1, Vertex& v2) {
	if (v1.y == v2.y)
		return v1.x < v2.x;
	return v1.y < v2.y;
}

void SetEdge(vector<Vertex>& vert, vector<Edge>& edges) {

	edges[0].InitEdge(vert[0], vert[1]);
	edges[1].InitEdge(vert[1], vert[2]);
	edges[2].InitEdge(vert[0], vert[2]);

	float dY;
	for (int i = 0; i < 3; i++) {
		if (i == 1) {
			dY = ceil(vert[1].y) - vert[1].y;
			edges[i].current.x = vert[i].x + edges[i].slopex * dY;
			edges[i].current.y = vert[i].y + dY;
			edges[i].current.z = vert[i].z + edges[i].slopez * dY;
			edges[i].current.nx = vert[i].nx + edges[i].sloper * dY;
			edges[i].current.ny = vert[i].ny + edges[i].slopeg * dY;
			edges[i].current.nz = vert[i].nz + edges[i].slopeb * dY;
			edges[i].current.tx = vert[i].tx + edges[i].slopetx * dY;
			edges[i].current.ty = vert[i].ty + edges[i].slopety * dY;
		}
		else {
			dY = ceil(vert[0].y) - vert[0].y;
			edges[i].current.x = vert[0].x + edges[i].slopex * dY;
			edges[i].current.y = vert[0].y + dY;
			edges[i].current.z = vert[0].z + edges[i].slopez * dY;
			edges[i].current.nx = vert[0].nx + edges[i].sloper * dY;
			edges[i].current.ny = vert[0].ny + edges[i].slopeg * dY;
			edges[i].current.nz = vert[0].nz + edges[i].slopeb * dY;
			edges[i].current.tx = vert[0].tx + edges[i].slopetx * dY;
			edges[i].current.ty = vert[0].ty + edges[i].slopety * dY;
		}
	}
}

int TriShape(vector<Vertex>& vert, vector<Edge>& edges) {
	if (edges[0].slopex < edges[2].slopex) // L
		return 0;
	else
		return 1;
}

void RasterTri(GzRender *render, vector<Vertex>& vert, vector<Edge>& edges, int trishape) {
	Edge spanline;

	//scan from vertex 1 >> vertex 2
	while (edges[0].current.y < vert[1].y) {
		if (trishape == 0) //L
			RastLine(render, spanline, edges, 0, 2);
		else if (trishape == 1) //R
			RastLine(render, spanline, edges, 2, 0);

		// move y : y++
		MoveEdgeY(edges[0]);
		MoveEdgeY(edges[2]);
	}

	//scan from vertex 2 >> vertex 3
	while (edges[1].current.y < vert[2].y) {
		if (trishape == 0) //L
			RastLine(render, spanline, edges, 1, 2);
		else if (trishape == 1) //R
			RastLine(render, spanline, edges, 2, 1);

		// move y : y++
		MoveEdgeY(edges[1]);
		MoveEdgeY(edges[2]);
	}
}

void InitSpanLine(Edge & span, vector<Edge>& edges, int i, int j) {
	float dX = ceil(edges[i].current.x) - edges[i].current.x;
	span.InitEdge(edges[i].current, edges[j].current);
	if (i > j) {
		span.slopez = (edges[i].current.z - edges[j].current.z) / (edges[i].current.x - edges[j].current.x);
		span.sloper = (edges[i].current.nx - edges[j].current.nx) / (edges[i].current.x - edges[j].current.x);
		span.slopeg = (edges[i].current.ny - edges[j].current.ny) / (edges[i].current.x - edges[j].current.x);
		span.slopeb = (edges[i].current.nz - edges[j].current.nz) / (edges[i].current.x - edges[j].current.x);
		span.slopetx = (edges[i].current.tx - edges[j].current.tx) / (edges[i].current.x - edges[j].current.x);
		span.slopety = (edges[i].current.ty - edges[j].current.ty) / (edges[i].current.x - edges[j].current.x);
	}
	else {
		span.slopez = (edges[j].current.z - edges[i].current.z) / (edges[j].current.x - edges[i].current.x);
		span.sloper = (edges[j].current.nx - edges[i].current.nx) / (edges[j].current.x - edges[i].current.x);
		span.slopeg = (edges[j].current.ny - edges[i].current.ny) / (edges[j].current.x - edges[i].current.x);
		span.slopeb = (edges[j].current.nz - edges[i].current.nz) / (edges[j].current.x - edges[i].current.x);
		span.slopetx = (edges[j].current.tx - edges[i].current.tx) / (edges[j].current.x - edges[i].current.x);
		span.slopety = (edges[j].current.ty - edges[i].current.ty) / (edges[j].current.x - edges[i].current.x);
	}
	span.current.x = edges[i].current.x + dX;
	span.current.y = edges[i].current.y;
	span.current.z = edges[i].current.z + span.slopez * dX;

	span.current.nx = edges[i].current.nx + span.sloper * dX;
	span.current.ny = edges[i].current.ny + span.slopeg * dX;
	span.current.nz = edges[i].current.nz + span.slopeb * dX;

	span.current.tx = edges[i].current.tx + span.slopetx * dX;
	span.current.ty = edges[i].current.ty + span.slopety * dX;
}

void MoveEdgeY(Edge & edge) {
	edge.current.x += edge.slopex;
	edge.current.y++;
	edge.current.z += edge.slopez;

	edge.current.nx += edge.sloper;
	edge.current.ny += edge.slopeg;
	edge.current.nz += edge.slopeb;

	edge.current.tx += edge.slopetx;
	edge.current.ty += edge.slopety;
}

void MoveSpanX(Edge & span) {
	span.current.x++;
	span.current.z += span.slopez;
	
	span.current.nx += span.sloper;
	span.current.ny += span.slopeg;
	span.current.nz += span.slopeb;

	span.current.tx += span.slopetx;
	span.current.ty += span.slopety;
}

void RastLine(GzRender *render, Edge & span, vector<Edge>& edges, int i, int j) {
	GzDepth fbZ;
	GzColor color;
	GzTextureIndex tex = {0, 0};

	InitSpanLine(span, edges, i, j);
	while (span.current.x < edges[j].current.x) {
		fbZ = GzGetFbufZ(render->display, (int)span.current.x, (int)edges[(i + j) % 2].current.y);
		if (span.current.z >= 0 && (fbZ == 0 || span.current.z < fbZ)) {
			if (render->interp_mode == GZ_COLOR) {
				uv_perspective(span.current, tex);
				render->tex_fun(tex[0], tex[1], color);
				GzPutDisplay(render->display, (int)span.current.x, (int)edges[(i + j) % 2].current.y, ctoi(color[0] * span.current.nx), ctoi(color[1] * span.current.ny), ctoi(color[2] * span.current.nz), 0, (GzDepth)span.current.z);
			}
			else if (render->interp_mode == GZ_NORMALS) {
				uv_perspective(span.current, tex);
				render->tex_fun(tex[0], tex[1], color);
				memcpy(&render->Ka, color, sizeof(GzColor));
				memcpy(&render->Kd, color, sizeof(GzColor));
				Color(render, Normal_N(span.current), color);
				GzPutDisplay(render->display, (int)span.current.x, (int)edges[(i + j) % 2].current.y, ctoi(color[0]), ctoi(color[1]), ctoi(color[2]), 0, (GzDepth)span.current.z);
			}
			else 
				GzPutDisplay(render->display, (int)span.current.x, (int)edges[(i + j) % 2].current.y, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 0, (GzDepth)span.current.z);
		}
		MoveSpanX(span);
	}
}

// T S R Matrix
int GzRotXMat(float degree, GzMatrix mat)
{
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value

	float radian = (degree / 180) * PI;
	mat[0][0] = 1;	mat[1][0] = 0;				mat[2][0] = 0;				mat[3][0] = 0;
	mat[0][1] = 0;	mat[1][1] = cos(radian);	mat[2][1] = sin(radian);	mat[3][1] = 0;
	mat[0][2] = 0;	mat[1][2] = -sin(radian);	mat[2][2] = cos(radian);	mat[3][2] = 0;
	mat[0][3] = 0;	mat[1][3] = 0;				mat[2][3] = 0;				mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRotYMat(float degree, GzMatrix mat)
{
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value

	float radian = (degree / 180) * PI;
	mat[0][0] = cos(radian);	mat[1][0] = 0;	mat[2][0] = -sin(radian);	mat[3][0] = 0;
	mat[0][1] = 0;				mat[1][1] = 1;	mat[2][1] = 0;				mat[3][1] = 0;
	mat[0][2] = sin(radian);	mat[1][2] = 0;	mat[2][2] = cos(radian);	mat[3][2] = 0;
	mat[0][3] = 0;				mat[1][3] = 0;	mat[2][3] = 0;				mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRotZMat(float degree, GzMatrix mat)
{
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value

	float radian = (degree / 180) * PI;
	mat[0][0] = cos(radian);	mat[1][0] = sin(radian);	mat[2][0] = 0;	mat[3][0] = 0;
	mat[0][1] = -sin(radian);	mat[1][1] = cos(radian);	mat[2][1] = 0;	mat[3][1] = 0;
	mat[0][2] = 0;				mat[1][2] = 0;				mat[2][2] = 1;	mat[3][2] = 0;
	mat[0][3] = 0;				mat[1][3] = 0;				mat[2][3] = 0;	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzTrxMat(GzCoord translate, GzMatrix mat)
{
	// Create translation matrix
	// Pass back the matrix using mat value
	mat[0][0] = 1;				mat[1][0] = 0;				mat[2][0] = 0;				mat[3][0] = 0;
	mat[0][1] = 0;				mat[1][1] = 1;				mat[2][1] = 0;				mat[3][1] = 0;
	mat[0][2] = 0;				mat[1][2] = 0;				mat[2][2] = 1;				mat[3][2] = 0;
	mat[0][3] = translate[0];	mat[1][3] = translate[1];	mat[2][3] = translate[2];	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzScaleMat(GzCoord scale, GzMatrix mat)
{
	// Create scaling matrix
	// Pass back the matrix using mat value

	mat[0][0] = scale[0];	mat[1][0] = 0;			mat[2][0] = 0;			mat[3][0] = 0;
	mat[0][1] = 0;			mat[1][1] = scale[1];	mat[2][1] = 0;			mat[3][1] = 0;
	mat[0][2] = 0;			mat[1][2] = 0;			mat[2][2] = scale[2];	mat[3][2] = 0;
	mat[0][3] = 0;			mat[1][3] = 0;			mat[2][3] = 0;			mat[3][3] = 1;

	return GZ_SUCCESS;
}

//----------------------------------------------------------
// Begin main functions

int GzPutCamera(GzRender *render, GzCamera *camera) {
	/*
	- overwrite renderer camera structure with new camera definition
	*/
	render->camera.FOV = camera->FOV;
	for (int i = 0; i < 3; i++) {
		render->camera.position[i] = camera->position[i];
		render->camera.lookat[i] = camera->lookat[i];
		render->camera.worldup[i] = camera->worldup[i];
	}

	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix) {
	/*
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	if (render->matlevel == 0) {
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = matrix[i][j];
				if(i == j)
					render->Xnorm[render->matlevel][i][j] = 1;
				else
					render->Xnorm[render->matlevel][i][j] = 0;
			}
	}
	else {
		GzMatrixMul(render->Ximage[render->matlevel - 1], matrix, render->Ximage[render->matlevel]);
		if (render->matlevel == 1)
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++){
					if (i == j)
						render->Xnorm[render->matlevel][i][j] = 1;
					else
						render->Xnorm[render->matlevel][i][j] = 0;
				}
		else{
			GzMatrix temp_M;
			float temp = 1 / sqrt(matrix[0][0] * matrix[0][0] + matrix[1][0] * matrix[1][0] + matrix[2][0] * matrix[2][0]);
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++) {
					temp_M[i][j] = matrix[i][j] * temp;
				}
				temp_M[i][3] = 0;
				temp_M[3][i] = 0;
			}	
			temp_M[3][3] = 1;
			GzMatrixMul(render->Xnorm[render->matlevel - 1], temp_M, render->Xnorm[render->matlevel]);
		}
	}
	render->matlevel++;

	return GZ_SUCCESS;
}

void GzMatrixMul(GzMatrix matrix1, GzMatrix matrix2, GzMatrix & matrix_out) {

	// for 4x4 GzMatrix Only
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			matrix_out[i][j] = 0;
			for (int k = 0; k < 4; k++)
				matrix_out[i][j] += matrix1[i][k] * matrix2[k][j];
		}
}

int GzPopMatrix(GzRender *render)
{
	/*
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (render->matlevel <= 0)
	{
		return GZ_FAILURE;
	}
	render->matlevel--;

	return GZ_SUCCESS;
}

float Dot(Vertex& v1, Vertex& v2){
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vertex Cross(Vertex& v1, Vertex& v2){
	Vertex v;
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}

Vertex Normal(Vertex v1) {
	Vertex v;
	float sqrt_sum = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
	v.x = v1.x / sqrt_sum;
	v.y = v1.y / sqrt_sum;
	v.z = v1.z / sqrt_sum;
	return v;
}

Vertex Normal_N(Vertex v1) {
	Vertex v;
	float sqrt_sum = sqrt(v1.nx * v1.nx + v1.ny * v1.ny + v1.nz * v1.nz);
	v.nx = v1.nx / sqrt_sum;
	v.ny = v1.ny / sqrt_sum;
	v.nz = v1.nz / sqrt_sum;
	return v;
}

void Normal_C(GzCoord & c1) {

	float sqrt_sum = sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
	for (int i = 0; i < 3; i++)
		c1[i] /= sqrt_sum;
}

int CreateXsp(GzRender *render)
{
	float radian = (render->camera.FOV / 180) * PI;
	float _d = tan(radian / 2);	// 1/d

	render->Xsp[0][0] = render->display->xres / 2;	render->Xsp[1][0] = 0;							render->Xsp[2][0] = 0;				render->Xsp[3][0] = 0;
	render->Xsp[0][1] = 0;							render->Xsp[1][1] = -render->display->yres / 2;	render->Xsp[2][1] = 0;				render->Xsp[3][1] = 0;
	render->Xsp[0][2] = 0;							render->Xsp[1][2] = 0;							render->Xsp[2][2] = INT_MAX * _d;	render->Xsp[3][2] = 0;
	render->Xsp[0][3] = render->display->xres / 2;	render->Xsp[1][3] = render->display->yres / 2;	render->Xsp[2][3] = 0;				render->Xsp[3][3] = 1;

	return GZ_SUCCESS;
}

int CreateXpi(GzRender *render)
{
	float radian = (render->camera.FOV / 180) * PI;
	float _d = tan(radian / 2); // 1/d

	render->camera.Xpi[0][0] = 1; render->camera.Xpi[1][0] = 0; render->camera.Xpi[2][0] = 0; render->camera.Xpi[3][0] = 0;
	render->camera.Xpi[0][1] = 0; render->camera.Xpi[1][1] = 1; render->camera.Xpi[2][1] = 0; render->camera.Xpi[3][1] = 0;
	render->camera.Xpi[0][2] = 0; render->camera.Xpi[1][2] = 0; render->camera.Xpi[2][2] = 1; render->camera.Xpi[3][2] = _d;
	render->camera.Xpi[0][3] = 0; render->camera.Xpi[1][3] = 0; render->camera.Xpi[2][3] = 0; render->camera.Xpi[3][3] = 1;

	return GZ_SUCCESS;
}

int CreateXiw(GzRender *render)
{
	Vertex cx;	Vertex cy;	Vertex cz;
	Vertex up = Vertex(render->camera.worldup);
	Vertex camera = Vertex(render->camera.position);

	cz = Normal(VertexMerge(1, render->camera.lookat, -1, render->camera.position));
	cy = Normal(VertexMerge(1, up, -Dot(up, cz), cz));
	cx = Cross(cy, cz);

	render->camera.Xiw[0][0] = cx.x;				render->camera.Xiw[1][0] = cy.x;				render->camera.Xiw[2][0] = cz.x;				render->camera.Xiw[3][0] = 0;
	render->camera.Xiw[0][1] = cx.y;				render->camera.Xiw[1][1] = cy.y;				render->camera.Xiw[2][1] = cz.y;				render->camera.Xiw[3][1] = 0;
	render->camera.Xiw[0][2] = cx.z;				render->camera.Xiw[1][2] = cy.z;				render->camera.Xiw[2][2] = cz.z;				render->camera.Xiw[3][2] = 0;
	render->camera.Xiw[0][3] = -Dot(cx, camera);	render->camera.Xiw[1][3] = -Dot(cy, camera);	render->camera.Xiw[2][3] = -Dot(cz, camera);	render->camera.Xiw[3][3] = 1;

	return GZ_SUCCESS;
}

int ScreenCheck(GzDisplay *display, vector<Vertex>& vert) {
	int counter = 0;
	for (int i = 0; i < 3; i++) {
		if (vert[i].z < 0)
			return 1;
		else if (vert[i].x < 0 || vert[i].x > display->xres || vert[i].y < 0 || vert[i].y > display->yres) {
			counter++;
			if (counter >= 3)
				return 1;
		}
	}
	return 0;
}

Vertex VertexMerge(float a, Vertex& v1, float b, Vertex& v2) {
	Vertex vert;
	vert.x = a * v1.x + b * v2.x;
	vert.y = a * v1.y + b * v2.y;
	vert.z = a * v1.z + b * v2.z;
	return vert;
}

Vertex VertexMerge(float a, GzCoord& c1, float b, GzCoord& c2) {
	Vertex vert;
	vert.x = a * c1[0] + b * c2[0];
	vert.y = a * c1[1] + b * c2[1];
	vert.z = a * c1[2] + b * c2[2];
	return vert;
}

int Color(GzRender *render, const Vertex v, GzColor color) {

	float nL, nE, rE;
	Vertex nVert, eye, reflect;
	Vertex lDir, lColor;
	Vertex sumS, sumD;

	eye.x = 0;		eye.y = 0;		eye.z = -1;
	nVert.x = v.nx;		nVert.y = v.ny;		nVert.z = v.nz;
	sumS.x = 0;		sumS.y = 0;		sumS.z = 0;
	sumD.x = 0;		sumD.y = 0;		sumD.z = 0;

	for (int i = 0; i < render->numlights; i++) {
		lDir.x = render->lights[i].direction[0];
		lDir.y = render->lights[i].direction[1];
		lDir.z = render->lights[i].direction[2];

		lColor.x = render->lights[i].color[0];
		lColor.y = render->lights[i].color[1];
		lColor.z = render->lights[i].color[2];

		nL = Dot(lDir, nVert);
		nE = Dot(eye, nVert);

		if (nL * nE > 0 && nL < 0) {
			nVert.x *= -1;
			nVert.y *= -1;
			nVert.z *= -1;
			nL *= -1;
		}

		else if (nL * nE < 0)
			continue;

		sumD.x += lColor.x * nL;
		sumD.y += lColor.y * nL;
		sumD.z += lColor.z * nL;

		reflect.x = 2 * nL * nVert.x - lDir.x;
		reflect.y = 2 * nL * nVert.y - lDir.y;
		reflect.z = 2 * nL * nVert.z - lDir.z;
		rE = Dot(Normal(reflect), eye);
		if (rE < 0)
			rE = 0;
		else if (rE > 1)
			rE = 1;

		rE = pow(rE, render->spec);
		sumS.x += lColor.x * rE;
		sumS.y += lColor.y * rE;
		sumS.z += lColor.z * rE;
		
	}

	color[0] = render->Ks[0] * sumS.x + render->Kd[0] * sumD.x + render->Ka[0] * render->ambientlight.color[0];
	color[1] = render->Ks[1] * sumS.y + render->Kd[1] * sumD.y + render->Ka[1] * render->ambientlight.color[1];
	color[2] = render->Ks[2] * sumS.z + render->Kd[2] * sumD.z + render->Ka[2] * render->ambientlight.color[2];

	for (int j = 0; j < 3; j++) {
		if (color[j] < 0)
			color[j] = 0;
		else if (color[j] > 1)
			color[j] = 1;
	}

	return GZ_SUCCESS;
}

void v_perspective(vector<Vertex> &v) {
	for (int i = 0; i < 3; i++)
	{
		v[i].tx = v[i].tx / (v[i].z / (INT_MAX - v[i].z) + 1);
		v[i].ty = v[i].ty / (v[i].z / (INT_MAX - v[i].z) + 1);
	}
}


void uv_perspective(Vertex &v, GzTextureIndex tex) {
	tex[0] = v.tx * (v.z / (INT_MAX - v.z) + 1);
	tex[1] = v.ty * (v.z / (INT_MAX - v.z) + 1);
}

void shiftVerts(GzRender *render, vector<Vertex>& v) {
	for (int i = 0; i < 3; i++) {
		v[i].x -= render->shift_offsetX;
		v[i].y -= render->shift_offsetY;
	}
}