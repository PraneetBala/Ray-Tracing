/************************************************************

Name : Praneet Bala
Student ID: 5393094
Date: 11/08/2019
Code: Assignment 1d
Objective: To implement reflection and refraction in solid objects

************************************************************/

// The total iterations for this code is 10; can be updated

// Instructions to run the code: 
// > make
// > ./out test#.txt

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <stack>
#define PI 3.14159265
#define epsilon 0.00009
#define iterations 10	// Set the number of iterations for the rays to bounce
using namespace std;


struct size1_t {
	float x;
} vfov;

struct size2_t {
	float w,h;
} imsize;

struct size3_t {
	float x,y,z;
} eye, viewdir, updir, bkgcolor, mtlcolor;

struct sphere_t {
	float x,y,z;
	float r;
	int tex;
	int m,t;
};

struct RayType {
	float x,y,z;
	float dx,dy,dz;
} ray;

struct ColorType {
	float odr,odg,odb;
	float osr,osg,osb;
	float ka,kd,ks,n;
	float alpha,eta;
};

struct LightType {
	float x,y,z;
	float w;
	float r,g,b;
};

struct FaceType {
	int x,y,z;
	int nx,ny,nz;
	int tx,ty,tz;
	int n,tex; // Flags for smooth shading and texture
	int m,t; // Array value of material color and texture
};

struct TexType {
	int w,h;
	vector<float> r;
	vector<float> g;
	vector<float> b;
	int t;
};

// Function used in determining the output image name
void FindAllOccurences(std::vector<size_t> & vec, std::string data, std::string toSearch) {
	size_t pos = data.find(toSearch);
	while( pos != std::string::npos) {
		vec.push_back(pos);
		pos = data.find(toSearch, pos + toSearch.size());
	}
}

// Function used in extracting information from the scene description file
void AssignAttributes(std::vector<size_t> & vec, std::string data) {
	string attrib;
	attrib = data.substr(0, vec[0]);

	if (attrib == "vfov") {
		vfov.x = stof(data.substr(vec[0]+1, data.length() - vec[0] - 1));
	}
	else if (attrib == "imsize") {
		imsize.w = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		imsize.h = stof(data.substr(vec[1]+1, data.length() - vec[1] - 1));
	}
	else if (attrib == "eye") {
		eye.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		eye.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
		eye.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	}
	else if (attrib == "viewdir") {
		viewdir.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		viewdir.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
		viewdir.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	}
	else if (attrib == "updir") {
		updir.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		updir.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
		updir.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	}
	else if (attrib == "bkgcolor") {
		bkgcolor.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		bkgcolor.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
		bkgcolor.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	}
}

// Function to extract the sphere information
sphere_t AssignSphere(std::vector<size_t> & vec, std::string data, int m, int t) {
	sphere_t sphere;
	sphere.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	sphere.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	sphere.z = stof(data.substr(vec[2]+1, vec[3] - vec[2] - 1));
	sphere.r = stof(data.substr(vec[3]+1, data.length() - vec[3] - 1));
	sphere.tex = 0;
	sphere.m = m;
	sphere.t = t;
	if (t != -1) 
		sphere.tex = 1;
	return sphere;
}

// Function to extract the color information
ColorType AssignMtlColor(std::vector<size_t> & vec, std::string data) {
	ColorType mtlC;
	mtlC.odr = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	mtlC.odg = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	mtlC.odb = stof(data.substr(vec[2]+1, vec[3] - vec[2] - 1));
	mtlC.osr = stof(data.substr(vec[3]+1, vec[4] - vec[3] - 1));
	mtlC.osg = stof(data.substr(vec[4]+1, vec[5] - vec[4] - 1));
	mtlC.osb = stof(data.substr(vec[5]+1, vec[6] - vec[5] - 1));
	mtlC.ka = stof(data.substr(vec[6]+1, vec[7] - vec[6] - 1));
	mtlC.kd = stof(data.substr(vec[7]+1, vec[8] - vec[7] - 1));
	mtlC.ks = stof(data.substr(vec[8]+1, vec[9] - vec[8] - 1));
	mtlC.n = stof(data.substr(vec[9]+1, vec[10] - vec[9] - 1));
	mtlC.alpha = stof(data.substr(vec[10]+1, vec[11] - vec[10] - 1));
	mtlC.eta = stof(data.substr(vec[11]+1, data.length() - vec[11] - 1));
	return mtlC;
}

// Function to extract the light information
LightType AssignLight(std::vector<size_t> & vec, std::string data) {
	LightType lt;
	lt.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	lt.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	lt.z = stof(data.substr(vec[2]+1, vec[3] - vec[2] - 1));
	lt.w = stof(data.substr(vec[3]+1, vec[4] - vec[3] - 1));
	lt.r = stof(data.substr(vec[4]+1, vec[5] - vec[4] - 1));
	lt.g = stof(data.substr(vec[5]+1, vec[6] - vec[5] - 1));
	lt.b = stof(data.substr(vec[6]+1, data.length() - vec[6] - 1));
	return lt;
}

// Function to extract vertex points
size3_t AssignPoints(std::vector<size_t> & vec, std::string data) {
	size3_t Pt;
	Pt.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	Pt.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	Pt.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	return Pt;
}

// Function to extract vertex normals
size3_t AssignNormals(std::vector<size_t> & vec, std::string data) {
	size3_t Pt;
	Pt.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	Pt.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	Pt.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	return Pt;
}

// Function to assign textures to the vertices of a face
size2_t AssignVtex(std::vector<size_t> & vec, std::string data) {
	size2_t Pt;
	Pt.w = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	Pt.h = stof(data.substr(vec[1]+1, data.length() - vec[1] - 1));
	return Pt;
}

// Function to save the texture information
TexType AssignTextures(std::vector<size_t> & vec, std::string data) {
	TexType tex;
	string str = data.substr(vec[0]+1, data.length() - vec[0] - 1);

	ifstream inTex;
	string attrib;
	inTex.open(str);
	if (!inTex) {
		cerr << "Unable to open ppm file\n";
		exit(1);
	}
	getline(inTex,str);
	if (inTex && !str.empty()) {
		std::vector<size_t> vc;
		FindAllOccurences(vc, str, " ");
		attrib = str.substr(0, vc[0]);
		tex.w = stoi(str.substr(vc[0]+1, vc[1] - vc[0] - 1));
		tex.h = stoi(str.substr(vc[1]+1, vc[2] - vc[1] - 1));
	}
	while (inTex) {
		getline(inTex,str);
		if (inTex && !str.empty())
			tex.r.push_back(stoi(str));
		getline(inTex,str);
		if (inTex && !str.empty())
			tex.g.push_back(stoi(str));
		getline(inTex,str);
		if (inTex && !str.empty())
			tex.b.push_back(stoi(str));
	}
	inTex.close();
	return tex;
}

// Function to handle the different cases of face input
FaceType AssignFaces(std::vector<size_t> & vec, std::string data, int m, int t) {
	FaceType Fc;
	vector<size_t> v1,v2,v3;
	string s;

	s = data.substr(vec[0]+1, vec[1] - vec[0] - 1);
	FindAllOccurences(v1, s, "/");
	if (v1.size() == 0) {	// When only vertex points are given
		Fc.n = 0;
		Fc.tex = 0;
		Fc.x = stoi(data.substr(vec[0]+1, vec[1] - vec[0] - 1)) - 1;
	}
	else if (v1.size() == 1) {	// When points and texture information is given
		Fc.n = 0;
		Fc.tex = 1;
		Fc.x = stoi(s.substr(0, v1[0])) - 1;
		Fc.tx = stoi(s.substr(v1[0]+1, s.length() - v1[0] - 1)) - 1;
	}
	else if (v1.size() == 2) {
		Fc.n = 1;
		if (v1[1] == v1[0]+1) {	// When points and normal information is given
			Fc.tex = 0;
			Fc.x = stoi(s.substr(0, v1[0])) - 1;
			Fc.nx = stoi(s.substr(v1[1]+1, s.length() - v1[1] - 1)) - 1;
		}
		else if (v1[1] != v1[0]+1) {	// When points, texture and normal information is given
			Fc.tex = 1;
			Fc.x = stoi(s.substr(0, v1[0])) - 1;
			Fc.tx = stoi(s.substr(v1[0]+1, v1[1] - v1[0] - 1)) - 1;
			Fc.nx = stoi(s.substr(v1[1]+1, s.length() - v1[1] - 1)) - 1;
		}
	}

	s = data.substr(vec[1]+1, vec[2] - vec[1] - 1);
	FindAllOccurences(v2, s, "/");
	if (v2.size() == 0) {
		Fc.y = stoi(data.substr(vec[1]+1, vec[2] - vec[1] - 1)) - 1;
	}
	else if (v2.size() == 1) {
		Fc.y = stoi(s.substr(0, v2[0])) - 1;
		Fc.ty = stoi(s.substr(v2[0]+1, s.length() - v2[0] - 1)) - 1;
	}
	else if (v2.size() == 2) {
		if (v2[1] == v2[0]+1) {
			Fc.y = stoi(s.substr(0, v2[0])) - 1;
			Fc.ny = stoi(s.substr(v2[1]+1, s.length() - v2[1] - 1)) - 1;
		}
		else if (v2[1] != v2[0]+1) {
			Fc.y = stoi(s.substr(0, v2[0])) - 1;
			Fc.ty = stoi(s.substr(v2[0]+1, v2[1] - v2[0] - 1)) - 1;
			Fc.ny = stoi(s.substr(v2[1]+1, s.length() - v2[1] - 1)) - 1;
		}
	}

	s = data.substr(vec[2]+1, data.length() - vec[2] - 1);
	FindAllOccurences(v3, s, "/");
	if (v3.size() == 0) {
		Fc.z = stoi(data.substr(vec[2]+1, data.length() - vec[2] - 1)) - 1;
	}
	else if (v3.size() == 1) {
		Fc.z = stoi(s.substr(0, v3[0])) - 1;
		Fc.tz = stoi(s.substr(v3[0]+1, s.length() - v3[0] - 1)) - 1;
	}
	else if (v3.size() == 2) {
		if (v3[1] == v3[0]+1) {
			Fc.z = stoi(s.substr(0, v3[0])) - 1;
			Fc.nz = stoi(s.substr(v3[1]+1, s.length() - v3[1] - 1)) - 1;
		}
		else if (v3[1] != v3[0]+1) {
			Fc.z = stoi(s.substr(0, v3[0])) - 1;
			Fc.tz = stoi(s.substr(v3[0]+1, v3[1] - v3[0] - 1)) - 1;
			Fc.nz = stoi(s.substr(v3[1]+1, s.length() - v3[1] - 1)) - 1;
		}
	}

	Fc.m = m;
	Fc.t = t;
	return Fc;
}

// Function that calculates the cross product between vectors A and B 
// ( i.e. cross_P = AxB )
void Cross(float A[], float B[], float cross_P[]) {
	cross_P[0] = A[1]*B[2] - A[2]*B[1];
	cross_P[1] = A[2]*B[0] - A[0]*B[2];
	cross_P[2] = A[0]*B[1] - A[1]*B[0];
}

// Function used for calculating dot product between vectors A and B
float Dot(float A[], float B[]) {
	float product = 0;
	for (int i = 0; i < 3; i++) {
		product = product + A[i]*B[i];
	}
	return product;
}

// Function used to normalize a given vector A
void Normalize(float A[], float N[]) {
	float norm = sqrt(Dot(A,A));
	for (int i = 0; i < 3; i++) {
		N[i] = A[i]/norm;
	}
}

// Function to assign color to background and foreground pixels
ColorType Shade_Ray(int tok, vector<ColorType> & mc, int tu, float N[], vector<LightType> & light, float inter[], vector<sphere_t> & sphere, int pa, vector<size3_t> & v_points, vector<FaceType> & face, int ij[], vector<TexType> & tex) {
	ColorType color;
	float z = 0;
	float NL, NH;
	float L[3], H[3];
	float view_d[] = {viewdir.x, viewdir.y, viewdir.z};	// Viewing direction
	float view[3];
	Normalize(view_d, view);
	float dist, t;
	float S = 1;
	float inter2[] = {0,0,0};

	// Token value for background color data
	if (tok == 0) {
		color.odr = bkgcolor.x;
		color.odg = bkgcolor.y;
		color.odb = bkgcolor.z;
	}
	// Token value for sphere color assignment
	else if (tok == 1) {
		if (sphere[pa].tex == 1) {	// Using texture data
			color.odr = tex[sphere[pa].t].r[ij[0]+tex[sphere[pa].t].w*ij[1]]/255;
			color.odg = tex[sphere[pa].t].g[ij[0]+tex[sphere[pa].t].w*ij[1]]/255;
			color.odb = tex[sphere[pa].t].b[ij[0]+tex[sphere[pa].t].w*ij[1]]/255;
		}
		else {
			color.odr = mc[tu].odr;
			color.odg = mc[tu].odg;
			color.odb = mc[tu].odb;
		}
		for (int i = 0; i < light.size(); i++) {
			S = 1;
			if (light[i].w == 0) {
				float L_n[] = {-light[i].x, -light[i].y, -light[i].z};
				dist = 100000.0;
				Normalize(L_n, L);
				float H_n[] = {L[0]-view[0], L[1]-view[1], L[2]-view[2]};
				Normalize(H_n, H);
				NL = Dot(N, L);
				NH = Dot(N, H);
			}
			else if (light[i].w == 1) {
				float L_n[] = {light[i].x-inter[0], light[i].y-inter[1], light[i].z-inter[2]};
				dist = sqrt(Dot(L_n,L_n));
				Normalize(L_n, L);
				float H_n[] = {L[0]-view[0], L[1]-view[1], L[2]-view[2]};
				Normalize(H_n, H);
				NL = Dot(N, L);
				NH = Dot(N, H);
			}

			for (int j = 0; j < sphere.size(); j++) {
				if (j == pa)
					continue;
				float Bi, Ci, t1, t2, discri;
				Bi = 2 * (L[0]*(inter[0]-sphere[j].x) + L[1]*(inter[1]-sphere[j].y) + L[2]*(inter[2]-sphere[j].z));
				Ci = pow((inter[0]-sphere[j].x),2) + pow((inter[1]-sphere[j].y),2) + pow((inter[2]-sphere[j].z),2) - pow(sphere[j].r,2);
				discri = Bi*Bi - 4*Ci;
				
				if (discri > 0) {
					t1 = (-Bi + sqrt(discri)) / 2;
					t2 = (-Bi - sqrt(discri)) / 2;
					if (t1 > 0 && t2 > 0)
						t = min(t1,t2);
					else if (t1 > 0 && t2 < 0)
						t = t1;
					else if (t2 > 0 && t1 < 0)
						t = t2;
					else
						t = 100000;
				}
				else if (discri == 0) {
					t = -Bi/2;
				}
				else if (discri < 0) {
					t = 100000;
				}
				
				if (t > 0) {
					if (t < dist) {
						S = S*(1-mc[sphere[j].m].alpha);					
					}
				}
			}
			float p1[3], p2[3], p3[3], e1[3], e2[3], n[3], d;
			for (int j = 0; j < face.size(); j++) {
				t = -1;
				p1[0] = v_points[face[j].x].x;
				p1[1] = v_points[face[j].x].y;
				p1[2] = v_points[face[j].x].z;

				p2[0] = v_points[face[j].y].x;
				p2[1] = v_points[face[j].y].y;
				p2[2] = v_points[face[j].y].z;

				p3[0] = v_points[face[j].z].x;
				p3[1] = v_points[face[j].z].y;
				p3[2] = v_points[face[j].z].z;

				for (int q = 0; q < 3; q++) {
					e1[q] = p2[q] - p1[q];
					e2[q] = p3[q] - p1[q]; 
				}
				Cross(e1, e2, n);
				d = -(n[0]*p1[0] + n[1]*p1[1] + n[2]*p1[2]);

				float num, den;
				num = n[0]*inter[0] + n[1]*inter[1] + n[2]*inter[2] + d;
				den = n[0]*L[0] + n[1]*L[1] + n[2]*L[2];
				if (den != 0) {
					t = -num/den;
				}
				if (t > 0 && t < dist) {
					inter2[0] = inter[0]+(t*L[0]); 
					inter2[1] = inter[1]+(t*L[1]);
					inter2[2] = inter[2]+(t*L[2]);
					Cross(e1,e2,n);
					float A1 = 0.5*sqrt(Dot(n,n));
					float e3[3], e4[3];
					for (int q = 0; q < 3; q++) {
						e3[q] = inter2[q] - p2[q];
						e4[q] = inter2[q] - p3[q]; 
					}
					Cross(e3,e4,n);
					float a1 = 0.5*sqrt(Dot(n,n));
					Cross(e4,e2,n);
					float b1 = 0.5*sqrt(Dot(n,n));
					Cross(e1,e3,n);
					float c1 = 0.5*sqrt(Dot(n,n));
					float alpha = a1/A1;
					float beta = b1/A1;
					float gamma = c1/A1;
					if (alpha+beta+gamma-1 < epsilon) {
						S = S*(1-mc[face[j].m].alpha);
					}
				}
			}
			color.odr = mc[tu].ka * color.odr + S * light[i].r * (mc[tu].kd * color.odr * max(z,NL) + mc[tu].ks * mc[tu].osr * max(z,pow(NH,mc[tu].n)));
			color.odg = mc[tu].ka * color.odg + S * light[i].g * (mc[tu].kd * color.odg * max(z,NL) + mc[tu].ks * mc[tu].osg * max(z,pow(NH,mc[tu].n)));
			color.odb = mc[tu].ka * color.odb + S * light[i].b * (mc[tu].kd * color.odb * max(z,NL) + mc[tu].ks * mc[tu].osb * max(z,pow(NH,mc[tu].n)));
		}
	}
	// Token value for face color assignment
	else if (tok == 2) {
		if (face[pa].tex == 1) {	// Using texture data
			color.odr = tex[face[pa].t].r[ij[0]+tex[face[pa].t].w*ij[1]]/255;
			color.odg = tex[face[pa].t].g[ij[0]+tex[face[pa].t].w*ij[1]]/255;
			color.odb = tex[face[pa].t].b[ij[0]+tex[face[pa].t].w*ij[1]]/255;
		}
		else {
			color.odr = mc[tu].odr;
			color.odg = mc[tu].odg;
			color.odb = mc[tu].odb;
		}
		for (int i = 0; i < light.size(); i++) {
			S = 1;
			if (light[i].w == 0) {
				float L_n[] = {-light[i].x, -light[i].y, -light[i].z};
				dist = 100000.0;
				Normalize(L_n, L);
				float H_n[] = {L[0]-view[0], L[1]-view[1], L[2]-view[2]};
				Normalize(H_n, H);
				NL = Dot(N, L);
				NH = Dot(N, H);
			}
			else if (light[i].w == 1) {
				float L_n[] = {light[i].x-inter[0], light[i].y-inter[1], light[i].z-inter[2]};
				dist = sqrt(Dot(L_n,L_n));
				Normalize(L_n, L);
				float H_n[] = {L[0]-view[0], L[1]-view[1], L[2]-view[2]};
				Normalize(H_n, H);
				NL = Dot(N, L);
				NH = Dot(N, H);
			}

			for (int j = 0; j < sphere.size(); j++) {
				float Bi, Ci, t1, t2, discri;
				Bi = 2 * (L[0]*(inter[0]-sphere[j].x) + L[1]*(inter[1]-sphere[j].y) + L[2]*(inter[2]-sphere[j].z));
				Ci = pow((inter[0]-sphere[j].x),2) + pow((inter[1]-sphere[j].y),2) + pow((inter[2]-sphere[j].z),2) - pow(sphere[j].r,2);
				discri = Bi*Bi - 4*Ci;
				
				if (discri > 0) {
					t1 = (-Bi + sqrt(discri)) / 2;
					t2 = (-Bi - sqrt(discri)) / 2;
					if (t1 > 0 && t2 > 0)
						t = min(t1,t2);
					else if (t1 > 0 && t2 < 0)
						t = t1;
					else if (t2 > 0 && t1 < 0)
						t = t2;
					else
						t = 100000;
				}
				else if (discri == 0) {
					t = -Bi/2;
				}
				else if (discri < 0) {
					t = 100000;
				}
				
				if (t > 0) {
					if (t < dist) {
						S = S*(1-mc[sphere[j].m].alpha);
					}
				}
			}
			float p1[3], p2[3], p3[3], e1[3], e2[3], n[3], d;
			for (int j = 0; j < face.size(); j++) {
				if (j == pa)
					continue;
				t = -1;
				p1[0] = v_points[face[j].x].x;
				p1[1] = v_points[face[j].x].y;
				p1[2] = v_points[face[j].x].z;

				p2[0] = v_points[face[j].y].x;
				p2[1] = v_points[face[j].y].y;
				p2[2] = v_points[face[j].y].z;

				p3[0] = v_points[face[j].z].x;
				p3[1] = v_points[face[j].z].y;
				p3[2] = v_points[face[j].z].z;

				for (int q = 0; q < 3; q++) {
					e1[q] = p2[q] - p1[q];
					e2[q] = p3[q] - p1[q]; 
				}
				Cross(e1, e2, n);
				d = -(n[0]*p1[0] + n[1]*p1[1] + n[2]*p1[2]);

				float num, den;
				num = n[0]*inter[0] + n[1]*inter[1] + n[2]*inter[2] + d;
				den = n[0]*L[0] + n[1]*L[1] + n[2]*L[2];
				if (den != 0) {
					t = -num/den;
				}
				if (t > 0 && t < dist) {
					inter2[0] = inter[0]+(t*L[0]); 
					inter2[1] = inter[1]+(t*L[1]);
					inter2[2] = inter[2]+(t*L[2]);
					Cross(e1,e2,n);
					float A1 = 0.5*sqrt(Dot(n,n));
					float e3[3], e4[3];
					for (int q = 0; q < 3; q++) {
						e3[q] = inter2[q] - p2[q];
						e4[q] = inter2[q] - p3[q]; 
					}
					Cross(e3,e4,n);
					float a1 = 0.5*sqrt(Dot(n,n));
					Cross(e4,e2,n);
					float b1 = 0.5*sqrt(Dot(n,n));
					Cross(e1,e3,n);
					float c1 = 0.5*sqrt(Dot(n,n));
					float alpha = a1/A1;
					float beta = b1/A1;
					float gamma = c1/A1;
					if (alpha+beta+gamma-1 < epsilon) {
						S = S*(1-mc[face[j].m].alpha);
					}
				}
			}
			color.odr = mc[tu].ka * color.odr + S * light[i].r * (mc[tu].kd * color.odr * max(z,NL) + mc[tu].ks * mc[tu].osr * max(z,pow(NH,mc[tu].n)));
			color.odg = mc[tu].ka * color.odg + S * light[i].g * (mc[tu].kd * color.odg * max(z,NL) + mc[tu].ks * mc[tu].osg * max(z,pow(NH,mc[tu].n)));
			color.odb = mc[tu].ka * color.odb + S * light[i].b * (mc[tu].kd * color.odb * max(z,NL) + mc[tu].ks * mc[tu].osb * max(z,pow(NH,mc[tu].n)));
		}
	}
	return color;
}

// Function to identify the intersection with sphere
ColorType Trace_Ray(RayType ray, vector<sphere_t> & sphere, vector<ColorType> & mtl, vector<LightType> & light, vector<size3_t> & v_points, vector<FaceType> & face, vector<size3_t> & v_normals, vector<size2_t> & v_tex, vector<TexType> & tex, vector <float> & eta123, int depth) {
	int pa, tu, pa1, pa2, token;
	float min_T = 100000;
	float t;
	float N[] = {0,0,0};
	float N1[] = {0,0,0};
	float inter[] = {0,0,0};
	float inter1[] = {0,0,0};
	float inter2[] = {0,0,0};
	int ij[] = {0,0};
	
	for (int i = 0; i < sphere.size(); i++) {
		float Bi, Ci, t1, t2, discri;
		Bi = 2 * (ray.dx*(ray.x-sphere[i].x) + ray.dy*(ray.y-sphere[i].y) + ray.dz*(ray.z-sphere[i].z));
		Ci = pow((ray.x-sphere[i].x),2) + pow((ray.y-sphere[i].y),2) + pow((ray.z-sphere[i].z),2) - pow(sphere[i].r,2);
		discri = Bi*Bi - 4*Ci;
		if (discri > 0) {
			t1 = (-Bi + sqrt(discri)) / 2;
			t2 = (-Bi - sqrt(discri)) / 2;
			if (t1 > 0 && t2 > 0)
				t = min(t1,t2);
			else if (t1 > 0 && t2 < 0)
				t = t1;
			else if (t2 > 0 && t1 < 0)
				t = t2;
			else 
				t = min_T;
		}
		else if (discri == 0) {
			t = -Bi/2;
		}
		else if (discri < 0) {
			t = 100000;
		}
		
		if (t > 0) {
			if (min(min_T,t) != min_T) {
				pa1 = i;
			}
			min_T = min(min_T, t);
		}
	}
	if (min_T != 100000) {
		inter1[0] = ray.x+(min_T*ray.dx); 
		inter1[1] = ray.y+(min_T*ray.dy);
		inter1[2] = ray.z+(min_T*ray.dz); // Point of intersection between ray and sphere[pa]
		N1[0] = (inter1[0]-sphere[pa1].x)/sphere[pa1].r;
		N1[1] = (inter1[1]-sphere[pa1].y)/sphere[pa1].r;
		N1[2] = (inter1[2]-sphere[pa1].z)/sphere[pa1].r;
		inter1[0] = inter1[0] + 0.002*N1[0];
		inter1[1] = inter1[1] + 0.002*N1[1];
		inter1[2] = inter1[2] + 0.002*N1[2];
	}

	float min2 = min_T;

	float p1[3], p2[3], p3[3], e1[3], e2[3], n[3], d;
	for (int j = 0; j < face.size(); j++) {
		t = -1;
		p1[0] = v_points[face[j].x].x;
		p1[1] = v_points[face[j].x].y;
		p1[2] = v_points[face[j].x].z;

		p2[0] = v_points[face[j].y].x;
		p2[1] = v_points[face[j].y].y;
		p2[2] = v_points[face[j].y].z;

		p3[0] = v_points[face[j].z].x;
		p3[1] = v_points[face[j].z].y;
		p3[2] = v_points[face[j].z].z;

		for (int q = 0; q < 3; q++) {
			e1[q] = p2[q] - p1[q];
			e2[q] = p3[q] - p1[q]; 
		}

		Cross(e1, e2, n);
		d = -(n[0]*p1[0] + n[1]*p1[1] + n[2]*p1[2]);

		float num, den;
		num = n[0]*ray.x + n[1]*ray.y + n[2]*ray.z + d;
		den = n[0]*ray.dx + n[1]*ray.dy + n[2]*ray.dz;
		if (den != 0) {
			t = -num/den;
		}

		if (t > 0) {
			inter2[0] = ray.x+(t*ray.dx); 
			inter2[1] = ray.y+(t*ray.dy);
			inter2[2] = ray.z+(t*ray.dz);
			Cross(e1,e2,n);
			float A1 = 0.5*sqrt(Dot(n,n));
			float e3[3], e4[3];
			for (int q = 0; q < 3; q++) {
				e3[q] = inter2[q] - p2[q];
				e4[q] = inter2[q] - p3[q]; 
			}
			Cross(e3,e4,n);
			float a1 = 0.5*sqrt(Dot(n,n));
			Cross(e4,e2,n);
			float b1 = 0.5*sqrt(Dot(n,n));
			Cross(e1,e3,n);
			float c1 = 0.5*sqrt(Dot(n,n));
			float alpha = a1/A1;
			float beta = b1/A1;
			float gamma = c1/A1;
			if ((min(min_T,t) != min_T) && (alpha >= 0) && (alpha <= 1) && (beta >= 0) && (beta <= 1) && (gamma >= 0) && (gamma <= 1) && (alpha+beta+gamma-1 < epsilon)) {
				pa2 = j;
				min_T = min(min_T, t);
			}
		}
	}

	// Token differentiates between background and foreground pixels
	if (min_T == 100000) {
		token = 0; // Background pixel
		tu = 0;
	}
	else {
		if (min_T != min2) {	// When triangle is closer
			token = 2;
			tu = face[pa2].m;
			int j = pa2;
			inter[0] = ray.x+(min_T*ray.dx); 
			inter[1] = ray.y+(min_T*ray.dy);
			inter[2] = ray.z+(min_T*ray.dz);

			p1[0] = v_points[face[j].x].x;
			p1[1] = v_points[face[j].x].y;
			p1[2] = v_points[face[j].x].z;

			p2[0] = v_points[face[j].y].x;
			p2[1] = v_points[face[j].y].y;
			p2[2] = v_points[face[j].y].z;

			p3[0] = v_points[face[j].z].x;
			p3[1] = v_points[face[j].z].y;
			p3[2] = v_points[face[j].z].z;

			for (int q = 0; q < 3; q++) {
				e1[q] = p2[q] - p1[q];
				e2[q] = p3[q] - p1[q]; 
			}

			if (face[pa2].n == 0) {
				float uv[] = {0,0};
				float vn[3], vn1[3], vn2[3], vn3[3];
				Cross(e1, e2, n);
				float A1 = 0.5*sqrt(Dot(n,n));
				float e3[3], e4[3];
				for (int q = 0; q < 3; q++) {
					e3[q] = inter[q] - p2[q];
					e4[q] = inter[q] - p3[q]; 
				}
				Cross(e3,e4,n);
				float a1 = 0.5*sqrt(Dot(n,n));
				Cross(e4,e2,n);
				float b1 = 0.5*sqrt(Dot(n,n));
				Cross(e1,e3,n);
				float c1 = 0.5*sqrt(Dot(n,n));
				float alpha = a1/A1;
				float beta = b1/A1;
				float gamma = c1/A1;

				Normalize(n,N);

				if (face[pa2].tex == 1) {
					vn1[0] = v_tex[face[j].tx].w;
					vn1[1] = v_tex[face[j].tx].h;

					vn2[0] = v_tex[face[j].ty].w;
					vn2[1] = v_tex[face[j].ty].h;

					vn3[0] = v_tex[face[j].tz].w;
					vn3[1] = v_tex[face[j].tz].h;

					uv[0] = alpha*vn1[0] + beta*vn2[0] + gamma*vn3[0];
					uv[1] = alpha*vn1[1] + beta*vn2[1] + gamma*vn3[1];

					ij[0] = round(uv[0]*(tex[face[j].t].w-1));
					ij[1] = round(uv[1]*(tex[face[j].t].h-1));
				}
			}
			else if (face[pa2].n == 1) {	// In case of smooth shading
				float uv[] = {0,0};
				float vn[3], vn1[3], vn2[3], vn3[3];
				vn1[0] = v_normals[face[j].nx].x;
				vn1[1] = v_normals[face[j].nx].y;
				vn1[2] = v_normals[face[j].nx].z;
				Normalize(vn1,vn1);

				vn2[0] = v_normals[face[j].ny].x;
				vn2[1] = v_normals[face[j].ny].y;
				vn2[2] = v_normals[face[j].ny].z;
				Normalize(vn2,vn2);

				vn3[0] = v_normals[face[j].nz].x;
				vn3[1] = v_normals[face[j].nz].y;
				vn3[2] = v_normals[face[j].nz].z;
				Normalize(vn3,vn3);

				Cross(e1,e2,n);
				float A1 = 0.5*sqrt(Dot(n,n));
				float e3[3], e4[3];
				for (int q = 0; q < 3; q++) {
					e3[q] = inter[q] - p2[q];
					e4[q] = inter[q] - p3[q]; 
				}
				Cross(e3,e4,n);
				float a1 = 0.5*sqrt(Dot(n,n));
				Cross(e4,e2,n);
				float b1 = 0.5*sqrt(Dot(n,n));
				Cross(e1,e3,n);
				float c1 = 0.5*sqrt(Dot(n,n));
				float alpha = a1/A1;
				float beta = b1/A1;
				float gamma = c1/A1;

				vn[0] = alpha*vn1[0] + beta*vn2[0] + gamma*vn3[0];
				vn[1] = alpha*vn1[1] + beta*vn2[1] + gamma*vn3[1];
				vn[2] = alpha*vn1[2] + beta*vn2[2] + gamma*vn3[2];
				Normalize(vn,N);

				if (face[pa2].tex == 1) {
					vn1[0] = v_tex[face[j].tx].w;
					vn1[1] = v_tex[face[j].tx].h;

					vn2[0] = v_tex[face[j].ty].w;
					vn2[1] = v_tex[face[j].ty].h;

					vn3[0] = v_tex[face[j].tz].w;
					vn3[1] = v_tex[face[j].tz].h;

					uv[0] = alpha*vn1[0] + beta*vn2[0] + gamma*vn3[0];
					uv[1] = alpha*vn1[1] + beta*vn2[1] + gamma*vn3[1];

					ij[0] = round(uv[0]*(tex[face[j].t].w-1));
					ij[1] = round(uv[1]*(tex[face[j].t].h-1));
				}
			}
			pa = pa2;
		}
		else if (min_T == min2) {	// When sphere is closer
			token = 1;
			pa = pa1;
			tu = sphere[pa].m;
			inter[0] = inter1[0];
			inter[1] = inter1[1];
			inter[2] = inter1[2];
			float uv[] = {0,0};
			N[0] = N1[0];
			N[1] = N1[1];
			N[2] = N1[2];

			if (sphere[pa].tex == 1) {
				float phi,theta;
				phi = acos(-(inter[0]-sphere[pa].x)/sphere[pa].r);
				theta = atan2(-(inter[2]-sphere[pa].z),(inter[1]-sphere[pa].y));
				
				if (theta < 0)
					theta = theta + 2*PI;
				uv[0] = 0.5*theta/PI;
				uv[1] = phi/PI;
				ij[0] = round(uv[0]*(tex[sphere[pa].t].w-1));
				ij[1] = round(uv[1]*(tex[sphere[pa].t].h-1));
			}
		}
	}

	// The reflection and refraction calculation
	ColorType color;
	if (depth < iterations && token != 0) {
		float N1[3];
		float I[] = {ray.dx,ray.dy,ray.dz};
		Normalize(I,I);
		I[0] = -I[0];
		I[1] = -I[1];
		I[2] = -I[2];
		Normalize(N,N);
		float ni = Dot(N,I);
		float eta_i, eta_t;
		float theta_i, theta_c;
		RayType T, R;
		float epsilo = 10*epsilon;
		
		float Fr, Fo;
		N1[0] = N[0];
		N1[1] = N[1];
		N1[2] = N[2];
		int not_T = 0;

		/*
		cout << "\nToken: " << token << endl;
		cout << min_T << endl;
		cout << "Ray: " << ray.x << " " << ray.y << " " << ray.z << endl;
		cout << "Ray_dir: " << ray.dx << " " << ray.dy << " " << ray.dz << endl;
		cout << "Intersection: " << inter[0] << ", " << inter[1] << ", " << inter[2] << endl;
		cout << "Normal: " << N[0] << ", " << N[1] << ", " << N[2] << endl;
		cout << "Incident: " << I[0] << ", " << I[1] << ", " << I[2] << endl;
		cout << "Normal1: " << N1[0] << ", " << N1[1] << ", " << N1[2] << endl;
		cout << "N.I: " << ni << endl;
		*/

		// Entry condition
		if (ni >= 0) {
			//cout << "1\n";
			eta_i = eta123[eta123.size()-1];
			eta_t = mtl[tu].eta;
			eta123.push_back(eta_t);

			T.x = inter[0] - 10*epsilo*N1[0];
			T.y = inter[1] - 10*epsilo*N1[1];
			T.z = inter[2] - 10*epsilo*N1[2];

			R.x = inter[0] + epsilo*N1[0];
			R.y = inter[1] + epsilo*N1[1];
			R.z = inter[2] + epsilo*N1[2];

			R.dx = 2*ni*N1[0] - I[0];
			R.dy = 2*ni*N1[1] - I[1];
			R.dz = 2*ni*N1[2] - I[2];
			
			Fo = pow((eta_t-eta_i)/(eta_t+eta_i),2);
			Fr = Fo + (1-Fo)*pow((1-ni),5);
			
			if (eta_i <= eta_t) {
				//cout << "1.1\n";
				T.dx = -N1[0]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[0]-I[0]);
				T.dy = -N1[1]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[1]-I[1]);
				T.dz = -N1[2]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[2]-I[2]);
			}
			else { 
				//cout << "1.2\n";
				theta_i = acos(ni);
				theta_c = asin(eta_t/eta_i);
				if (theta_c < theta_i && theta_i < PI/2) { //Total Internal Reflection
					//cout << "1.2.1\n";
					T.dx = 0;
					T.dy = 0;
					T.dz = 0;
					not_T = 1;
				}
				else {
					//cout << "1.2.2\n";
					T.dx = -N1[0]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[0]-I[0]);
					T.dy = -N1[1]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[1]-I[1]);
					T.dz = -N1[2]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[2]-I[2]);
				}
			}
			
		}
		// Exit condition
		else if (ni < 0) {
			//cout << "2\n";
			eta_i = mtl[tu].eta;
			eta_i = eta123[eta123.size()-1];
			eta_t = 1;
			eta123.pop_back();
			eta_t = eta123[eta123.size()-1];
			N1[0] = -N[0];
			N1[1] = -N[1];
			N1[2] = -N[2];
			//N[0] = -N[0];
			//N[1] = -N[1];
			//N[2] = -N[2];

			ni = Dot(N1,I);
			Fo = pow((eta_t-eta_i)/(eta_t+eta_i),2);
			Fr = Fo + (1-Fo)*pow((1-ni),5);

			T.x = inter[0] - 10*epsilo*N1[0];
			T.y = inter[1] - 10*epsilo*N1[1];
			T.z = inter[2] - 10*epsilo*N1[2];

			R.x = inter[0] + epsilo*N1[0];
			R.y = inter[1] + epsilo*N1[1];
			R.z = inter[2] + epsilo*N1[2];

			R.dx = 2*ni*N1[0] - I[0];
			R.dy = 2*ni*N1[1] - I[1];
			R.dz = 2*ni*N1[2] - I[2];
			
			if (eta_i <= eta_t) {
				//cout << "2.1\n";
				T.dx = -N1[0]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[0]-I[0]);
				T.dy = -N1[1]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[1]-I[1]);
				T.dz = -N1[2]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[2]-I[2]);		
			}
			else { 
				//cout << "2.2\n";
				theta_i = acos(ni);
				theta_c = asin(eta_t/eta_i);
				if (theta_c < theta_i && theta_i < PI/2) { //Total Internal Reflection
					//cout << "2.2.1\n";
					T.dx = 0;
					T.dy = 0;
					T.dz = 0;
					not_T = 1;
				}
				else {
					//cout << "2.2.2\n";
					T.dx = -N1[0]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[0]-I[0]);
					T.dy = -N1[1]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[1]-I[1]);
					T.dz = -N1[2]*sqrt(1-(pow(eta_i/eta_t,2)*(1-pow(ni,2)))) + (eta_i/eta_t)*(ni*N1[2]-I[2]);
				}
			}
			
		}	
		ColorType R_color, T_color;
		depth++;
		/*
		cout << "Depth: " << depth << endl;
		cout << "Intersection: " << inter[0] << ", " << inter[1] << ", " << inter[2] << endl;
		cout << "Normal: " << N[0] << ", " << N[1] << ", " << N[2] << endl;
		cout << "Incident: " << I[0] << ", " << I[1] << ", " << I[2] << endl;
		cout << "Normal1: " << N1[0] << ", " << N1[1] << ", " << N1[2] << endl;
		cout << "N.I: " << ni << endl;
		cout << "Eta_i: " << eta_i << "		Eta_t: " << eta_t << endl;
		
		cout << "T: " << T.x << " " << T.y << " " << T.z << endl;
		cout << "T_dir: " << T.dx << " " << T.dy << " " << T.dz << endl;
		cout << " \n";
		*/
		if (not_T != 1)
			T_color = Trace_Ray(T, sphere, mtl, light, v_points, face, v_normals, v_tex, tex, eta123, depth);
		else
			T_color = Shade_Ray(token, mtl, tu, N, light, inter, sphere, pa, v_points, face, ij, tex);
		/*
		cout << "Depth: " << depth << endl;
		cout << "R: " << R.x << " " << R.y << " " << R.z << endl;
		cout << "R_dir: " << R.dx << " " << R.dy << " " << R.dz << endl;
		*/
		R_color = Trace_Ray(R, sphere, mtl, light, v_points, face, v_normals, v_tex, tex, eta123, depth);
		color = Shade_Ray(token, mtl, tu, N, light, inter, sphere, pa, v_points, face, ij, tex);
		
		//color.odr += Fr*R_color.odr;
		//color.odg += Fr*R_color.odg;
		//color.odb += Fr*R_color.odb;
		
		color.odr += Fr*R_color.odr + (1-Fr)*(1-mtl[tu].alpha)*T_color.odr;
		color.odg += Fr*R_color.odg + (1-Fr)*(1-mtl[tu].alpha)*T_color.odg;
		color.odb += Fr*R_color.odb + (1-Fr)*(1-mtl[tu].alpha)*T_color.odb;
	}
	else {
		//cout << "Bkg\n";
		//cout << "Depth: " << depth << endl << endl;
		//depth = iterations;
		color = Shade_Ray(token, mtl, tu, N, light, inter, sphere, pa, v_points, face, ij, tex);
	}
	return color;
}

// Main function that accepts scene description file as an input
int main(int argc, char** argv) {
	string filename = argv[1];
	ifstream inFile;

	// Check to see if file exists in the directory
	inFile.open(filename);
	if (!inFile) {
		cerr << "Unable to open file\n";
		exit(1);
	}

	// Check to see if file is empty
	if (inFile.peek() == std::ifstream::traits_type::eof()) {
		cerr << "File is empty\n";
		exit(1);
	}

	// Find the last occurence of '.' in the filename
	std::vector<size_t> vec;
	FindAllOccurences(vec, filename, ".");
	int position = vec[vec.size()-1];

	// New filename will include every occurence of '.' up till the file extension
	string newfile;
	newfile = filename.substr(0, position);
	newfile = newfile + ".ppm";

	// Reading the contents of scene description file
	string str, attrib;
	vector<sphere_t> sphere;
	vector<ColorType> mtl;
	vector<LightType> light;
	vector<size3_t> v_points;
	vector<FaceType> face;
	vector<size3_t> v_normals;
	vector<size2_t> v_tex;
	vector<TexType> tex;
	int m = -1;
	int t = -1;

	cout << "Reading file\n";

	while (inFile) {
		getline(inFile,str);
		if (inFile && !str.empty()) {
			std::vector<size_t> vec;
			FindAllOccurences(vec, str, " ");
			attrib = str.substr(0, vec[0]);

			if (attrib == "sphere") {
				sphere.push_back(AssignSphere(vec, str, m, t));
			}
			else if (attrib == "mtlcolor") {
				m++;
				mtl.push_back(AssignMtlColor(vec, str));
			}
			else if (attrib == "light") {
				light.push_back(AssignLight(vec, str));
			}
			else if (attrib == "v") {
				v_points.push_back(AssignPoints(vec, str));
			}
			else if (attrib == "f") {
				face.push_back(AssignFaces(vec, str, m, t));
			}
			else if (attrib == "vn") {
				v_normals.push_back(AssignNormals(vec, str));
			}
			else if (attrib == "vt") {
				v_tex.push_back(AssignVtex(vec, str));
			}
			else if (attrib == "texture") {
				t++;
				tex.push_back(AssignTextures(vec, str));
			}
			else {
				AssignAttributes(vec, str);
			}
		}
	}

	cout << "Done reading\n";

	/*
	cout << "Total material: " << mtl.size() << endl;
	for (int t = 0; t < mtl.size(); t++) {
		cout << mtl[t].eta << endl;
	}

	*/

	// Extracted parameters from the scene description file
	float view_d[] = {viewdir.x, viewdir.y, viewdir.z};	// Viewing direction
	float up_d[] = {updir.x, updir.y, updir.z};	// Up direction
	float eye_d[] = {eye.x, eye.y, eye.z}; // Eye position
	float width = imsize.w;	// Width of the output image
	float height = imsize.h; // Height of the output image
	float fov = vfov.x;	// Vertical Feild of View
	float aspect_ratio = width/height; // Aspect ratio of the output image

	float cross_d[3]; // Vector that stores any cross product
	float dist;
	float UL[3], UR[3], LL[3], LR[3], xp[3];

	// Calculating 'u' vector
	float u[3];
	Cross(view_d, up_d, cross_d);
	Normalize(cross_d, u);

	if (Dot(view_d,view_d) == 0) {
		cerr << "View direction vector is zero\n";
		exit(1);
	}

	if (Dot(up_d,up_d) == 0) {
		cerr << "Up vector is zero\n";
		exit(1);
	}

	if (fov >= 180) {
		cerr << "Field of View should not be greater than 180 degrees\n";
		exit(1);
	}

	if (isnan(u[0])) {
		cerr << "Up and Viewing direction are parallel. In correct orientation\n";
		exit(1);
	}

	if ((width <= 0) || (height <= 0)) {
		cerr << "Negative value for width/height detected\n";
		exit(1);
	}

	// Calculating 'v' vector
	float v[3];
	Cross(u, view_d, cross_d);
	Normalize(cross_d, v);

	// Calculating viewing distance
	dist = height/(2*tan(fov*PI/(2*180)));

	// Normalize viewing vector
	float view_n[3];
	Normalize(view_d, view_n);

	// Calculating 3D coordinates of the corners of the viewing window
	for (int i = 0; i < 3; i++) {
		UL[i] = eye_d[i] + dist*view_n[i] - (width/2)*u[i] + (height/2)*v[i];
		UR[i] = eye_d[i] + dist*view_n[i] + (width/2)*u[i] + (height/2)*v[i];
		LL[i] = eye_d[i] + dist*view_n[i] - (width/2)*u[i] - (height/2)*v[i];
		LR[i] = eye_d[i] + dist*view_n[i] + (width/2)*u[i] - (height/2)*v[i];
	}

	// Open output file to write data into
	ofstream outFile;
	outFile.open(newfile);
	outFile << "P3" << endl;
	outFile << width << " " << height << endl;
	outFile << 255 << endl;

	float dh, dv, dch, dcv;
	float dir[3], dir_n[3];
	int r,g,b;
	int depth = 0;
	ColorType color;
	vector<float> eta;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			//cout << i << " " << j << endl;
			for (int k = 0; k < 3; k++) {
				dh = (UR[k]-UL[k]) / width;
				dv = (LL[k]-UL[k]) / height;
				dch = dh / 2;
				dcv = dv / 2;
				xp[k] = UL[k] + j*dh + i*dv + dch + dcv; // Point on the viewing window
				dir[k] = xp[k] - eye_d[k];
			}
			Normalize(dir,dir_n); // Normalized direction vector
			ray.x = eye_d[0];
			ray.y = eye_d[1];
			ray.z = eye_d[2];
			ray.dx = dir_n[0];
			ray.dy = dir_n[1];
			ray.dz = dir_n[2];
			depth = 0;
			eta.push_back(1.0);
			// Tracing the vector trajectory to find intersection
			color = Trace_Ray(ray, sphere, mtl, light, v_points, face, v_normals, v_tex, tex, eta, depth);
			eta.clear();
			r = int(255*color.odr);
			g = int(255*color.odg);
			b = int(255*color.odb);
			// Consider overflow condition
			if (r > 255)
				r = 255;
			if (g > 255)
				g = 255;
			if (b > 255)
				b = 255;
			// Writing to output file
			outFile << r << " " << g << " " << b << " ";
		}
		outFile << endl;
	}
	outFile.close();
	inFile.close();
	return 0;
}