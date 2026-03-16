/************************************************************

Name : Praneet Bala
Student ID: 5393094
Date: 09/20/2019
Code: Assignment 1a
Objective: To generate a scene from given scene description file

************************************************************/

// To run: ./trial hello.txt

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#define PI 3.14159265
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
	int m;
};

struct ellipsoid_t {
	float x,y,z;
	float x1,y1,z1;
} ellipse;

struct RayType {
	float x,y,z;
	float dx,dy,dz;
} ray;

struct ColorType {
	float cr,cg,cb;
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

	if (vec.size() == 1) {
		vfov.x = stof(data.substr(vec[0]+1, data.length() - vec[0] - 1));
	}
	else if (vec.size() == 2) {
		imsize.w = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		imsize.h = stof(data.substr(vec[1]+1, data.length() - vec[1] - 1));
	}
	else if (vec.size() == 3) {
		size3_t var1;
		var1.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
		var1.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
		var1.z = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));

		if (attrib == "eye") {
			eye = var1;
		}
		else if (attrib == "viewdir") {
			viewdir = var1;
		}
		else if (attrib == "updir") {
			updir = var1;
		}
		else if (attrib == "bkgcolor") {
			bkgcolor = var1;
		}
	}
}

// Function to extract the sphere information
sphere_t AssignSphere(std::vector<size_t> & vec, std::string data, int m) {
	sphere_t sphere;
	sphere.x = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	sphere.y = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	sphere.z = stof(data.substr(vec[2]+1, vec[3] - vec[2] - 1));
	sphere.r = stof(data.substr(vec[3]+1, data.length() - vec[3] - 1));
	sphere.m = m;
	return sphere;
}

// Function to extract the color information
ColorType AssignMtlColor(std::vector<size_t> & vec, std::string data) {
	ColorType mtlC;
	mtlC.cr = stof(data.substr(vec[0]+1, vec[1] - vec[0] - 1));
	mtlC.cg = stof(data.substr(vec[1]+1, vec[2] - vec[1] - 1));
	mtlC.cb = stof(data.substr(vec[2]+1, data.length() - vec[2] - 1));
	return mtlC;
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
	for (int i = 0; i < 3; i++) {
		N[i] = A[i]/sqrt(Dot(A,A));
	}
}

// Function to assign color to background and foreground pixels
ColorType Shade_Ray(int tok, ColorType mtlcol) {
	ColorType color;
	if (tok == 0) {
		color.cr = bkgcolor.x;
		color.cg = bkgcolor.y;
		color.cb = bkgcolor.z;
	}
	else if (tok == 1) {
		color.cr = mtlcol.cr;
		color.cg = mtlcol.cg;
		color.cb = mtlcol.cb;
	}
	return color;
}

// Function to identify the intersection with sphere
float Trace_Ray(RayType ray, sphere_t sphere) {
	float Bi, Ci, t1, t2, t, discri;
	Bi = 2 * (ray.dx*(ray.x-sphere.x) + ray.dy*(ray.y-sphere.y) + ray.dz*(ray.z-sphere.z));
	Ci = pow((ray.x-sphere.x),2) + pow((ray.y-sphere.y),2) + pow((ray.z-sphere.z),2) - pow(sphere.r,2);
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
	}
	else if (discri == 0) {
		t = -Bi/2;
	}
	else if (discri < 0) {
		t = 100000;
	}
	return t;
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
	string st1 = "sphere";
	string st2 = "mtlcolor";
	int m = -1;
	struct sphere_t sphere[20];
	struct ColorType mtl[20];
	int t_sphere = -1;
	while (inFile) {
		getline(inFile,str);
		if (inFile) {
			std::vector<size_t> vec;
			FindAllOccurences(vec, str, " ");
			attrib = str.substr(0, vec[0]);

			if (st1.compare(attrib) == 0) {
				t_sphere++;
				sphere[t_sphere] = AssignSphere(vec, str, m);
			}
			else if (st2.compare(attrib) == 0) {
				m++;
				mtl[m] = AssignMtlColor(vec, str);
			}
			else {
				AssignAttributes(vec, str);
			}
		}
	}

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
	float min_T, t;
	int token;
	int pa, tu;
	ColorType color;
	
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			min_T = 100000;
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

			// Tracing the vector trajectory to find intersection
			for (int g = 0; g < t_sphere+1; g++) {
				float tt = Trace_Ray(ray, sphere[g]);
				//gt[g] = tt;
				if (tt > 0) {
					if (min(min_T,tt) != min_T) {
						pa = g;
					}
					min_T = min(min_T, tt);
				}
			}
			// Token differentiates between background and foreground pixels
			if (min_T == 100000) {
				token = 0; // Background pixel
				tu = 0;
			}
			else {
				token = 1; // Foreground pixel
				tu = sphere[pa].m;
			}

			color = Shade_Ray(token, mtl[tu]);

			// Varying the background color
			if (token == 0) {
				color.cr = color.cr*(i*2/width);
				color.cg = color.cg*(j/width);
				color.cb = color.cb*((width-i)/width);
			}

			// Writing to output file
			outFile << int(255*color.cr) << " " << int(255*color.cg) << " " << int(255*color.cb) << " ";
		}
		outFile << endl;
	}
	outFile.close();
	inFile.close();
	return 0;
}