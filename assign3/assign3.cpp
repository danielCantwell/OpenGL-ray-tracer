/*
CSCI 420
Assignment 3 Raytracer

Name: <Daniel Cantwell>
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>
#include <math.h>
#include <iostream>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

#define PI 3.14159265

char *filename = 0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode = MODE_DISPLAY;

#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

double dot(double v1[3], double v2[3]);
double area(double v1[2], double v2[2], double v3[2]);

#define cross(result, v1, v2) \
	result[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]); \
	result[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]); \
	result[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);

struct Vertex
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

struct Pixel
{
	double x, y, z;
};

typedef struct _Light
{
	double position[3];
	double color[3];
} Light;

typedef struct _Ray
{
	double origin[3];
	double direction[3];
} Ray;

class Sphere {
public:
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;

	bool rayIntersect(Ray ray, double &minT, bool &inside) {

		double x = ray.origin[0] - position[0];
		double y = ray.origin[1] - position[1];
		double z = ray.origin[2] - position[2];

		// a = 1 = x_direction^2 + y_direction^2 + z_direction^2
		double b = 2 * ((ray.direction[0] * x) + (ray.direction[1] * y) + (ray.direction[2] * z));
		double c = (x * x) + (y * y) + (z * z) - (radius * radius);

		double g = (b * b) - (4 * c);
		if (g < 0) return false;	// no intersection

		// otherwise, find the intersection
		// if g = 0, the ray is tangent to the sphere

		double t0 = (-b - sqrt(g)) / 2;	// point along ray which enters/exits the sphere
		double t1 = (-b + sqrt(g)) / 2;	// point along ray which exits/enters the sphere

		// determine which point is closer
		if (t0 < 0) {
			minT = t1;
			inside = true;
		}
		else if (t1 < 0) {
			minT = t0;
			inside = true;
		}
		else {
			minT = min(t0, t1);
		}

		return true;
	}
};

class Triangle {
public:
	struct Vertex vertex[3];

	bool rayIntersect(Ray ray, double &t, double &alpha, double &beta, double &gamma) {

		double edge1[3], edge2[3];
		/* calculate edges */
		for (int i = 0; i < 3; i++) {
			edge1[i] = vertex[1].position[i] - vertex[0].position[i];
			edge2[i] = vertex[2].position[i] - vertex[0].position[i];
		}

		/* calculate normal of plane containing triangle */
		double planeNormal[3];
		/* planeNormal = edge1 x edge2 */
		cross(planeNormal, edge1, edge2);
		/* normalize normal */
		double mag = sqrt((planeNormal[0] * planeNormal[0])
			+ (planeNormal[1] * planeNormal[1])
			+ (planeNormal[2] * planeNormal[2]));
		planeNormal[0] /= mag;
		planeNormal[1] /= mag;
		planeNormal[2] /= mag;

		/* calculate intersection of ray with plane */
		/* t = - (origin - vertex) .dot. normal / direction .dot. normal */
		double ov[3];	// origin - vertex
		ov[0] = ray.origin[0] - vertex[0].position[0];
		ov[1] = ray.origin[1] - vertex[0].position[1];
		ov[2] = ray.origin[2] - vertex[0].position[2];

		// ray = O + tD
		double DdotN = dot(ray.direction, planeNormal);
		if (DdotN == 0) return false;	// ray parallel to plane

		t = -dot(ov, planeNormal) / DdotN;
		// QUESTIONABLE >>>
		if (t < 0.00001) return false;	// intersection is behind origin

		/* calculate intersection point */
		double point[3];
		point[0] = ray.origin[0] + (t * ray.direction[0]);
		point[1] = ray.origin[1] + (t * ray.direction[1]);
		point[2] = ray.origin[2] + (t * ray.direction[2]);

		/* test is point is inside triangle using barycentric coordinates */
		double pointProjected[2], v0Projected[2], v1Projected[2], v2Projected[2];

		for (int i = 0; i < 2; i++) {
			pointProjected[i] = point[i];
			v0Projected[i] = vertex[0].position[i];
			v1Projected[i] = vertex[1].position[i];
			v2Projected[i] = vertex[2].position[i];
		}

		double triArea = area(v0Projected, v1Projected, v2Projected);

		alpha = area(pointProjected, v1Projected, v2Projected) / triArea;
		beta = area(v0Projected, pointProjected, v2Projected) / triArea;
		gamma = area(v0Projected, v1Projected, pointProjected) / triArea;

		if (alpha * beta >= 0 && beta * gamma >= 0) {
			/* point is in triangle */
			return true;
		}
		else {
			return false;
		}
	}
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void save_jpg();

double aspectRatio;
// image corner coordinates
double xLeft;	// left
double xRight;	// right
double yBottom;	// bottom
double yTop;	// top
double z;

Pixel pixels[WIDTH][HEIGHT];

void calculateCorners() {
	aspectRatio = (double)WIDTH / (double)HEIGHT;

	xLeft = -aspectRatio * tan(fov * PI / 360); //	fov/2 * PI/180
	xRight = -xLeft;

	yBottom = -tan(fov * PI / 360);
	yTop = -yBottom;

	z = -1;
}

/* Calculate the area of a triangle in 2D */
double area(double v1[2], double v2[2], double v3[2]) {

	double abx = v2[0] - v1[0];
	double aby = v2[1] - v1[1];

	double acx = v3[0] - v1[0];
	double acy = v3[1] - v1[1];

	return ((abx * acy) - (acx * aby)) / 2;
}

/* Calculate the Dot Product of two vertices */
double dot(double v1[3], double v2[3]) {
	double result = 0;
	for (int i = 0; i < 3; i++) {
		result += v1[i] * v2[i];
	}
	return result;
}

Ray createShadowRay(double p[3], Light l) {
	Ray shadowRay;

	// origin
	shadowRay.origin[0] = p[0];
	shadowRay.origin[1] = p[1];
	shadowRay.origin[2] = p[2];
	// direction
	shadowRay.direction[0] = l.position[0] - p[0];
	shadowRay.direction[1] = l.position[1] - p[1];
	shadowRay.direction[2] = l.position[2] - p[2];
	// magnitude
	double xShadowMag = shadowRay.direction[0];
	double yShadowMag = shadowRay.direction[1];
	double zShadowMag = shadowRay.direction[2];
	double mag = sqrt((xShadowMag * xShadowMag) + (yShadowMag * yShadowMag) + (zShadowMag * zShadowMag));
	// normalize
	shadowRay.direction[0] /= mag;
	shadowRay.direction[1] /= mag;
	shadowRay.direction[2] /= mag;

	return shadowRay;
}

Pixel calculateSphereColor(Sphere s, Ray shadowRay, Light l, double pIntersect[3], bool insideSphere) {

	double nIntersect[3];
	/* calculate normal of intersection point on sphere */
	nIntersect[0] = (pIntersect[0] - s.position[0]) / s.radius;
	nIntersect[1] = (pIntersect[1] - s.position[1]) / s.radius;
	nIntersect[2] = (pIntersect[2] - s.position[2]) / s.radius;

	/* negate normal if ray starts inside of sphere */
	if (insideSphere) {
		nIntersect[0] *= -1;
		nIntersect[1] *= -1;
		nIntersect[2] *= -1;
	}

	Ray L = shadowRay;
	double LdotN = dot(L.direction, nIntersect);
	LdotN = (LdotN < 0) ? 0 : ((LdotN > 1) ? 1 : LdotN);

	Ray R = L;
	R.direction[0] = 2 * LdotN * nIntersect[0] - L.direction[0];
	R.direction[1] = 2 * LdotN * nIntersect[1] - L.direction[1];
	R.direction[2] = 2 * LdotN * nIntersect[2] - L.direction[2];
	// normalize
	double xRMag = R.direction[0];
	double yRMag = R.direction[1];
	double zRMag = R.direction[2];
	double Rmag = sqrt((xRMag * xRMag) + (yRMag * yRMag) + (zRMag * zRMag));
	R.direction[0] /= Rmag;
	R.direction[1] /= Rmag;
	R.direction[2] /= Rmag;

	// pIntersect normalize
	double pMag = sqrt((pIntersect[0] * pIntersect[0]) + (pIntersect[1] * pIntersect[1]) + (pIntersect[2] * pIntersect[2]));
	pIntersect[0] /= -pMag;
	pIntersect[1] /= -pMag;
	pIntersect[2] /= -pMag;

	double RdotV = dot(R.direction, pIntersect);
	RdotV = (RdotV < 0) ? 0 : ((RdotV > 1) ? 1 : RdotV);

	/* calculate light intensity */
	double colorX = l.color[0] * ((s.color_diffuse[0] * LdotN) + s.color_specular[0] * pow(RdotV, s.shininess));
	double colorY = l.color[1] * ((s.color_diffuse[1] * LdotN) + s.color_specular[1] * pow(RdotV, s.shininess));
	double colorZ = l.color[2] * ((s.color_diffuse[2] * LdotN) + s.color_specular[2] * pow(RdotV, s.shininess));

	/* clamp to 0,1 */
	/* if color < 0, color = 0, else if color > 1, color = 1, else no change*/
	colorX = (colorX < 0) ? 0 : ((colorX > 1) ? 1 : colorX);
	colorY = (colorY < 0) ? 0 : ((colorY > 1) ? 1 : colorY);
	colorZ = (colorZ < 0) ? 0 : ((colorZ > 1) ? 1 : colorZ);

	Pixel p;
	p.x = colorX;
	p.y = colorY;
	p.z = colorZ;

	return p;
}

Pixel calculateTriangleColor(Triangle t, Ray shadowRay, Light l, double diffuse_color[3], double pIntersect[3], 
	double normal[3], double alpha, double beta, double gamma) {

	Ray L = shadowRay;
	double LdotN = dot(L.direction, normal);
	LdotN = (LdotN < 0) ? 0 : ((LdotN > 1) ? 1 : LdotN);

	Ray R = L;
	R.direction[0] = 2 * LdotN * normal[0] - L.direction[0];
	R.direction[1] = 2 * LdotN * normal[1] - L.direction[1];
	R.direction[2] = 2 * LdotN * normal[2] - L.direction[2];
	// normalize
	double xRMag = R.direction[0];
	double yRMag = R.direction[1];
	double zRMag = R.direction[2];
	double Rmag = sqrt((xRMag * xRMag) + (yRMag * yRMag) + (zRMag * zRMag));
	R.direction[0] /= Rmag;
	R.direction[1] /= Rmag;
	R.direction[2] /= Rmag;

	// pIntersect normalize
	double pMag = sqrt((pIntersect[0] * pIntersect[0]) + (pIntersect[1] * pIntersect[1]) + (pIntersect[2] * pIntersect[2]));
	pIntersect[0] /= -pMag;
	pIntersect[1] /= -pMag;
	pIntersect[2] /= -pMag;

	double RdotV = dot(R.direction, pIntersect);
	RdotV = (RdotV < 0) ? 0 : ((RdotV > 1) ? 1 : RdotV);

	/* calculate light intensity */
	double colorX = l.color[0] * ((diffuse_color[0] * LdotN));// +specular_color[0] * pow(RdotV, shininess));
	double colorY = l.color[1] * ((diffuse_color[1] * LdotN));// +specular_color[1] * pow(RdotV, shininess));
	double colorZ = l.color[2] * ((diffuse_color[2] * LdotN));// +specular_color[2] * pow(RdotV, shininess));

	/* clamp to 0,1 */
	/* if color < 0, color = 0, else if color > 1, color = 1, else no change*/
	colorX = (colorX < 0) ? 0 : ((colorX > 1) ? 1 : colorX);
	colorY = (colorY < 0) ? 0 : ((colorY > 1) ? 1 : colorY);
	colorZ = (colorZ < 0) ? 0 : ((colorZ > 1) ? 1 : colorZ);

	Pixel p;
	p.x = colorX;
	p.y = colorY;
	p.z = colorZ;

	return p;
}

Pixel rayTrace(Ray ray) {

	Pixel pixel;
	pixel.x = 1;
	pixel.y = 1;
	pixel.z = 1;

	Sphere s;
	double minSphereDistance = -1;
	bool insideSphere = false;

	/* test if the ray intersects a sphere */
	for (int i = 0; i < num_spheres; i++) {
		Sphere sphere = spheres[i];
		double minT;
		if (sphere.rayIntersect(ray, minT, insideSphere)) {
			if (minSphereDistance == -1 || minT < minSphereDistance) {
				minSphereDistance = minT;
				s = sphere;
			}
		}
	}

	Triangle t;
	double alpha, beta, gamma;
	double minTriangleDistance = -1;

	/* test if the ray intersects a triangle */
	for (int i = 0; i < num_triangles; i++) {
		Triangle triangle = triangles[i];
		double minT;
		if (triangle.rayIntersect(ray, minT, alpha, beta, gamma)) {
			if (minTriangleDistance == -1 || minT < minTriangleDistance) {
				minTriangleDistance = minT;
				t = triangle;
			}
		}
	}

	/* if a ray intersects both a sphere and a triangle, determine which is closer */
	if (minSphereDistance != -1 && minTriangleDistance != -1) {
		if (minSphereDistance <= minTriangleDistance) {
			minTriangleDistance = -1;
		}
		else {
			minSphereDistance = -1;
		}
	}

	/********************************************************************************/
	/*********************    if the ray intersects a sphere    *********************/
	/********************************************************************************/
	if (minSphereDistance != -1) {

		pixel.x = ambient_light[0];
		pixel.y = ambient_light[1];
		pixel.z = ambient_light[2];

		double pIntersect[3];

		/* calculate intersection point on sphere */
		pIntersect[0] = ray.origin[0] + ray.direction[0] * minSphereDistance;
		pIntersect[1] = ray.origin[1] + ray.direction[1] * minSphereDistance;
		pIntersect[2] = ray.origin[2] + ray.direction[2] * minSphereDistance;

		/* for each light source, fire shadow ray */
		for (int i = 0; i < num_lights; i++) {
			Light l = lights[i];

			/* unit vector to the light */
			Ray shadowRay = createShadowRay(pIntersect, l);

			///////////// this should not be here ////////////////////////////////////////////////////////////
			Pixel pixelColor = calculateSphereColor(s, shadowRay, l, pIntersect, insideSphere);
			pixel.x += pixelColor.x;
			pixel.y += pixelColor.y;
			pixel.z += pixelColor.z;

			if (true) {
				/* for each unblocked shadow ray, evaluate local phong model for */
				/* that light, and add result to pixel color */
				for (int j = 0; j < num_spheres; j++) {
					Sphere sphere = spheres[j];
					/* determine if shadow ray is blocked by any objects */
					double dummy1;
					bool dummy2;
					if (!(sphere.rayIntersect(shadowRay, dummy1, dummy2))) {

						for (int k = 0; k < num_triangles; k++) {
							Triangle triangle = triangles[k];
							double d1, d2, d3, d4;
							if (!(triangle.rayIntersect(shadowRay, d1, d2, d3, d4))) {

								/* if no objects are blocking the shadow ray, calculate color */

								Pixel pixelColor = calculateSphereColor(s, shadowRay, l, pIntersect, insideSphere);

								pixel.x += pixelColor.x;
								pixel.y += pixelColor.y;
								pixel.z += pixelColor.z;
							}
						}

					}
				}
			}
		}

	}

	/********************************************************************************/
	/********************    if the ray intersects a triangle    ********************/
	/********************************************************************************/
	else if (minTriangleDistance != -1) {

		pixel.x = ambient_light[0];
		pixel.y = ambient_light[1];
		pixel.z = ambient_light[2];

		double point[3];
		point[0] = ray.origin[0] + (minTriangleDistance * ray.direction[0]);
		point[1] = ray.origin[1] + (minTriangleDistance * ray.direction[1]);
		point[2] = ray.origin[2] + (minTriangleDistance * ray.direction[2]);

		//double pIntersect[3];
		///* calculate intersection point on triangle */
		//for (int i = 0; i < 3; i++) {
		//	pIntersect[i] = (alpha * t.vertex[0].position[i])
		//		+ (beta * t.vertex[1].position[i])
		//		+ (gamma * t.vertex[2].position[i]);
		//}

		/* calculate interpolated diffuse color */
		double diffuse_color[3];
		for (int i = 0; i < 3; i++) {
			diffuse_color[i] = ((alpha * t.vertex[0].color_diffuse[i])
				+ (beta * t.vertex[1].color_diffuse[i])
				+ (gamma * t.vertex[2].color_diffuse[i]));
		}

		/*std::cout << std::endl;
		std::cout << "point[0] : " << point[0] << std::endl;
		std::cout << "point[1] : " << point[1] << std::endl;
		std::cout << "point[2] : " << point[2] << std::endl;
		std::cout << "intersect[0] : " << pIntersect[0] << std::endl;
		std::cout << "intersect[1] : " << pIntersect[1] << std::endl;
		std::cout << "intersect[2] : " << pIntersect[2] << std::endl;
		std::cout << std::endl;*/

		/* calculate interpolated normal of intersection point on triangle */
		double normal[3];
		for (int i = 0; i < 3; i++) {
			normal[i] = ((alpha * t.vertex[0].normal[i])
				+ (beta * t.vertex[1].normal[i])
				+ (gamma * t.vertex[2].normal[i]));
		}

		/* calculate interpolated specular color */
		double specular_color[3];
		for (int i = 0; i < 3; i++) {
			specular_color[i] = ((alpha * t.vertex[0].color_specular[i])
				+ (beta	* t.vertex[1].color_specular[i])
				+ (gamma * t.vertex[2].color_specular[i]));
		}

		/* calculate interpolated shininess */
		double shininess = (alpha * t.vertex[0].shininess) + (beta	* t.vertex[1].shininess) + (gamma * t.vertex[2].shininess);

		/* for each light source, fire shadow ray */
		for (int i = 0; i < num_lights; i++) {
			Light l = lights[i];

			/* unit vector to the light */
			Ray shadowRay = createShadowRay(point, l);

			//Pixel pixelColor = calculateTriangleColor(t, shadowRay, l, diffuse_color, point, normal, alpha, beta, gamma);

			//pixel.x += pixelColor.x;
			//pixel.y += pixelColor.y;
			//pixel.z += pixelColor.z;

			if (true) {
				/* for each unblocked shadow ray, evaluate local phong model for */
				/* that light, and add result to pixel color */
				for (int j = 0; j < num_spheres; j++) {
					Sphere sphere = spheres[j];
					/* determine if shadow ray is blocked by any objects */
					double dummy1;
					bool dummy2;
					if (!(sphere.rayIntersect(shadowRay, dummy1, dummy2))) {
						for (int j = 0; j < num_triangles; j++) {
							Triangle triangle = triangles[j];
							/* determine if shadow ray is blocked by any objects */
							double dummy3, dummy4, dummy5, dummy6;
							if (!(triangle.rayIntersect(shadowRay, dummy3, dummy4, dummy5, dummy6))) {
								/* if no objects are blocking the shadow ray, calculate color */

								Pixel pixelColor = calculateTriangleColor(t, shadowRay, l, diffuse_color, point, normal, alpha, beta, gamma);

								pixel.x += pixelColor.x;
								pixel.y += pixelColor.y;
								pixel.z += pixelColor.z;
							}
						}
					}
				}
			}
		}
	}

	return pixel;
}

void castRays() {

	std::cout << "\ncasting rays\n";

	// ray origin = camera origin = 0,0,0
	Ray ray;
	ray.origin[0] = 0;
	ray.origin[1] = 0;
	ray.origin[2] = 0;

	// calculate image width, height and pixel separation
	double imageWidth = xRight - xLeft;	// right - left
	double imageHeight = yTop - yBottom;	// top - bottom
	double stepWidth = imageWidth / WIDTH;
	double stepHeight = imageHeight / HEIGHT;

	std::cout << "image width : " << stepWidth << std::endl;
	std::cout << "image height : " << imageHeight << std::endl;

	// j = top; while j > bottom; decrement j
	for (int n = 0; n < HEIGHT; n++){
		for (int m = 0; m < WIDTH; m++)
		{
			double x = xLeft + (1.0 * m + 0.5) * stepWidth;
			double y = yBottom + (1.0 * n + 0.5) * stepHeight;
		
			double magnitude = sqrt(x*x + y*y + 1); // z = -1, so z*z = 1

			ray.direction[0] = x / magnitude;
			ray.direction[1] = y / magnitude;
			ray.direction[2] = -1 / magnitude;

			pixels[m][n] = rayTrace(ray);
		}
	}

}

void draw_scene()
{
	calculateCorners();
	castRays();

	unsigned int x, y;
	//simple output
	for (x = 0; x < WIDTH; x++)
	{
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (y = 0; y < HEIGHT; y++)
		{
			//plot_pixel(x, y, x % 256, y % 256, (x + y) % 256);
			glColor3f(pixels[x][y].x, pixels[x][y].y, pixels[x][y].z);
			glVertex2f(x, y);
		}
		glEnd();
		glFlush();
	}

	printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	glColor3f(((double)r) / 256.f, ((double)g) / 256.f, ((double)b) / 256.f);
	glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	buffer[HEIGHT - y - 1][x][0] = r;
	buffer[HEIGHT - y - 1][x][1] = g;
	buffer[HEIGHT - y - 1][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	plot_pixel_display(x, y, r, g, b);
	if (mode == MODE_JPEG)
		plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
	Pic *in = NULL;

	in = pic_alloc(640, 480, 3, NULL);
	printf("Saving JPEG file: %s\n", filename);

	memcpy(in->pix, buffer, 3 * WIDTH*HEIGHT);
	if (jpeg_write(filename, in))
		printf("File saved Successfully\n");
	else
		printf("Error in Saving\n");

	pic_free(in);

}

void parse_check(char *expected, char *found)
{
	if (stricmp(expected, found))
	{
		char error[100];
		printf("Expected '%s ' found '%s '\n", expected, found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}

}

void parse_doubles(FILE*file, char *check, double p[3])
{
	char str[100];
	fscanf(file, "%s", str);
	parse_check(check, str);
	fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
	printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE*file, double *r)
{
	char str[100];
	fscanf(file, "%s", str);
	parse_check("rad:", str);
	fscanf(file, "%lf", r);
	printf("rad: %f\n", *r);
}

void parse_shi(FILE*file, double *shi)
{
	char s[100];
	fscanf(file, "%s", s);
	parse_check("shi:", s);
	fscanf(file, "%lf", shi);
	printf("shi: %f\n", *shi);
}

int loadScene(char *argv)
{
	FILE *file = fopen(argv, "r");
	int number_of_objects;
	char type[50];
	int i;
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file, "%i", &number_of_objects);

	printf("number of objects: %i\n", number_of_objects);
	char str[200];

	parse_doubles(file, "amb:", ambient_light);

	for (i = 0; i < number_of_objects; i++)
	{
		fscanf(file, "%s\n", type);
		printf("%s\n", type);
		if (stricmp(type, "triangle") == 0)
		{

			printf("found triangle\n");
			int j;

			for (j = 0; j < 3; j++)
			{
				parse_doubles(file, "pos:", t.vertex[j].position);
				parse_doubles(file, "nor:", t.vertex[j].normal);
				parse_doubles(file, "dif:", t.vertex[j].color_diffuse);
				parse_doubles(file, "spe:", t.vertex[j].color_specular);
				parse_shi(file, &t.vertex[j].shininess);
			}

			if (num_triangles == MAX_TRIANGLES)
			{
				printf("too many triangles, you should increase MAX_TRIANGLES!\n");
				exit(0);
			}
			triangles[num_triangles++] = t;
		}
		else if (stricmp(type, "sphere") == 0)
		{
			printf("found sphere\n");

			parse_doubles(file, "pos:", s.position);
			parse_rad(file, &s.radius);
			parse_doubles(file, "dif:", s.color_diffuse);
			parse_doubles(file, "spe:", s.color_specular);
			parse_shi(file, &s.shininess);

			if (num_spheres == MAX_SPHERES)
			{
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
		else if (stricmp(type, "light") == 0)
		{
			printf("found light\n");
			parse_doubles(file, "pos:", l.position);
			parse_doubles(file, "col:", l.color);

			if (num_lights == MAX_LIGHTS)
			{
				printf("too many lights, you should increase MAX_LIGHTS!\n");
				exit(0);
			}
			lights[num_lights++] = l;
		}
		else
		{
			printf("unknown type in scene description:\n%s\n", type);
			exit(0);
		}
	}
	return 0;
}

void display()
{
}

void init()
{
	glMatrixMode(GL_PROJECTION);
	glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
	//hack to make it only draw once
	static int once = 0;
	if (!once)
	{
		draw_scene();
		if (mode == MODE_JPEG)
			save_jpg();
	}
	once = 1;
}

int main(int argc, char ** argv)
{
	if (argc<2 || argc > 3)
	{
		printf("usage: %s <scenefile> [jpegname]\n", argv[0]);
		exit(0);
	}
	if (argc == 3)
	{
		mode = MODE_JPEG;
		filename = argv[2];
	}
	else if (argc == 2)
		mode = MODE_DISPLAY;

	glutInit(&argc, argv);
	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIDTH, HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	init();
	glutMainLoop();
}
