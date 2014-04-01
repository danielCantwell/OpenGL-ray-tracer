/*
CSCI 480
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

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
//#define WIDTH 320
//#define HEIGHT 240

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

struct point
{
	double x, y, z;
};

struct normal
{
	double x, y, z;
};

struct Pixel
{
	double x, y, z;
};

typedef struct _Triangle
{
	struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;
} Sphere;

typedef struct _Light
{
	double position[3];
	double color[3];
} Light;

typedef struct _Ray
{
	point origin;
	point direction;
} Ray;

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

double aspectRatio;
// image corner coordinates
// x
double x00;	// bottom left
double x10;	// bottom right
double x01;	// top left
double x11;	// top right
// y
double y00;	// bottom left
double y10;	// bottom right
double y01;	// top left
double y11;	// top right
// z
double z00;	// bottom left
double z10;	// bottom right
double z01;	// top left
double z11;	// top right

Pixel pixels[WIDTH][HEIGHT];

void calculateCorners() {
	aspectRatio = (double)WIDTH / (double)HEIGHT;

	x00 = -aspectRatio * tan(fov * PI / 360); //	fov/2 * PI/180
	x10 = -x00;
	x01 = x00;
	x11 = x10;

	y00 = -tan(fov * PI / 360);
	y10 = y00;
	y01 = -y00;
	y11 = y01;

	z00 = -1;
	z10 = -1;
	z01 = -1;
	z11 = -1;

	std::cout << std::endl << std::endl;
	std::cout << "x00 : " << x00 << std::endl;
	std::cout << "x10 : " << x10 << std::endl;
	std::cout << "x01 : " << x01 << std::endl;
	std::cout << "x11 : " << x11 << std::endl;
	std::cout << "y00 : " << y00 << std::endl;
	std::cout << "y10 : " << y10 << std::endl;
	std::cout << "y01 : " << y01 << std::endl;
	std::cout << "y11 : " << y11 << std::endl;
	std::cout << "z00 : " << z00 << std::endl;
	std::cout << "z10 : " << z10 << std::endl;
	std::cout << "z01 : " << z01 << std::endl;
	std::cout << "z11 : " << z11 << std::endl;
}

bool findIntersectionWithSpheres(Ray ray, point &pIntersect, normal &nIntersect, Sphere &sphere, int &sIndex) {

	bool intersects = false;
	double minDistance = 1000000000;	// arbitrary large number

	// determine closest sphere intersection
	for (int i = 0; i < num_spheres; i++) {
		Sphere s = spheres[i];

		double x = ray.origin.x - s.position[0];
		double y = ray.origin.y - s.position[1];
		double z = ray.origin.z - s.position[2];

		// a = 1 = x_direction^2 + y_direction^2 + z_direction^2
		double b = 2 * ((ray.direction.x * x) + (ray.direction.y * y) + (ray.direction.z * z));
		double c = (x * x) + (y * y) + (z * z) - (s.radius * s.radius);

		double g = (b * b) - (4 * c);
		if (g < 0) continue;	// no intersection

		// otherwise, find the intersection
		// if g = 0, the ray is tangent to the sphere

		double t0 = (-b - sqrt(g)) / 2;	// point along ray which enters the sphere
		double t1 = (-b + sqrt(g)) / 2;	// point along ray which exits the sphere

		// determine which point is closer
		double minT = min(t0, t1);

		if (minT < minDistance) {
			sphere = s;
			sIndex = i;
			minDistance = minT;
			intersects = true;
		}

	}

	if (intersects) {
		pIntersect.x = ray.origin.x + ray.direction.x*minDistance;
		pIntersect.y = ray.origin.y + ray.direction.y*minDistance;
		pIntersect.z = ray.origin.z + ray.direction.z*minDistance;

		nIntersect.x = (pIntersect.x - sphere.position[0]) / sphere.radius;
		nIntersect.y = (pIntersect.y - sphere.position[1]) / sphere.radius;
		nIntersect.z = (pIntersect.z - sphere.position[2]) / sphere.radius;
	}

	if (minDistance < 0) {
		nIntersect.x *= -1;
		nIntersect.y *= -1;
		nIntersect.z *= -1;
	}

	return intersects;
}

bool intersectsObject(Ray ray, int fromIndex) {

	// determine sphere intersections
	for (int i = 0; i < num_spheres; i++) {

		if (i == fromIndex) continue;

		Sphere s = spheres[i];

		double x = ray.origin.x - s.position[0];
		double y = ray.origin.y - s.position[1];
		double z = ray.origin.z - s.position[2];

		// a = 1 = x_direction^2 + y_direction^2 + z_direction^2
		double b = 2 * ((ray.direction.x * x) + (ray.direction.y * y) + (ray.direction.z * z));
		double c = (x * x) + (y * y) + (z * z) - (s.radius * s.radius);

		double g = (b * b) - (4 * c);
		if (g < 0) continue;	// no intersection

		// otherwise, there is an intersection
		return true;
	}

	return false;
}

Pixel rayTrace(Ray ray) {

	Pixel pixel;
	pixel.x = 0;
	pixel.y = 0;
	pixel.z = 0;

	point pIntersect;
	normal nIntersect;
	Sphere s;
	int sIndex;
	
	// if the ray intersects a sphere
	if (findIntersectionWithSpheres(ray, pIntersect, nIntersect, s, sIndex)) {
		pixel.x = s.color_diffuse[0];
		pixel.y = s.color_diffuse[1];
		pixel.z = s.color_diffuse[2];
		// for each light source, fire shadow ray
		for (int i = 0; i < num_lights; i++) {
			Light l = lights[i];
			Ray shadowRay;
			shadowRay.origin.x = pIntersect.x;
			shadowRay.origin.y = pIntersect.y;
			shadowRay.origin.z = pIntersect.z;

			shadowRay.direction.x = l.position[0] - pIntersect.x;
			shadowRay.direction.y = l.position[1] - pIntersect.y;
			shadowRay.direction.z = l.position[2] - pIntersect.z;

			double mag = sqrt((shadowRay.direction.x * shadowRay.direction.x) +
				(shadowRay.direction.y * shadowRay.direction.y) +
				(shadowRay.direction.z * shadowRay.direction.z));

			// normalize
			shadowRay.direction.x /= mag;
			shadowRay.direction.y /= mag;
			shadowRay.direction.z /= mag;

			// for each unblocked shadow ray, evaluate local phong model for
			// that light, and add result to pixel color
			if (!intersectsObject(shadowRay, sIndex)) {

				Ray L = shadowRay;

				point LdotN;
				LdotN.x = L.direction.x * nIntersect.x;
				LdotN.y = L.direction.y * nIntersect.y;
				LdotN.z = L.direction.z * nIntersect.z;

				/*
				Ray R = L;
				R.direction.x -= pIntersect.x;
				R.direction.y -= pIntersect.y;
				R.direction.z -= pIntersect.z;
				*/				
				Ray R = L;
				R.direction.x = 2 * LdotN.x * nIntersect.x - L.direction.x;
				R.direction.y = 2 * LdotN.y * nIntersect.y - L.direction.y;
				R.direction.z = 2 * LdotN.z * nIntersect.z - L.direction.z;
				double Rmag = sqrt((R.direction.x * R.direction.x)
					+ (R.direction.y * R.direction.y)
					+ (R.direction.z * R.direction.z));
				R.direction.x /= Rmag;
				R.direction.y /= Rmag;
				R.direction.z /= Rmag;

				point RdotV;
				RdotV.x = R.direction.x * -pIntersect.x;
				RdotV.y = R.direction.y * -pIntersect.y;
				RdotV.z = R.direction.z * -pIntersect.z;
				//std::cout << "  x == " << shadowRay.direction.x + shadowRay.direction.y + shadowRay.direction.z;

				double colorX = l.color[0] * ((s.color_diffuse[0] * LdotN.x) + (s.color_specular[0] * RdotV.x));
				double colorY = l.color[1] * ((s.color_diffuse[1] * LdotN.y) + (s.color_specular[1] * RdotV.y));
				double colorZ = l.color[2] * ((s.color_diffuse[2] * LdotN.z) + (s.color_specular[2] * RdotV.z));

				pixel.x += colorX;
				pixel.y += colorY;
				pixel.z += colorZ;
			}
			else {
				pixel.x = s.color_diffuse[0];
				pixel.y = s.color_diffuse[1];
				pixel.z = s.color_diffuse[2];
			}
		}
	}
	// if the ray does not intersect a sphere
	else {
		pixel.x = 0.4;
		pixel.y = 0.4;
		pixel.z = 0.4;
	}

	
	return pixel;
}

void castRays() {

	std::cout << "\ncasting rays\n";

	// ray origin = camera origin = 0,0,0
	Ray ray;
	ray.origin.x = 0;
	ray.origin.y = 0;
	ray.origin.z = 0;

	
	// calculate image width, height and pixel separation
	double imageWidth = x10 - x00;	// right - left
	double imageHeight = y01 - y00;	// top - bottom
	double stepWidth = imageWidth / WIDTH;
	double stepHeight = imageHeight / HEIGHT;

	std::cout << "image width : " << stepWidth << std::endl;
	std::cout << "image height : " << imageHeight << std::endl;

	// j = top, while j > bottom, decrement j
	int m = 0, n = 0;
	for (double j = y01 - stepHeight / 2; j > y00 + stepHeight; j -= stepHeight, n++) {
		// i = left, while i < right, increment i
		m = 0;
		for (double i = x00 + stepWidth / 2; i < x10 - stepWidth; i += stepWidth, m++) {
			double magnitude = sqrt(i*i + j*j + 1); // z = -1, so z*z = 1
			ray.direction.x = i / magnitude;
			ray.direction.y = j / magnitude;
			ray.direction.z = -1 / magnitude;

			// trace the ray and assign color
			pixels[m][n] = rayTrace(ray);
		}
	}
	
}

//MODIFY THIS FUNCTION
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
				parse_doubles(file, "pos:", t.v[j].position);
				parse_doubles(file, "nor:", t.v[j].normal);
				parse_doubles(file, "dif:", t.v[j].color_diffuse);
				parse_doubles(file, "spe:", t.v[j].color_specular);
				parse_shi(file, &t.v[j].shininess);
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
