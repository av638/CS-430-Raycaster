// CS 430 Project 2 Raycaster
// In this project we are tasked with not only being able to take in a JSON file and parse it into
// a structure that will accommodate  all valid objects, but we must also take those objects and visualize them
// From here we must use the ppm converter that we originally created in the first project and using the same
// idea we will have to print out the scene/objects to a ppm file.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "JSONparser.c"
#include "ppmWrite.c"

Object objects[128];

static inline double sqr(double v)
{
  return v*v;
}


static inline void normalize(double* v)
{
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

// Takes in Ro - ray origin, Rd - ray direction, p - position of plane, and n - normal of the plane
// If there is no intersection with the plane then we will return -1 otherwise return a t
double planeIntersection(double* p, double* n, double* Rd, double* Ro){

	double top, bottom, t;

	top = (n[0] * Ro[0] - n[0] * p[0] + n[1] * Ro[1] - n[1] * p[1] + n[2] * Ro[2] - n[2] * p[2]);
	bottom = (- 1 * ((n[0] * Rd[0] + n[1] * Rd[1] + n[2] * Rd[2])));
	t = top / bottom;

	if(t > 0) return t;

	return (-1);
}


// Takes in Ro - ray origin, Rd - ray direction, p - position or center sphere, and r - radius
// If at any point we find that there is no intersection with the sphere then we
// must make sure to return -1 that way we know whether or not to color in the pixel or not
double sphereIntersection(double* p, double r,
                           double* Rd, double* Ro)
{
    double a, b, c, t, det;

    a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
    b =  2 * (Rd[0] * (Ro[0] - p[0]) + Rd[1] * (Ro[1] - p[1]) + Rd[2] * (Ro[2] - p[2]));
    c = sqr(Ro[0]) - (2 * Ro[0] * p[0] + sqr(p[0]) + sqr(Ro[1])) - (2 * Ro[1] * p[1] + sqr(p[1]) + sqr(Ro[2])) - (2 * Ro[2] * p[2] + sqr(p[2])) - sqr(r);
    det = b*b-4*c;

    t = det;
    //printf("%d", b);
    //exit(0);
    //det = sqrt(det);
    //t = (-b - det) / ( 2 * a);

    if(t < 0) return -1;

    return t;
}
// Raycaster will take in the array of objects from the JSON file, the width and height desired
// as well as the actually number of objects that were read in in order to loop through the array
// of objects
// Will return the buffer holding the image that we want to print out to the ppm file
Pixmap* rayCaster(Object objects[], int width, int height, int numObjects)
{
	double cx, cy, h, w, pixelHeight, pixelWidth;
	int i, t, x, y;

	PixelColor pixel;
    Pixmap *buffer = (Pixmap *)malloc(sizeof(Pixmap));

	cx = 0;
	cy = 0;

	buffer->width = width;
	buffer->height = height;
	buffer->color = 255;
	// Allocated memory size for image
	buffer->image = malloc(sizeof(Pixmap) * buffer->width * buffer->height);

	for(i = 0; i < numObjects; i++)
    {
		if(strcmp((objects[i].type), "camera") == 0)
        {
			h = objects[i].properties.camera.height;
			w = objects[i].properties.camera.width;
		}
	}

	pixelHeight = h / height;
    pixelWidth = w / width;

    // Start going row by column and checking whether each pixel intersects with an object or not
	for (y = 0; y < height; y++)
    {
		for (x = 0; x < width; x++)
		{   // rd = normalize(P - ro)
            double Ro[3] = {0, 0, 0};
			double Rd[3] = {cx - (w/2) + pixelWidth * (x + 0.5), cy - (h/2) + pixelHeight * (y + 0.5), 1};
			normalize(Rd);
            double best_t = INFINITY;

            //printf("%d", Rd);
            // Go through each of the objects at each pixel and find out if they will intersect
			for (i = 0; i < numObjects; i++)
            {
				t = 0;
				if(strcmp((objects[i].type), "sphere")){

					t = sphereIntersection(objects[i].properties.sphere.position, objects[i].properties.sphere.radius,Rd, Ro);

				} else if(strcmp((objects[i].type), "plane")){

					t = planeIntersection(objects[i].properties.plane.position, objects[i].properties.plane.normal, Rd, Ro);

				}
				else
                {
					fprintf(stderr, "Error, object type is not valid.\n");
					exit(-1);
				}

				if (t > 0 && t < best_t)
                {
					best_t = t;
                    //printf("%d", t);
				}

				if(best_t > 0 && best_t != INFINITY)
                {
                    if(strcmp((objects[i].type), "sphere"))
                    {
                        //printf("Hey we are changing the color!");
                        pixel.r = objects[i].properties.sphere.color[0];
                        pixel.g = objects[i].properties.sphere.color[1];
                        pixel.b = objects[i].properties.sphere.color[2];
                        buffer->image[y * width + x] = pixel;
                    }
                    else if(strcmp((objects[i].type), "plane"))
                    {   //printf("plane color");
                        pixel.r = objects[i].properties.plane.color[0];
                        pixel.g = objects[i].properties.plane.color[1];
                        pixel.b = objects[i].properties.plane.color[2];
                        buffer->image[y * width + x] = pixel;
                    }
                    else
                    {
                        fprintf(stderr, "ERROR, Object type has no color property");
                        exit(-1);
                    }
				}
				// If no intersection let white be the background
				else
                { // printf("no intersection!!\n");

					pixel.r =  0;
					pixel.g = 0;
					pixel.b = 0;
					buffer->image[y * width + x] = pixel;

				}

			}

		}
	}

	return buffer;
}

int main(int argc, char *argv[])
{
	FILE *json;
	int numObjects, i;
    Pixmap* picbuffer;
    picbuffer = (Pixmap *)malloc(sizeof(Pixmap));

	json = fopen("objects.json", "r");
	if(json == NULL)
    {
		fprintf(stderr, "Error: could not open file.\n");
		fclose(json);
		exit(-1);

	}
	else
	{
		numObjects = readScene(json, objects);
        picbuffer = rayCaster(objects, 200, 200, numObjects);

        int size = picbuffer->height * picbuffer->width;
        //printf("%d", size);

        ppmWriter(picbuffer, "test2.ppm", size , 3);

	}
    printf("Number of objects read from the JSON file: %d\n", numObjects);



	return 0;
}
