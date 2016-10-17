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
#include "dotMath.h"

Object objects[128];

static inline double sqr(double v)
{
  return v*v;
}
static inline double v3_len(double* x)
{
    return sqrt(sqr(x[0]) + sqr(x[1]) + sqr(x[2]));
}

static inline void normalize(double* v)
{
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

static inline void v3_reflect(double* x, double* y,double* z)
{
    normalize(y);

    int perp = 2.0 * v3_dot(y, x);

     double* tmp;

    v3_scale(y, perp, tmp);
    //printf("\n testing perp scale");


    v3_subtract(x, tmp, z);

}


double clamp(double colorVal){
    double max = 1;
    double min = 0;

	if(colorVal > max) return max;

    else if (colorVal < min) return min;

    else   return colorVal;
}



// Takes in Ro - ray origin, Rd - ray direction, p - position of plane, and n - normal of the plane
// If there is no intersection with the plane then we will return -1 otherwise return a t
double planeIntersection(double* p, double* n, double* Rd, double* Ro){
    //double c = {p[0] - Ro[0], p[1] - Ro[1], p[2] - Ro[2]};
    normalize(n);
    double denom = v3_dot(n, Rd);
    double c[3];
    v3_subtract(p, Ro, c);
    //printf("working");exit (0);
    //double centMinRo = v3_subtract(p, Ro);

    if (fabs(denom) < 0.0001f) return -1;// using epsilon

    double t = v3_dot(c, n) / denom;

    // no intersection
    if(t < 0.0) return -1;

    return t;
}


// Takes in Ro - ray origin, Rd - ray direction, p - position or center sphere, and r - radius
// If at any point we find that there is no intersection with the sphere then we
// must make sure to return -1 that way we know whether or not to color in the pixel or not
double sphereIntersection(double* p, double r, double* Rd, double* Ro)
{

    double a, b, c;
    // calculate quadratic formula
    // First find a, b, c

    a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
    //printf("value of a: %lf\n", a);
    b = 2 * (Rd[0]*(Ro[0]-p[0]) + Rd[1]*(Ro[1]-p[1]) + Rd[2]*(Ro[2]-p[2]));
    //printf("value of b: %lf\n", b);
    c = sqr(Ro[0]-p[0]) + sqr(Ro[1]-p[1]) + sqr(Ro[2]-p[2]) - sqr(r);
    //printf("value of c: %lf\n", c);

    // make sure a is 1 (means the ray direction was normalized)
    if (a > 1.0001 || a < .9999) {

        //printf("a = %lf\n", a);
        fprintf(stderr, "Ray direction was not normalized\n");
        exit(-1);

    }

    // Check that discriminant is <, =, or > 0
    double disc = sqr(b) - 4*a*c;

    double t0, t1;  // solutions

    // No intersection
    if (disc < 0) return -1;

    // Single solution
    else if (disc == 0)
    {
        t0 = -1*(b / (2*a));
        return t0;
        //printf("t0 = %lf\n", t0);

    }

    // 2 solutions: find the smaller
    else
    {
        t0 = (-1*b - sqrt(sqr(b) - 4*c))/2;
        t1 = (-1*b + sqrt(sqr(b) - 4*c))/2;
        //printf("t0 = %lf\n", t0);
        //printf("t1 = %lf\n", t1);
    }

    // No intersection
    if (t0 < 0 && t1 < 0) return -1;


    else if (t0 < 0 && t1 > 0) return t1;

    else if (t0 > 0 && t1 < 0) return t0;

    // If they were both positives then return the smaller one
    else
    {

        if (t0 <= t1) return t0;

        else return t1;
    }
    //printf("value of t0: %lf\n", t0);
    //printf("value of t1: %lf\n", t1);
    //return -1;

}

int shadeCheck(Object objects[], double* newRd, double* newRo, int numObjects, int closestObject, double maxDistance)
{
     // Do intersections with new Ron and Rdn of all other objects in the scene(planes or spheres)
    int k;
    int newBest_o = 1;
    for(k = 0; k < numObjects; k++)
    {   // skip the closest object
        double newBestT = INFINITY;
        double newT = 0;
        if (k == closestObject) continue;
        else if(strcmp(objects[k].type, "sphere") == 0){
            newT = sphereIntersection(objects[k].properties.sphere.position, objects[k].properties.sphere.radius, newRd, newRo);

        } else if(strcmp(objects[k].type, "plane") == 0){
            newT = planeIntersection(objects[k].properties.plane.position, objects[k].properties.plane.normal, newRd, newRo);
        }
        // Shade the pixel because there is another object in the way of the light source

        if (newT > 0 && newT < newBestT)
        {
            newBestT = newT;
            //best_i = i;

        }
        if(newBestT > 0 && newBestT != INFINITY) return -1;

    }
    // found no intersections
    return 0;

}

// Find the diffuse for the closest object
// takes in the object normal, light vector, light color, and the objects diffuse color
// places the result into the outColor
void calculateDiffuse(double *objNormal, double *light, double *illumColor, double *objDiffuse, double *outColor) {

    // K_a*I_a should be added to the beginning of this whole thing, which is a constant and ambient light
    double normDotLight = v3_dot(objNormal, light);

    if (normDotLight > 0)
    {
        double diffuseProduct[3];

        diffuseProduct[0] = objDiffuse[0] * illumColor[0];
        diffuseProduct[1] = objDiffuse[1] * illumColor[1];
        diffuseProduct[2] = objDiffuse[2] * illumColor[2];

        // multiply by n_dot_l and store in out_color
        v3_scale(diffuseProduct, normDotLight, outColor);
        /*printf("Color is: %lf", outColor[1]);
        printf("Color is: %lf", outColor[2]);*/
    }

    else
    {
        // would normally return K_a*I_a here...
        outColor[0] = 0;
        outColor[1] = 0;
        outColor[2] = 0;

    }
}

// Find the specular for the closest object
// Takes in the ns which as shown in class we can set to one for testing at least
// light vector, light reflection, object normal, ray direction(V), object specular, light color
// Returns the proper color in the out color
void calculateSpecular(double ns, double *light, double *lightRef, double *objNormal,
                        double *V, double *objSpecular, double *illumColor, double *outColor) {

    double rayDotLight = v3_dot(V, lightRef);
    double normDotLight = v3_dot(objNormal, light);

    if (rayDotLight > 0 && normDotLight > 0)
    {
        double rayPowerOfNS = pow(rayDotLight, ns);

        double specularProduct[3];
        specularProduct[0] = objSpecular[0] * illumColor[0];
        specularProduct[1] = objSpecular[1] * illumColor[1];
        specularProduct[2] = objSpecular[2] * illumColor[2];

        v3_scale(specularProduct, rayPowerOfNS, outColor);

    }

    else
    {
        outColor[0] = 0;
        outColor[1] = 0;
        outColor[2] = 0;
    }

}

double calculateAngularAt(Object objects[], double rayToObject[3], int numObjects, int currLight)
{
    //TODO: check if light is not a spotlight. if not, return 1.0
            double directionDotLight = v3_dot(objects[currLight].properties.light.direction, rayToObject);
            double fang = pow(directionDotLight, objects[currLight].properties.light.angularAZero);
            return fang;

}
// Raycaster will take in the array of objects from the JSON file, the width and height desired
// as well as the actually number of objects that were read in in order to loop through the array
// of objects
// Will return the buffer holding the image that we want to print out to the ppm file
int rayCaster(Object objects[], Pixmap * buffer, double width, double height, int numObjects)
{
	double cx, cy, h, w, pixelHeight, pixelWidth;
	int i, x, y;
    double Ro[3] = {0, 0, 0};
	double Rd[3] = {0, 0, 0};
	double point[3] = {0,0, 1};
	double view[2] = {0,0};

	cx = 0;
	cy = 0;

	buffer->width = width;
	buffer->height = height;
	buffer->color = 255;

    //Grab the size of the view plane
	for(i = 0; i < numObjects; i++)
    {
		if(strcmp(objects[i].type, "camera") == 0)
        {
			h = objects[i].properties.camera.height;
			w = objects[i].properties.camera.width;
		}
	}

	pixelHeight = h / height;
    pixelWidth = w / width;

    // Start going row by column and checking whether each pixel intersects with an object or not
	for (y = 0; y < width; y++)
    {
        point[1] = -(view[1] - h/2.0 + pixelHeight*(y + 0.5));
		for (x = 0; x < height; x++)
		{
            point[0] = view[0] - w/2.0 + pixelWidth*(x + 0.5);
            normalize(point);
			Rd[0] = point[0];
            Rd[1] = point[1];
            Rd[2] = point[2];
            //normalize(Rd);

            double best_t = INFINITY;

            //printf("%d", Rd);
            // Go through each of the objects at each pixel and find out if they will intersect
			int best_i = 0;
			for (i = 0; i < numObjects; i++)
            {
				double t = 0;
				if(strcmp(objects[i].type, "sphere") == 0){

					t = sphereIntersection(objects[i].properties.sphere.position, objects[i].properties.sphere.radius,Rd, Ro);
                    //printf("This is i:%d", i);
				} else if(strcmp(objects[i].type, "plane") == 0){

					t = planeIntersection(objects[i].properties.plane.position, objects[i].properties.plane.normal, Rd, Ro);
                    //printf("%d", i);

				}

                // If they do set the t as well as the object that is intersecting
				if (t > 0 && t < best_t)
                {
					best_t = t;
                    best_i = i;
                    //printf("%d", t);
				}
				int l,k;
                double* color = malloc(sizeof(double)*3);


                // place the color of the current intersection into the image buffer
				if(best_t > 0 && best_t != INFINITY)
                {
                    if(strcmp(objects[best_i].type, "sphere") == 0)
                    {  // printf("hello there!");
                        color[0] = 0;//ambientColor[0];
                        color[1] = 0;//ambientColor[1];
                        color[2] = 0;//ambientColor[2];
                        // In order to find a shadow
                        for (l = 0; l < numObjects; l++)
                        {
                            // Look for a light to see if that object has a shadow casted on it by a light
                            if(strcmp(objects[l].type, "light") == 0)
                            {   // calc new ray origin and direction
                                double temp[3];
                                double Ron[3];
                                double Rdn[3];
                                v3_scale(Rd, best_t, temp);
                                v3_add(temp, Ro, Ron);
                                v3_subtract(objects[l].properties.light.position, Ron, Rdn);

                                double distanceTLight = v3_len(Rdn);
                                double otherObjIntersect = shadeCheck(objects, Rdn, Ron, numObjects, best_i, distanceTLight);
                                // there is an object in the way so shade in
                                if(otherObjIntersect == -1)
                                {   //printf("hi there");
                                    buffer->image[y*3 * buffer->width + x*3].r = 0*255 ;
                                    buffer->image[y*3 * buffer->width + x*3+1].g = 0*255;
                                    buffer->image[y*3 * buffer->width + x*3+2].b = 0*255;
                                    break;
                                }
                                // otherwise there is no object in between us and the light source
                                else
                                {
                                    //printf("hello");
                                    double spherePos[3] = {objects[best_i].properties.sphere.position[0],objects[best_i].properties.sphere.position[1],objects[best_i].properties.sphere.position[2]};

                                    double n[3] = {Ron[0] - spherePos[0], Ron[1]-spherePos[1], Ron[2]-spherePos[2]}; //objects[best_i]->normal;
                                    //n = v3_subtract(Ron, objects[best_i].properties.sphere.position, n);
                                    double lVector[3] = {Rdn[0], Rdn[1], Rdn[2]}; // light vector is just Rdn
                                    double lReflection[3];
                                    double V[3] = {Rd[0], Rd[1], Rd[2]};
                                    double diffuseColor[3];
                                    double specularColor[3];
                                    double diffPlusSpec[3];
                                    double lightRayToClosestObj[3];
                                    v3_reflect(lVector, n, lReflection); // in the book

                                    //printf("%d", objects[best_i].properties.sphere.diffuseColor[0]);
                                    calculateDiffuse(n, lVector, objects[l].properties.light.color, objects[best_i].properties.sphere.diffuseColor, diffuseColor);
                                    calculateSpecular(1, lVector, lReflection, n, V, objects[best_i].properties.sphere.specularColor, objects[l].properties.light.color, specularColor);

                                    v3_add(diffuseColor, specularColor, diffPlusSpec);

                                    v3_scale(lVector, -1, lightRayToClosestObj);

                                    double fang = calculateAngularAt(objects, lightRayToClosestObj, numObjects, l);
                                    //printf("Color is: %lf", diffuseColor[1]);
                                    //printf("Color is: %lf", diffPlusSpec[2]);

                                    // still need frad * rest
                                    color[0] += fang * diffPlusSpec[0];
                                    color[1] += fang * diffPlusSpec[1];
                                    color[2] += fang * diffPlusSpec[2];

                                    buffer->image[y*3 * buffer->width + x*3].r = clamp(color[0]) *255;
                                    buffer->image[y*3 * buffer->width + x*3+1].g = clamp(color[1]) *255;
                                    buffer->image[y*3 * buffer->width + x*3+2].b = clamp(color[2]) *255;
                                }

                            }

                        }

                        //printf("Hey we are changing the color!");
                        //1exit(0);
                        // Change object color for simply color

                        /*buffer->image[y*3 * buffer->width + x*3].r = 1 *255;
                        buffer->image[y*3 * buffer->width + x*3+1].g = 0 *255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = 0 *255;*/
                        //printf("hi");
                    }
                    else if(strcmp(objects[best_i].type, "plane") == 0)
                    {   //printf("plane color");
                        /*buffer->image[y*3 * buffer->width + x*3].r = objects[best_i].properties.plane.diffuseColor[0]*255 ;
                        buffer->image[y*3 * buffer->width + x*3+1].g = objects[best_i].properties.plane.diffuseColor[1]*255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = objects[best_i].properties.plane.diffuseColor[2]*255;*/
                        buffer->image[y*3 * buffer->width + x*3].r = 0 *255;
                        buffer->image[y*3 * buffer->width + x*3+1].g = 1 *255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = 0 *255;
                    }

				}
				// If no intersection let black be the background
				else
                { // printf("no intersection!!\n");

                    buffer->image[y*3 * buffer->width + x*3].r = 1*255 ;
                    buffer->image[y*3 * buffer->width + x*3+1].g = 1*255;
                    buffer->image[y*3 * buffer->width + x*3+2].b = 1*255;

				}

			}

		}
	}

	return 0;
}

int main(int argc, char *argv[])
{
	FILE *json;
	int numObjects, i;
	int ppmFormat = 3;
	double width, height;
	width = atof(argv[1]);
	height = atof(argv[2]);

	// Create the image buffer that will be used for the raycaster
	// malloc enough space for the image to account for the size and rgb
    Pixmap picbuffer;
    picbuffer.image = (PixelColor*)malloc(sizeof(PixelColor)*width* (height*3));

	json = fopen(argv[4], "r");
	if(json == NULL)
    {
		fprintf(stderr, "Error: could not open file.\n");
		fclose(json);
		exit(-1);

	}
	// Place the objects into an array, and then call the raycaster.
	// Finally take the image write out the file, and close the input file and free memory.
	else
	{
		numObjects = readScene(json, objects);
        rayCaster(objects, &picbuffer, width, height, numObjects);

        int size = height * width;
        //printf("%d", size);

        ppmWriter(&picbuffer, argv[3], size , ppmFormat);

	}
    //printf("Number of objects read from the JSON file: %d\n", numObjects);
    fclose(json);
    free(picbuffer.image);
    printf("Json scene has been created!");


	return 0;
}
