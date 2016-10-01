#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
double planeIntersection(double* p, double* n, double* Ro, double* Rd){

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
    c = sqr(Ro[0] - p[0]) + sqr(Ro[1] - p[1]) + sqr(Ro[2] - p[2]) - sqr(r);
    det = sqr(b) - 4 * a * c;

    if(det < 0) return (-1);

    det = sqrt(det);
    t = (-b - det) / ( 2 * a);

    if(t > 0) return t;

    t = (-b + det) / (2 * a);

    if(t > 0) return (-1);

    return (-1);
}

