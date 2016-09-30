#ifndef JSONOBJECTS_H_INCLUDED
#define JSONOBJECTS_H_INCLUDED

typedef struct Camera {
	double width;
	double height;
} Camera;


typedef struct Plane {
	double color[3];
	double position[3];
	double normal[3];
} Plane;


typedef struct Sphere {
	double color[3];
	double position[3];
	double radius;
} Sphere;

#endif // JSONOBJECTS_H_INCLUDED
