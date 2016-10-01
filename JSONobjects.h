// This is a header file which includes all of the object types that we will find int the
// JSON file. Each object has been separated but unionized in order to keep the
// Attributes separate but still be able to only have to allocate memory once when
// Creating the object array

#ifndef JSONOBJECTS_H_INCLUDED
#define JSONOBJECTS_H_INCLUDED

typedef struct Camera
{
	double width;
	double height;
} Camera;


typedef struct Plane
{
	double color[3];
	double position[3];
	double normal[3];

} Plane;


typedef struct Sphere
{
	double color[3];
	double position[3];
	double radius;

} Sphere;

typedef struct Object
{
	char *type;

	union properties
	{
		Camera camera;
		Plane plane;
		Sphere sphere;

	} properties;

} Object;

#endif // JSONOBJECTS_H_INCLUDED
