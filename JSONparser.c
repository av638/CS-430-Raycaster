#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "JSONobjects.h"

// Line number to report errors
int lineNumber = 1;


// Checks to see if the current character in the file is a new line
// If it is then we want to increase the overall line count
int nextChar(FILE* json) {
    int c = fgetc(json);
    // If new line add one to the line number counter
    if (c == '\n')
    {
        lineNumber += 1;
    }

    // If end of file exit program and close the file
    if (c == EOF)
    {
        fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", lineNumber);
        fclose(json);
        exit(-1);
    }

    return c;
}


// Checks that the next character is what is expected
// If it is not then we will exit and return an error
void expect_c(FILE* json, int d) {
    int c = nextChar(json);
    if (c == d) return;

    fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, lineNumber);
    fclose(json);
    exit(-1);
}


// Checks the current character in the file and if it is a space then we will skip it
void skipWhiteSpace(FILE* json) {
    int c = nextChar(json);

    while (isspace(c))
    {
      c = nextChar(json);
    }
    ungetc(c, json);
}


// Validation which checks the string for escape codes, if it is larger than ascii
// and if it is larger than 128
char* nextString(FILE* json)
{
    char buffer[129];
    int c = nextChar(json), i = 0;

    if (c != '"')
    {
        fprintf(stderr, "Error: Expected string on line %d.\n", lineNumber);
        fclose(json);
        exit(-1);
    }
    c = nextChar(json);

    while (c != '"')
    {
        if (i >= 128)
        {
          fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
          fclose(json);
          exit(-1);
        }
        if (c == '\\')
        {
          fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
          fclose(json);
          exit(-1);
        }
        if (c < 32 || c > 126)
        {
          fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
          fclose(json);
          exit(-1);
        }
        // Add the character to the buffer and move on to the next character in the string
        buffer[i] = c;
        i += 1;
        c = nextChar(json);
    }
    // Null terminate the buffer and return the string for use
    buffer[i] = 0;
    return strdup(buffer);
}

// Checks to see that the current character in the file is in fact a double
double nextNumber(FILE* json)
{
    double value;
    if(fscanf(json, "%lf", &value) == 0)
    {
        fprintf(stderr, "Error, line number %d; expected numeric value.\n", lineNumber);
        fclose(json);
        exit(-1);
    }
    return value;
}

// Takes in an array containing a vector and returns a proper array of doubles
double* nextVector(FILE* json)
{
    double* vector = malloc(3*sizeof(double));

    expect_c(json, '[');
    skipWhiteSpace(json);
    vector[0] = nextNumber(json);

    skipWhiteSpace(json);
    expect_c(json, ',');
    skipWhiteSpace(json);
    vector[1] = nextNumber(json);

    skipWhiteSpace(json);
    expect_c(json, ',');
    skipWhiteSpace(json);
    vector[2] = nextNumber(json);

    skipWhiteSpace(json);
    expect_c(json, ']');

    return vector;
}

// Reads in a scene of objects taken from a JSON file
// And accounts for varience in the way that the JSON file is created
// EX: empty scenes, objects, comma separated and non-comma separated  objects/value pairs
int readScene(FILE *json, Object objects[])
{
    int c, numObjects;
	double *vector;
	char *name, *value;
	numObjects = 0;

	// Skip whitespace(s) read in the first character
	skipWhiteSpace(json);

	c = nextChar(json);

	// Check to see if the first character is an opening bracket
	if(c != '[')
    {
		fprintf(stderr, "Error, line number %d; invalid scene definition '%c'\n", lineNumber, c);
		fclose(json);
		exit(-1);
	}

	skipWhiteSpace(json);
	c = nextChar(json);

	// Check to see if the json file has no objects in it
	if(c != ']') ungetc(c, json);

	// While we have not run out of objects
	while(c != ']')
    {
		skipWhiteSpace(json);
		c = nextChar(json);

		// Check to see if this is the start of an object
		if(c != '{')
        {
			fprintf(stderr, "Error, line number %d; invalid object definition '%c'\n", lineNumber, c);
			fclose(json);
			exit(-1);
		}

		skipWhiteSpace(json);
		c = nextChar(json);

        // While still in the same object
		while(c != '}')
        {
			// If we encounter a string then unget and pass it into our nextString function
			if(c == '"') ungetc(c, json);

			name = nextString(json);

			if(strcmp(name, "type") == 0)
            {
				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);

				}
				// Grab the type of the object and place it into our object array
				else
                {
					skipWhiteSpace(json);
					value = nextString(json);
					objects[numObjects].type = value;
					// if we have a light then lets set all of the attributes like radial to 0 just to be safe
					// This is just in case the light has some but not all of the properties that we account for in the struct
					if(strcmp(value, "light") == 0)
                    {
                        objects[numObjects].properties.light.angularAZero = 0;
                        objects[numObjects].properties.light.radialAOne = 0;
                        objects[numObjects].properties.light.radialATwo = 0;
                        objects[numObjects].properties.light.radialAZero = 0;

                        objects[numObjects].properties.light.direction[0] = 0;
                        objects[numObjects].properties.light.direction[1] = 0;
                        objects[numObjects].properties.light.direction[2] = 0;

                        //objects[numObjects].properties.light.theta = NULL;

                    }
				}

			} else if(strcmp(name, "width") == 0) {

				skipWhiteSpace(json);
				//printf("%c", c);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-3);
				}
				// If it is the width then lets place it in the camera structure in the object array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.camera.width = nextNumber(json);
				}

			} else if(strcmp(name, "height") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the height then lets place it in the camera structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.camera.height = nextNumber(json);
				}

			} else if(strcmp(name, "radius") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
                // If it is the radius then lets place it in the sphere structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.sphere.radius = nextNumber(json);
				}

			} else if(strcmp(name, "radial-a2") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the radialATwo then lets place it in the camera structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.light.radialATwo = nextNumber(json);
				}

			} else if(strcmp(name, "radial-a1") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the radialAOne then lets place it in the camera structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.light.radialAOne = nextNumber(json);
				}

			} else if(strcmp(name, "radial-a0") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the radialAZero then lets place it in the camera structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.light.radialAZero = nextNumber(json);
				}

			} else if(strcmp(name, "angular-a0") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the angularAZero then lets place it in the camera structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.light.angularAZero = nextNumber(json);
				}

			} else if(strcmp(name, "theta") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the angularAZero then lets place it in the camera structure in the objects array
				else
                {
					skipWhiteSpace(json);
					objects[numObjects].properties.light.theta = nextNumber(json);
				}

			}else if(strcmp(name, "color") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If we have reached the color then lets place it into the right object in an array
				else
                {
					skipWhiteSpace(json);
					vector = nextVector(json);

                        if(strcmp(objects[numObjects].type, "light") == 0){
                        objects[numObjects].properties.light.color[0] = vector[0];
						objects[numObjects].properties.light.color[1] = vector[1];
						objects[numObjects].properties.light.color[2] = vector[2];
					}

				}

			}else if(strcmp(name, "diffuse_color") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If we have reached the diffuseColor then lets place it into the right object in an array
				else
                {
					skipWhiteSpace(json);
					vector = nextVector(json);
					if(strcmp(objects[numObjects].type, "plane") == 0) {

						objects[numObjects].properties.plane.diffuseColor[0] = vector[0];
						objects[numObjects].properties.plane.diffuseColor[1] = vector[1];
						objects[numObjects].properties.plane.diffuseColor[2] = vector[2];

					}else if(strcmp(objects[numObjects].type, "sphere") == 0) {

                    objects[numObjects].properties.sphere.diffuseColor[0] = vector[0];
                    objects[numObjects].properties.sphere.diffuseColor[1] = vector[1];
                    objects[numObjects].properties.sphere.diffuseColor[2] = vector[2];
					}
                }

            } else if(strcmp(name, "specular_color") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If we have reached the specularColor then lets place it into the right object in an array
				else
                {
					skipWhiteSpace(json);
					vector = nextVector(json);
					if(strcmp(objects[numObjects].type, "plane") == 0) {

						objects[numObjects].properties.plane.specularColor[0] = vector[0];
						objects[numObjects].properties.plane.specularColor[1] = vector[1];
						objects[numObjects].properties.plane.specularColor[2] = vector[2];

					}else if(strcmp(objects[numObjects].type, "sphere") == 0) {

                        objects[numObjects].properties.sphere.specularColor[0] = vector[0];
                        objects[numObjects].properties.sphere.specularColor[1] = vector[1];
                        objects[numObjects].properties.sphere.specularColor[2] = vector[2];
					}
                }

            }else if(strcmp(name, "position") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If we have reached the position then lets place it in an array of the proper object
				else
                {
					skipWhiteSpace(json);
					vector = nextVector(json);

					if(strcmp(objects[numObjects].type, "sphere") == 0)
                    {
						objects[numObjects].properties.sphere.position[0] = vector[0];
						objects[numObjects].properties.sphere.position[1] = vector[1];
						objects[numObjects].properties.sphere.position[2] = vector[2];

					} else if(strcmp(objects[numObjects].type, "plane") == 0) {
						objects[numObjects].properties.plane.position[0] = vector[0];
						objects[numObjects].properties.plane.position[1] = vector[1];
						objects[numObjects].properties.plane.position[2] = vector[2];

					} else if(strcmp(objects[numObjects].type, "light") == 0){
                        objects[numObjects].properties.light.position[0] = vector[0];
						objects[numObjects].properties.light.position[1] = vector[1];
						objects[numObjects].properties.light.position[2] = vector[2];
					}

				}

			} else if(strcmp(name, "normal") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; unexpected character '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it we have found the normal attribute then lets place it in the plane structure in the object array
				else
                {
					skipWhiteSpace(json);
					vector = nextVector(json);

					objects[numObjects].properties.plane.normal[0] = vector[0];
					objects[numObjects].properties.plane.normal[1] = vector[1];
					objects[numObjects].properties.plane.normal[2] = vector[2];
				}

			} else if(strcmp(name, "direction") == 0) {

				skipWhiteSpace(json);
				c = nextChar(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; unexpected character '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it we have found the normal attribute then lets place it in the plane structure in the object array
				else
                {
					skipWhiteSpace(json);
					vector = nextVector(json);

					objects[numObjects].properties.light.direction[0] = vector[0];
					objects[numObjects].properties.light.direction[1] = vector[1];
					objects[numObjects].properties.light.direction[2] = vector[2];
				}

			}

			// Otherwise if there was some other object attribute then lets exit the program
			else
            {

				fprintf(stderr, "Error, line number %d; invalid type '%s'.\n", name);
				fclose(json);
				exit(-1);
			}

			skipWhiteSpace(json);
			c = nextChar(json);

			if(c == ',')
            {
				skipWhiteSpace(json);
                c = nextChar(json);
			}
		}  // End of object parsing

		skipWhiteSpace(json);
		c = nextChar(json);

		if(c == '{') ungetc(c, json);

		if(c == ',')
        {
			skipWhiteSpace(json);
			c = nextChar(json);

			if(c == '{') ungetc(c, json);
		}

		// Increment array index counter
		numObjects += 1;

	} // End of JSON file

	// Return the total number of objects that were in the JSON file
	return numObjects;
}

