#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <JSONobjects.h>

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


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
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
    double* v = malloc(3*sizeof(double));

    expect_c(json, '[');
    skipWhiteSpace(json);
    v[0] = nextNumber(json);
    skipWhiteSpace(json);
    expect_c(json, ',');
    skipWhiteSpace(json);
    v[1] = nextNumber(json);
    skipWhiteSpace(json);
    expect_c(json, ',');
    skipWhiteSpace(json);
    v[2] = nextNumber(json);
    skipWhiteSpace(json);
    expect_c(json, ']');
    return v;
}


void read_scene(char* filename) {
  int c;
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }

  skipWhiteSpace(json);

  // Find the beginning of the list
  expect_c(json, '[');

  skipWhiteSpace(json);

  // Find the objects

  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skipWhiteSpace(json);

      // Parse the object
      char* key = nextString(json);
      if (strcmp(key, "type") != 0) {
	fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", lineNumber);
	exit(1);
      }

      skipWhiteSpace(json);

      expect_c(json, ':');

      skipWhiteSpace(json);

      char* value = nextString(json);

      if (strcmp(value, "camera") == 0) {
      } else if (strcmp(value, "sphere") == 0) {
      } else if (strcmp(value, "plane") == 0) {
      } else {
	fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, lineNumber);
	exit(1);
      }

      skipWhiteSpace(json);

      while (1) {
	// , }
	c = nextChar(json);
	if (c == '}') {
	  // stop parsing this object
	  break;
	} else if (c == ',') {
	  // read another field
	  skipWhiteSpace(json);
	  char* key = nextString(json);
	  skipWhiteSpace(json);
	  expect_c(json, ':');
	  skipWhiteSpace(json);
	  if ((strcmp(key, "width") == 0) ||
	      (strcmp(key, "height") == 0) ||
	      (strcmp(key, "radius") == 0)) {
	    double value = nextNumber(json);
	  } else if ((strcmp(key, "color") == 0) ||
		     (strcmp(key, "position") == 0) ||
		     (strcmp(key, "normal") == 0)) {
	    double* value = nextVector(json);
	  } else {
	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
		    key, lineNumber);
	    //char* value = next_string(json);
	  }
	  skipWhiteSpace(json);
	} else {
	  fprintf(stderr, "Error: Unexpected value on line %d\n", lineNumber);
	  exit(1);
	}
      }
      skipWhiteSpace(json);
      c = nextChar(json);
      if (c == ',') {
	// noop
	skipWhiteSpace(json);
      } else if (c == ']') {
	fclose(json);
	return;
      } else {
	fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", lineNumber);
	exit(1);
      }
    }
  }
}

int main(int c, char** argv) {
  read_scene(argv[1]);
  return 0;
}
