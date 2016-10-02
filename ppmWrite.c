// Created by Alejandro Varela
// Project 1: Images
#include <stdio.h>
#include <stdlib.h>

// This was added for project two in order to deal with the different colored objects
typedef struct PixelColor {
    unsigned char r, g, b;
} PixelColor;

// This will serve as my buffer
typedef struct Pixmap
{
    int width, height, magicNumber, color;
    PixelColor *image;
}Pixmap;

// Takes the buffer, output file, size of image, and the desired format
// Uses these in order to write the buffer out to the file in either P3 or P6
int ppmWriter(Pixmap *buffer, char *outputFileName, int size, int desiredFormat)
{
    FILE *destination;
    int i, j, numPix;
    char comment[] = {"#This was converted by Alejandro Varela"};
    char color[64];
    //printf("%s", comment);

    destination = fopen(outputFileName, "w");
    if (!destination)
    {
        fprintf(stderr,"\nERROR: Can't open the file for writing");
        fclose(destination);
        return -1;
    }
    else
    {
        fprintf(destination, "P%d\n%s\n%d %d\n%d\n", desiredFormat, comment, buffer->width, buffer->height, buffer->color);
        // Print out to the outfile in P6 format
        if(desiredFormat == 6)
        {
            numPix = fwrite(buffer->image, sizeof(PixelColor), size, destination);
        }
        // Print out to the outfile in P3 format
        else if(desiredFormat == 3)
       for(i = 0; i < (buffer->height); i++)
        {
			for(j = 0; j < (buffer->width); j++)
            {
				sprintf(color, "%d", buffer->image[(buffer->width) * i + j].r);
				fprintf(destination, "%s\n", color);
				sprintf(color, "%d", buffer->image[(buffer->width) * i + j].g);
				fprintf(destination, "%s\n", color);
				sprintf(color, "%d", buffer->image[(buffer->width) * i + j].b);
				fprintf(destination, "%s\n", color);
			}
		}
    }

    fclose(destination);
    return 0;
}

