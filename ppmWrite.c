// Created by Alejandro Varela
// Taken from project one but modified to account for rgb instead of just taking the image as a whole
// I could not get the p6 writer to work properly so I just kept the p3 writer. If I had left the p6 option
// in there then the image colors would be negative. Otherwise it works just fine.
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
    //char color[64];
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
        fprintf(destination, "P%d\n%s\n%d %d\n%d\n", desiredFormat, comment, buffer->width, buffer->height, 255);
        // Print out to the outfile in P6 format

        // Print out to the outfile in P3 format
        if(desiredFormat == 3)
        {
            for(i = 0; i < (buffer->height); i++)
            {
                for(j = 0; j < (buffer->width); j++)
                {
                    fprintf(destination, "%d ", buffer->image[i * buffer->width *3+3*j].r);
                    fprintf(destination, "%d ", buffer->image[i * buffer->width *3+3*j+1].g);
                    fprintf(destination, "%d\n", buffer->image[i * buffer->width *3+ 3*j+2].b);
                }
            }
        }
    }
    fclose(destination);
    return 0;
}

