# CS-430-Raycaster
In this second project we are tasked with creating a simple raycaster.

The main has the actual ray casting and the intersection functions. Everything else helps with the parsing
and the writing. The only thing that gives you a slight warning is that I am returning a local variable when
doing the subtraction of vectors. I do this because I was having trouble saving a pointer to the vector that 
I got from the subtraction function. Otherwise it still works just fine.

USE: Command Line Args
width height outfile name(remember to use quotes!)

ex: 200 200 "hi.ppm"