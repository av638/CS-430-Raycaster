all:	main.c
	gcc main.c -o raycaster

clean:
	rm -rf raycaster*~
