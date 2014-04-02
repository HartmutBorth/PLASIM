most.x: most.c
	gcc -o most.x most.c -lm -L/usr/X11/lib -lX11
clean:
	rm -f *.o *.x F90* most_*
