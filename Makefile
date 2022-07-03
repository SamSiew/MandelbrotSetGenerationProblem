ALL: mandelbrot_serial mandelbrot_distributed

mandelbrot_serial : mandelbrot_serial.c
	gcc mandelbrot_serial.c -o mandelbrot_serial -lm

mandelbrot_distributed : mandelbrot_distributed.c
	mpicc mandelbrot_distributed.c -o mandelbrot_distributed  -lm

run:
	./mandelbrot_serial
	mpirun -np 4 mandelbrot_distributed

clean :
	/bin/rm -f mandelbrot_serial *.o
	/bin/rm -f mandelbrot_distributed *.o
