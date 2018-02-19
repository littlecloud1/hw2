Lai Man Tang
Yuan Qin

There are three files : mandelbrot_joe mandelbrot_susie mandelbrot_ms
to run this file, you have to edit the mandelbrot.sh

Examples:
eg for Joe in size 1000 * 1000
mpirun -np 1  ./mandelbrot_joe 1000 1000

eg for Susie in size 1000 * 1000
mpirun -np 1  ./mandelbrot_susie 1000 1000

eg for Master/Slave in size 1000 * 1000
mpirun   ./mandelbrot_ms 1000 1000