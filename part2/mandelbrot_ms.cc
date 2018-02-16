/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>

int
main (int argc, char* argv[])
{
    double minX = -2.1;
    double maxX = 0.7;
    double minY = -1.25;
    double maxY = 1.25;
    
    int height, width;
    if (argc == 3) {
        height = atoi (argv[1]);
        width = atoi (argv[2]);
        assert (height > 0 && width > 0);
    } else {
        fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
        fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
        return -1;
    }
    
    double it = (maxY - minY)/height;
    double jt = (maxX - minX)/width;
    double x, y;
    
    
    gil::rgb8_image_t img(height, width);
    auto img_view = gil::view(img);
    
  /* Lucky you, you get to write MPI code */
    int curr_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size (MPI_COMM_WORLD, &size);
    
    
}

/* eof */
