/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#define DIE_TAG 2
#define WORK_TAG 1
void slave(int width, double minX, double jt, double minY, double it);

int
main(int argc, char* argv[])
{
    double minX = -2.1;
    double maxX = 0.7;
    double minY = -1.25;
    double maxY = 1.25;
    
    int height, width;
    if (argc == 3) {
        height = atoi(argv[1]);
        width = atoi(argv[2]);
        assert(height > 0 && width > 0);
    }
    else {
        fprintf(stderr, "usage: %s <height> <width>\n", argv[0]);
        fprintf(stderr, "where <height> and <width> are the dimensions of the image.\n");
        return -1;
    }
    
    double it = (maxY - minY) / height;
    double jt = (maxX - minX) / width;
    double x, y;
    
    
    /*MPI section*/
    int p, P；
    double *recvdata;
    //double *sentdata ;
    
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    
    N = height / P;
    double Stime = MPI_Wtime();
    
    if (p == 0) {
        master();
    }else{
        slave(width,minX,jt,minY,it);
        
    }
    
    MPI_Finalize();
    return 0;
}
void master(){
    s
}
void slave(int width, double minX, double jt, double minY, double it){
    int sentdata[width+1];
    MPI_Status status;
    int nextrow;
    while(true){
        MPI_Recv(&nextrow, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        if(status.MPI_TAG == DIE_TAG){return;}
        else{
            double y = minY +nextrow*it;
            for(int i =0,double x =minX;i<width; i++){
                sentdata[i] = mandelbrot(x,y);
                x = x+jt;
            }
            sentdata[width] = nextrow;
            MPI_Send(sentdata, width+1,MPI_INT,0,WORK_TAG,MPI_COMM_WORLD);
        }
    }
}

/* eof */
