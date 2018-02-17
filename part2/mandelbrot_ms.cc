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
        master(P, height, width, minX, minY);
    }else{
        slave(width,minX,jt,minY,it);
        
    }
    
    MPI_Finalize();
    return 0;
}
void master(int P, int height, int width, double minX, double minY){
    
	double MStime = MPI_Wtime();

	int* resultBuff = (int*)malloc(height*weight*sizeof(int));

	MPI_Status status;
	int* recvdata = (int*)malloc((width+1) * sizeof(int));
	int  nextRow = 0;

	/*
	* Seed the slaves.
	*/
	for (int i = 1; i < P; i++) {
		MPI_Send(&nextRow, 1, MPI_INT, i, WORK_TAGm MPI_COMM_WORLD);
		nextRow++;
	}

	/*
	* Receive a result from any slave and dispatch a new work
	* request work requests have been exhausted.
	*/
	while (nexRow < height) {
		MPI_Recv(recvdata, width + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int p = status.MPI_SOURCE;
		MPI_Send(&nextRow, 1, MPI_INT, p, WORK_TAG, MPI_COMM_WORLD);
		int currRow = recvdata[width];
		memcpy(resultBuff + currRow * width, recvdata, width* sizeof(int));
		nextRow++;

	}

	/*
	* Receive results for outstanding work requests and Tell all the slaves to exit.
	*/

	for (int i = 1; i < P; i++) {
		MPI_Recv(recvdata, width + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int p = status.MPI_SOURCE;
		MPI_Send(&nextRow, 1, MPI_INT, p, DIE_TAG, MPI_COMM_WORLD);
		int currRow = recvdata[width];
		memcpy(resultBuff + currRow* width, recvdata, width* sizeof(int));
	
	}

	printf("MS's algorithm time: %lf :", MPI_Wtime() - MStime);

	//rendering image
	gil::rgb8_image_t img(height, width);
	auto img_view = gil::view(img);


	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			img_view(j, i) = render(resultBuff[i*width + j] / 512.0);
		}
	}



	gil::png_write_view("mandelbrot_susie.png", const_view(img));


}



void slave(int width, double minX, double jt, double minY, double it){
    int sentdata[width+1];
    MPI_Status status;
    int nextRow;
    while(true){
        MPI_Recv(&nextRow, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        if(status.MPI_TAG == DIE_TAG){return;}
        else{
            double y = minY +nextRow*it;
            for(int i =0,double x =minX;i<width; i++){
                sentdata[i] = mandelbrot(x,y);
                x = x+jt;
            }
            sentdata[width] = nextRow;
            MPI_Send(sentdata, width+1,MPI_INT,0,WORK_TAG,MPI_COMM_WORLD);
        }
    }
}

/* eof */
