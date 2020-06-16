#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

	int number_dims=3;
	int num_proc, my_rank, cart_rank;
	MPI_Comm comm3D;
	int namber=3;
	int coord[namber];
	int bc[namber];
	int reorder, nrows,ncols;
	int p,rank_source[number_dims],rank_dest[number_dims];
	int check, dims[namber];
	dims[0],dims[1],dims[2]=0;
   /* if(strcmp(arg, "-particles") == 0){
      if(++i < argc){
        long long n_particles = atoll(argv[i]);
        if(n_particles < 1)
          print_error("Number of particles must be greater than 0");
        parameters->n_particles = n_particles;
      }
      else print_error("Error reading command line input '-particles'");
    }*/

	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	for(int j=1; j<argc; j++){
		if(strcmp(argv[j], "x")==0){
			dims[0]=atoi(argv[j+1]);
		}
		if(strcmp(argv[j], "y")==0){
			dims[1]=atoi(argv[j+1]);
		}
		if(strcmp(argv[j],"z")==0){
			dims[2]=atoi(argv[j+1]);
		}
	}
//	printf("%d,%d,%d,check:%d,%s,%s,%s,\n",dims[0],dims[1],dims[2],check,argv[1],argv[2],argv[3]);
	check=dims[0]*dims[1]*dims[2];
	if(num_proc!=check&&my_rank==0){
		printf("dimensions don't match the number of processes\n");
		MPI_Abort(MPI_COMM_WORLD,0);
	}


	clock_t start=clock();
	while(clock() < start+3000)
		;


	nrows=ncols=3;
	for(int i=0; i<3;i++){
//		dims[i]=2;
		bc[i]=1;
	}

//	MPI_Dims_create(num_proc, number_dims, dims);	
	reorder=1;
	MPI_Cart_create(MPI_COMM_WORLD, number_dims, dims, bc, reorder, &comm3D);
//int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[], int reorder, MPI_Comm * comm_cart)
	MPI_Cart_rank(comm3D, coord, &cart_rank);
//
	MPI_Cart_coords(comm3D, cart_rank, number_dims, coord);
	//this appears to not be right for the moment ^^
	MPI_Cart_shift(comm3D, 0,1,&rank_source[0], &rank_dest[0]);
//int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)    one needed for each dimension
	MPI_Cart_shift(comm3D, 1,1, &rank_source[1], &rank_dest[1]);
	MPI_Cart_shift(comm3D, 2,1, &rank_source[2], &rank_dest[2]);
//	printf("source: %d, dest: %d \n", rank_source[1], rank_dest[1]);
//	printf("my cartesian rank: %d; my coordinates: (%d,%d,%d)\n",cart_rank, coord[0],coord[1],coord[2]);
	namber=my_rank*12-13;
	MPI_Send(&namber, 1, MPI_INT, rank_dest[1], 0, comm3D);
//int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
	int hold=namber;
	MPI_Recv(&namber, 1, MPI_INT, rank_source[1], 0, comm3D, MPI_STATUS_IGNORE);
//int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status * status)
	printf("\n%d sent %d to %d and recieved %d from %d\n", cart_rank, hold, rank_dest[1], namber, rank_source[1]);
	MPI_Finalize();


}
