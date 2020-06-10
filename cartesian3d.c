#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

int main(int argc, char *argv[]){

	int number_dims=3;
	int num_proc, my_rank, cart_rank;
	MPI_Comm comm3D;
	int namber=3;
	int dims[namber],coord[namber];
	int bc[namber];
	int reorder, nrows,ncols;
	int p,rank_source[number_dims],rank_dest[number_dims];

	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	nrows=ncols=3;
	for(int i=0; i<3;i++){
		dims[i]=2;
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
