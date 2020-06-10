#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>


int main(int argc, char** argv){
	//this is something to do with initialization
	int p, my_rank, cart_rank;
	MPI_Comm comm2D; //not entirely sure what this means
	int number_dims=2;
	int dims[3], coordinates[3];
	int boundary_conditions[3];
	int reorder, nrows, ncols, ierr;
	//declaration of various things

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	nrows=ncols=2;
	dims[0]=dims[1]=dims[2]=0;
	boundary_conditions[0]=boundary_conditions[1]=boundary_conditions[3]=0;//i think these are BC
	reorder=1;//this allows reordering
	ierr=0;
	MPI_Dims_create(p,number_dims, dims);//topology	
	//this is the mapping?
	ierr=MPI_Cart_create(MPI_COMM_WORLD, number_dims, dims, boundary_conditions, reorder, &comm2D);	
	if(ierr!=0){printf("error with creation\n");}
	
	MPI_Cart_coords(comm2D, my_rank, number_dims, coordinates); 
		//the above does coordinates?
	MPI_Cart_rank(comm2D, coordinates, &cart_rank);
		//this is rank 
	printf("my rank: %d; cart rank: %d; my coordinates are: (%d,%d)\n",my_rank, cart_rank, coordinates[0],coordinates[1]);
		//print statement

//	printf("p:  %d",p); p is the number of processors
	MPI_Finalize();
	}
