#include <stdio.h>
#include <mpi.h>
#include <time.h>

//I will now attempt to implenment MPI in a simple code

int main(int argc, char** argv){
	MPI_Init(NULL,NULL);
	//this is something to do with initialization
	int number_of_processors;
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
	//declares the number of processors
	int rank_of_processors;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_processors);
	//rank is the identification number of each processor
	char name_of_processor[MPI_MAX_PROCESSOR_NAME];
	int length_of_processor_name;
	MPI_Get_processor_name(name_of_processor, &length_of_processor_name);


	printf("Hello world from processor %s, rank %d\n", name_of_processor, rank_of_processors);


	clock_t start=clock();
	while(clock() < start+3000)
		;
//test comment!

	 
	int namber;
	if (rank_of_processors == 0) {
		namber=0;
		MPI_Send(&namber, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		printf("Rank %d: sent %d\n", rank_of_processors, namber);
	}
	if (rank_of_processors == 1) {
		MPI_Recv(&namber, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 	namber++;
		MPI_Send(&namber, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("Rank %d: sent %d\n", rank_of_processors, namber);	
	}
	
	
	if (rank_of_processors == 2) {
		MPI_Recv(&namber, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 	namber++;
		MPI_Send(&namber, 1, MPI_INT, 3, 0, MPI_COMM_WORLD);
		printf("Rank %d: sent %d\n", rank_of_processors, namber);
	}
	
	if (rank_of_processors == 3) {
		MPI_Recv(&namber, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 	namber++;
		MPI_Send(&namber, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("Rank %d: sent %d\n", rank_of_processors, namber);
	}
	
	MPI_Finalize(); //I read that this frees up memory. sorta like an fclose i guess
}
