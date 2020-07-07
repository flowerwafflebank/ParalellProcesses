#include "header.h"
parameters->g.bank
grid
MPI_Status status;

void disperse(list of args)
{
	for(int i=1;i<length/dimension; i++){
		//here i will loop over the different grids?
		for(int j=1; j<=parameters->p; j++){
		//here I will loop over particles in each grid
		}
	}
}


/*
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)*/
