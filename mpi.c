#include"header.h"
//put the MPI_Comm_cart here.
//probably copy and paste from cart3d

	int num_proc, my_rank, cart_rank;
//	MPI_Comm comm3D; i think i put this in the header
	int coord[namber];
	int bc[namber];
	int reorder, nrows,ncols;
	int p,rank_source[number_dims],rank_dest[number_dims];
	int check, dims[namber];
	//refer to params file i think
	//automatically set all to 1

	MPI_Init(&argc, &argv);
	int number_dims=3;
	int namber=3;
	dims[0]=num_proc_x;
	dims[1]=num_proc_y;
	dims[2]=num_proc_z;	
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	for(int j=1; j<argc; j++){
		if(strcmp(argv[j], "-dimsx")==0){
			dims[0]=atoi(argv[j+1]);
		}
		if(strcmp(argv[j], "-dimsy")==0){
			dims[1]=atoi(argv[j+1]);
		}
		if(strcmp(argv[j],"-dimsz")==0){
			dims[2]=atoi(argv[j+1]);
		}
	}

//dimension checker
	check=dims[0]*dims[1]*dims[2];
	if(num_proc!=check&&my_rank==0){
		printf("dimensions don't match the number of processes\n");
		MPI_Abort(MPI_COMM_WORLD,0);
	}

	bc[0],bc[1],bc[2]=1;
	MPI_Cart_create(MPI_COMM_WORLD, number_dims, dims, bc, 1, &comm3D);
	
	MPI_Barrier(comm3D);
	//all wait
	MPI_Cart_coords(comm3D, cart_rank, number_dims, coord);
	//gives cartesian coordinates to mpi group

	MPI_Barrier(comm3D); //so they all have same time
	MPI_Wtime(); //gets the time

//i dont think we want to run this yet?? because we still need to use mpi; ask for clarification about what ordder the links run in
//	 MPI_Finalize():
