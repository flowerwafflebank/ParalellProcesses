#include "header.h"

void disperse(Parameters *params, Geometry *geo, Bank *sb)
{
	for(int j=0; j<3;j++){//loops over dimensions
		if(params->num_proc[j]==1)
			continue;
		unsigned long ct=0,sen=0, numrec=0;
		unsigned long cx=sb->n;
		Particle *mh=(Particle*)malloc(sb->n*sizeof(Particle));
		sb->n=0;
		for(int i=0;i<cx; i++){//here i will loop over the particles in the source bank
			int c1=sb->p[i].coord[j]-params->coords[j]; //find difference between this and the coordinate in this rank
			if(c1)//preparing send buffer
				memcpy(&mh[ct++],&sb->p[i],sizeof(Particle));
			else 
				memcpy(&sb->p[sb->n++],&sb->p[i],sizeof(Particle));
		}
			MPI_Allreduce(&ct, &sen,1,MPI_UNSIGNED_LONG,MPI_MAX, params->my_comm);
			Particle *rb=(Particle*)malloc(sizeof(Particle)*sen);
			for(int k=0; k<params->num_proc[j]; k++){//here I will loop over particles in each grid
				MPI_Sendrecv(&ct, 1, MPI_UNSIGNED_LONG,geo->my_naybor[2*j+1],0,&numrec,1, MPI_UNSIGNED_LONG,geo->my_naybor[2*j],0,params->my_comm, MPI_STATUS_IGNORE); 
				MPI_Sendrecv(&mh, ct, params->my_mpi_data, geo->my_naybor[2*j+1],0, &rb, numrec, params->my_mpi_data, geo->my_naybor[2*j], 0, params->my_comm, MPI_STATUS_IGNORE);
				ct=0;
				for(int l=0; l<numrec; l++){
					if(rb[l].coord[j]!=params->coords[j])
						memcpy(&mh[ct++],&rb[l],sizeof(Particle));
					else{
						if(sb->n==sb->sz)
							sb->resize(sb);
						memcpy(&sb->p[sb->n++], &rb[l], sizeof(Particle));
					}
				}
			}
			//here we want to set up the sending algorithm to move particles from one to next process
		free(mh);
		free(rb);
	}
	return;
}


/*
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)*/
