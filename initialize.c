#include "header.h"

Parameters *init_parameters(void)
{
  Parameters *p = malloc(sizeof(Parameters));

  p->n_particles = 1000000;
  p->n_batches = 10;
  p->n_generations = 1;
  p->n_active = 10;
  p->bc = REFLECT;
  p->n_nuclides = 1;
  p->tally = TRUE;
  p->n_bins = 16;
  p->seed = 1;
  p->nu = 2.5;
  p->xs_f = 0.012;
  p->xs_a = 0.03;
  p->xs_s = 0.27;
  p->gx = 400;
  p->gy = 400;
  p->gz = 400;
  p->write_tally = FALSE;
  p->write_keff = FALSE;
  p->tally_file = NULL;
  p->keff_file = NULL;
	p->dead=0;
	MPI_Datatype atom, base[2]={MPI_INT, MPI_DOUBLE};
	int blocks[2]={6,11};
	MPI_Aint offset[2],lb,extent;
	MPI_Type_get_extent(MPI_INT, &lb, &extent);
	offset[0]=lb; //how many bytes into the particle structure do the integers
	offset[1]=blocks[0]*extent+lb;//how many bytes in do the doubles begin
	MPI_Type_create_struct(2,blocks,offset, base, &atom);
	MPI_Type_commit(&atom);
	p->my_mpi_data=atom;
  MPI_Comm_size(MPI_COMM_WORLD, &p->size);
	MPI_Dims_create(p->size, 3,p->num_proc);
  return p;
}

Geometry *init_geometry(Parameters *parameters)
{
  Geometry *g = malloc(sizeof(Geometry));

  g->x = parameters->gx;
  g->xl=g->x/parameters->num_proc[0];
  g->y = parameters->gy;
	g->yl=g->y/parameters->num_proc[1];
  g->z = parameters->gz;
	g->zl=g->z/parameters->num_proc[2];
	for(int j=0; j<3;j++)
		MPI_Cart_shift(parameters->my_comm,j,1, &g->my_naybor[2*j], &g->my_naybor[2*j+1]); 
  g->bc = parameters->bc;

  return g;
}

Tally *init_tally(Parameters *parameters)
{
  Tally *t = malloc(sizeof(Tally));

  t->tallies_on = FALSE;
  t->nx = parameters->n_bins/parameters->num_proc[0];
  t->ny = parameters->n_bins/parameters->num_proc[1];
  t->nz = parameters->n_bins/parameters->num_proc[2];
  t->dx = parameters->gx/parameters->n_bins;
  t->dy = parameters->gy/parameters->n_bins;
  t->dz = parameters->gz/parameters->n_bins;
  t->flux = calloc(t->nx*t->ny*t->nz, sizeof(double));

  return t;
}

Material *init_material(Parameters *parameters)
{
  int i;
  Nuclide sum = {0, 0, 0, 0, 0};

  // Hardwire the material macroscopic cross sections for now to produce a keff
  // close to 1 (fission, absorption, scattering, total, atomic density)
  Nuclide macro = {parameters->xs_f, parameters->xs_a, parameters->xs_s,
     parameters->xs_f + parameters->xs_a + parameters->xs_s, 1.0};

  Material *m = malloc(sizeof(Material));
  m->n_nuclides = parameters->n_nuclides;
  m->nuclides = malloc(m->n_nuclides*sizeof(Nuclide));

  // Generate some arbitrary microscopic cross section values and atomic
  // densities for each nuclide in the material such that the total macroscopic
  // cross sections evaluate to what is hardwired above
  for(i=0; i<m->n_nuclides; i++){
    if(i<m->n_nuclides-1){
      m->nuclides[i].atom_density = rn()*macro.atom_density;
      macro.atom_density -= m->nuclides[i].atom_density;
    }
    else{
      m->nuclides[i].atom_density = macro.atom_density;
    }
    m->nuclides[i].xs_a = rn();
    sum.xs_a += m->nuclides[i].xs_a * m->nuclides[i].atom_density;
    m->nuclides[i].xs_f = rn();
    sum.xs_f += m->nuclides[i].xs_f * m->nuclides[i].atom_density;
    m->nuclides[i].xs_s = rn();
    sum.xs_s += m->nuclides[i].xs_s * m->nuclides[i].atom_density;
  }
  for(i=0; i<m->n_nuclides; i++){
    m->nuclides[i].xs_a /= sum.xs_a/macro.xs_a;
    m->nuclides[i].xs_f /= sum.xs_f/macro.xs_f;
    m->nuclides[i].xs_s /= sum.xs_s/macro.xs_s;
    m->nuclides[i].xs_t = m->nuclides[i].xs_a + m->nuclides[i].xs_s;
  }

  m->xs_f = parameters->xs_f;
  m->xs_a = parameters->xs_a;
  m->xs_s = parameters->xs_s;
  m->xs_t = parameters->xs_a + parameters->xs_s;

  return m;
}

Bank *init_source_bank(Parameters *parameters, Geometry *geometry)
{
  unsigned long i_p; // index over particles
  Bank *source_bank;

  // Initialize source bank
	if(parameters->rank==0)
 	 source_bank = init_bank(parameters->n_particles);
	else
 	 source_bank = init_bank(parameters->n_particles/parameters->size);

  // Sample source particles
  for(i_p=0; i_p<parameters->n_particles/parameters->size; i_p++){
    sample_source_particle(geometry, &(source_bank->p[i_p]));
    source_bank->n++;
  }

  return source_bank;
}

Bank *init_fission_bank(Parameters *parameters)
{
  Bank *fission_bank;
 if(parameters->rank==0)
	 fission_bank = init_bank(2*parameters->n_particles);
 else
	 fission_bank = init_bank(2*parameters->n_particles/parameters->size);
  return fission_bank;
}

Bank *init_bank(unsigned long n_particles)
{
  Bank *b = malloc(sizeof(Bank));
  b->p = malloc(n_particles*sizeof(Particle));
  b->sz = n_particles;
  b->n = 0;
  b->resize = resize_particles;

  return b;
}

void sample_source_particle(Geometry *geometry, Particle *p)
{
  p->alive = TRUE;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = rn()*geometry->x;
  p->y = rn()*geometry->y;
  p->z = rn()*geometry->z;
  p->hit=FALSE;
	p->coord[0]=p->x/geometry->xl;
	p->coord[1]=p->y/geometry->yl;
	p->coord[2]=p->z/geometry->zl;
	p->lx=p->x-p->coord[0]*geometry->xl;
	p->ly=p->y-p->coord[1]*geometry->yl;
	p->lz=p->z-p->coord[2]*geometry->zl;
  return;
}

void sample_fission_particle(Particle *p, Particle *p_old)
{
  p->alive = TRUE;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = p_old->x;
  p->y = p_old->y;
  p->z = p_old->z;
	p->hit=FALSE;

  return;
}

void resize_particles(Bank *b)
{
  b->p = realloc(b->p, sizeof(Particle)*2*b->sz);
  b->sz = 2*b->sz;

  return;
}

void free_bank(Bank *b)
{
  free(b->p);
  b->p = NULL;
  free(b);
  b = NULL;

  return;
}

void free_material(Material *material)
{
  free(material->nuclides);
  material->nuclides = NULL;
  free(material);
  material = NULL;

  return;
}

void free_tally(Tally *tally)
{
  free(tally->flux);
  tally->flux = NULL;
  free(tally);
  tally = NULL;

  return;
}
/*
MPI_DOUBLE doobie;
//this is me trying to declare an MPI_datatype not sure bout it

MPI_Aint holder;
//holding the mem address

MPI_Cart_rank(comm3D, coord, &cart_rank);

MPI_Cart_shift(comm3D, 0,1,&rank_source[0], &rank_dest[0]);
//int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest)    one needed for each dimension
MPI_Cart_shift(comm3D, 1,1, &rank_source[1], &rank_dest[1]);
MPI_Cart_shift(comm3D, 2,1, &rank_source[2], &rank_dest[2]);

MPI_Type_get_extent // gets and lower bound and extent of a datatype, which is basically what the sizeof() operator returns
MPI_Type_struct // once you have the types and count of the members within your new MPI datatype, pass them to this function to define the MPI datatype's structure
MPI_Type_commit // commit the new datatype to the MPI runtime
MPI_Cart_shift // call this three times, once per dimension, to store the neighbors in each direction
*/
