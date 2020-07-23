#include "header.h"

// Main logic to move particle
void transport(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, Particle *p)
{
  double d_b;
  double d_c;
  double d;

  while(p->alive){

    // Find distance to boundary
    d_b = distance_to_boundary(geometry, p);

    // Find distance to collision
    d_c = distance_to_collision(material);

    // Take smaller of two distances
    d = d_b < d_c ? d_b : d_c;

    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;
		p->coord[0] = p->x / geometry->xl;
    p->coord[1] = p->y / geometry->yl;
    p->coord[2] = p->z / geometry->zl;
    p->lx = p->x - (p->coord[0] * geometry->xl);
    p->ly = p->y - (p->coord[1] * geometry->yl);
    p->lz = p->z - (p->coord[2] * geometry->zl);
//update the particles coordinate and the particles local x,y,z (same lines in init.c)
    // Case where particle crosses boundary
    if(d_b < d_c){
      cross_surface(geometry, p, parameters);
    }
    // Case where particle has collision
    else{
      collision(material, fission_bank, parameters->nu, p);
    }
		if(p->alive==FALSE)
			source_bank->dead++;
  }
  return;
}

// Returns the distance to the nearest boundary for a particle traveling in a
// certain direction
double distance_to_boundary(Geometry *geometry, Particle *p)
{
  int i;
  double dist;
  double d = D_INF;
  int    surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
  double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
  double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  double s_coords[6] = {0, geometry->x, 0, geometry->y, 0, geometry->z};
  
  for(i=0; i<6; i++){
    if(p_angles[i] == 0){
      dist = D_INF;
    }
    else{
      dist = (s_coords[i] - p_coords[i])/p_angles[i];
      if(dist <= 0){
        dist = D_INF;
      }
    }
    if(dist < d){
      d = dist;
      p->surface_crossed = surfaces[i];
    }
  }

  return d;
}

// Returns the distance to the next collision for a particle
double distance_to_collision(Material *material)
{
  double d;

  if(material->xs_t == 0){
    d = D_INF;
  }
  else{
    d = -log(rn())/material->xs_t;
  }

  return d;
}

// Handles a particle crossing a surface in the geometry
void cross_surface(Geometry *geometry, Particle *p, Parameters *parameters)
{
  // Handle vacuum boundary conditions (particle leaks out)
  if(geometry->bc == VACUUM){
    p->alive = FALSE;
  }
//update particle lx,y,z and coord appropriately
  // Handle reflective boundary conditions
  else if(geometry->bc == REFLECT){
    if(p->surface_crossed == X0){
      p->u = -p->u;
      p->x = 0.0;
			p->coord[0]=0;
			p->lx=0.0;
    }
    else if(p->surface_crossed == X1){
      p->u = -p->u;
      p->x = geometry->x;
			p->coord[0]=parameters->num_proc[0]-1;//not positive on this one
			p->lx=geometry->xl;
    }
    else if(p->surface_crossed == Y0){
      p->v = -p->v;
      p->y = 0.0;
			p->coord[1]=0;
			p->ly = 0.0; 
    }
    else if(p->surface_crossed == Y1){
      p->v = -p->v;
      p->y = geometry->y;
			p->coord[1]=parameters->num_proc[1]-1;
			p->ly=geometry->yl;
    }
    else if(p->surface_crossed == Z0){
      p->w = -p->w;
      p->z = 0.0;
			p->coord[2]=0;
			p->lz=0.0;
    }
    else if(p->surface_crossed == Z1){
      p->w = -p->w;
      p->z = geometry->z;
			p->coord[2]=parameters->num_proc[2]-1;
			p->lz=geometry->zl;
    }
  }
  
  // Handle periodic boundary conditions
  else if(geometry->bc == PERIODIC){
    if(p->surface_crossed == X0){
      p->x = geometry->x;
			p->lx=geometry->xl; ///update below accordingly
			p->coord[0]=parameters->num_proc[0]-1;
    }
    else if(p->surface_crossed == X1){
      p->x = 0;
			p->lx=geometry->xl;
			p->coord[0]=0; //originally i had +1 cuz we are going in the other direction, but changed it to 0
    }
    else if(p->surface_crossed == Y0){
      p->y = geometry->y;
			p->ly=geometry->yl;
			p->coord[1]=parameters->num_proc[1]-1;
    }
    else if(p->surface_crossed == Y1){
      p->y = 0;
			p->ly=geometry->yl;
			p->coord[1]=0;
    }
    else if(p->surface_crossed == Z0){
      p->z = geometry->z;
			p->lz=geometry->zl;
			p->coord[2]=parameters->num_proc[2]-1;
    }
    else if(p->surface_crossed == Z1){
      p->z = 0;
			p->lz=geometry->zl;
			p->coord[2]=0;
    }
  }

  return;
}

void collision(Material *material, Bank *fission_bank, double nu, Particle *p)
{
  int n;
  int i = 0;
  double prob = 0.0;
  double cutoff;
  Nuclide nuc = {0, 0, 0, 0, 0};

  // Cutoff for sampling nuclide
  cutoff = rn()*material->xs_t;

  // Sample which nuclide particle has collision with
  while(prob < cutoff){
    nuc = material->nuclides[i];
    prob += nuc.atom_density*nuc.xs_t;
    i++;
  }

  // Cutoff for sampling reaction
  cutoff = rn()*nuc.xs_t;

  // Sample fission
  if(nuc.xs_f > cutoff){

    // Sample number of fission neutrons produced
    if(rn() > nu - (int)nu){
      n = nu;
    }
    else{
      n = nu + 1;
    }

    // Sample n new particles from the source distribution but at the current
    // particle's location
    if(fission_bank->n+n >= fission_bank->sz){
      fission_bank->resize(fission_bank);
    }
    for(i=0; i<n; i++){
      sample_fission_particle(&(fission_bank->p[fission_bank->n]), p);
      fission_bank->n++;
    }
    p->alive = FALSE;
  }

  // Sample absorption (disappearance)
  else if(nuc.xs_a > cutoff){
    p->alive = FALSE;
  }

  // Sample scattering
  else{
    p->mu = rn()*2 - 1;
    p->phi = rn()*2*PI;
    p->u = p->mu;
    p->v = sqrt(1 - p->mu*p->mu) * cos(p->phi);
    p->w = sqrt(1 - p->mu*p->mu) * sin(p->phi);
  }
	p->hit=TRUE;
  return;
}
