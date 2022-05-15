/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_cac_buck.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "neigh_list.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "asa_user.h"

#define MAXNEIGHOUT  50
#define MAXNEIGHIN  10
#define EXPAND 10
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCACBuck::PairCACBuck(LAMMPS *lmp) : PairCAC(lmp)
{
  restartinfo = 0;
  nmax = 0;
  outer_neighflag = 0;
  flux_enable = 1;
}

/* ---------------------------------------------------------------------- */

PairCACBuck::~PairCACBuck() {
  if (allocated) {
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(cut);
  memory->destroy(a);
  memory->destroy(rho);
  memory->destroy(c);
  memory->destroy(rhoinv);
  memory->destroy(buck1);
  memory->destroy(buck2);
  memory->destroy(offset);
  memory->destroy(inner_neighbor_coords);
  memory->destroy(inner_neighbor_types);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCACBuck::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut_lj");
  memory->create(a, n + 1, n + 1, "pair:a");
  memory->create(rho, n + 1, n + 1, "pair:rho");
  memory->create(c, n + 1, n + 1, "pair:c");
  memory->create(rhoinv, n + 1, n + 1, "pair:rhoinv");
  memory->create(buck1, n + 1, n + 1, "pair:buck1");
  memory->create(buck2, n + 1, n + 1, "pair:buck2");
  memory->create(offset, n + 1, n + 1, "pair:offset");
  PairCAC::allocate();
}

/* ----------------------------------------------------------------------
global settings
------------------------------------------------------------------------- */
void PairCACBuck::settings(int narg, char **arg) {
  if (narg <1 || narg>2) error->all(FLERR, "Illegal pair_style cac/buck command");

  force->newton_pair = 0;
  cut_global_s = utils::numeric(FLERR, arg[0],false,lmp);

   if (narg == 2) {
    if (strcmp(arg[1], "one") == 0) atom->one_layer_flag=one_layer_flag = 1;
    else error->all(FLERR, "Unexpected argument in pair cac/buck invocation; only accepts cutoff and the 'one' keyword");
  }
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global_s;
  }

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACBuck::coeff(int narg, char **arg) {
  if (narg < 5 || narg > 6)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a_one = utils::numeric(FLERR, arg[2],false,lmp);
  double rho_one = utils::numeric(FLERR, arg[3],false,lmp);
  if (rho_one <= 0) error->all(FLERR, "Incorrect args for pair coefficients");
  double c_one = utils::numeric(FLERR, arg[4],false,lmp);

  double cut_one = cut_global_s;
  if (narg == 6) cut_one = utils::numeric(FLERR, arg[5],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      c[i][j] = c_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACBuck::init_one(int i, int j) {

  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  rhoinv[i][j] = 1.0 / rho[i][j];
  buck1[i][j] = a[i][j] / rho[i][j];
  buck2[i][j] = 6.0*c[i][j];

  if (offset_flag && (cut[i][j] > 0.0)) {
    double rexp = exp(-cut[i][j] / rho[i][j]);
    offset[i][j] = a[i][j] * rexp - c[i][j] / pow(cut[i][j], 6.0);
  }
  else offset[i][j] = 0.0;

  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double rho1 = rho[i][j];
    double rho2 = rho1*rho1;
    double rho3 = rho2*rho1;
    double rc = cut[i][j];
    double rc2 = rc*rc;
    double rc3 = rc2*rc;
    etail_ij = 2.0*MY_PI*all[0] * all[1] *
      (a[i][j] * exp(-rc / rho1)*rho1*(rc2 + 2.0*rho1*rc + 2.0*rho2) -
        c[i][j] / (3.0*rc3));
    ptail_ij = (-1 / 3.0)*2.0*MY_PI*all[0] * all[1] *
      (-a[i][j] * exp(-rc / rho1)*
      (rc3 + 3.0*rho1*rc2 + 6.0*rho2*rc + 6.0*rho3) + 2.0*c[i][j] / rc3);
  }

  return cut_global_s;
}

/* ---------------------------------------------------------------------- */


void PairCACBuck::init_style()
{
  PairCAC::init_style();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;

}

/* ---------------------------------------------------------------------- */


void PairCACBuck::force_densities(int iii, double s, double t, double w, double coefficients,
  double &force_densityx, double &force_densityy, double &force_densityz) {

  double delx,dely,delz;
  double shape_func;
  double shape_func2;
  int timestep=update->ntimestep;
  double *special_lj = force->special_lj;
  double  fpair;
  int *type = atom->type;
  double distancesq;
  double scan_position[3];
  double rcut;
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;

//scan the surrounding unit cell locations in a cartesian grid
//of isoparametric space until the cutoff is exceeded
//for each grid scan

  int distanceflag=0;
  rcut = cut_global_s;
  int origin_type = type_array[poly_counter];
  int listtype;
  int listindex;
  int poly_index;
  int scan_type;
  int element_index;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int neigh_max = inner_quad_lists_counts[pqi];
  int **node_types = atom->node_types;
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  double ****nodal_positions = atom->nodal_positions;
  
  //allocate arrays that store neighbor information around just this quadrature point
  allocate_quad_memory();
  //set virtual neighbor types, etc.
  init_quad_arrays();
  //interpolate virtual atom coordinates from shape functions corresponding to unit cells
  interpolation(iii,s,t,w);

  //compute force at quadrature point
  for (int l = 0; l < neigh_max; l++) {
    scan_type = inner_neighbor_types[l];
    scan_position[0] = inner_neighbor_coords[l][0];
    scan_position[1] = inner_neighbor_coords[l][1];
    scan_position[2] = inner_neighbor_coords[l][2];
    delx = current_position[0] - scan_position[0];
    dely = current_position[1] - scan_position[1];
    delz = current_position[2] - scan_position[2];
    distancesq = delx*delx + dely*dely + delz*delz;
    if (distancesq < cut[origin_type][scan_type]* cut[origin_type][scan_type]) {
      fpair = pair_interaction(distancesq, origin_type, scan_type);
      force_densityx += delx*fpair;
      force_densityy += dely*fpair;
      force_densityz += delz*fpair;
      if(atom->CAC_virial){
        virial_density[0] += 0.5*delx*delx*fpair;
        virial_density[1] += 0.5*dely*dely*fpair;
        virial_density[2] += 0.5*delz*delz*fpair;
        virial_density[3] += 0.5*delx*dely*fpair;
        virial_density[4] += 0.5*delx*delz*fpair;
        virial_density[5] += 0.5*dely*delz*fpair;
      }
      if (quad_eflag) {
        quadrature_energy += (a[origin_type][scan_type]*rexp - c[origin_type][scan_type]*r6inv -
        offset[origin_type][scan_type])/2;
      }
      //cac flux contribution due to current quadrature point and neighbor pair interactions
      if(quad_flux_flag){
        current_quad_flux(l,delx*fpair,dely*fpair,delz*fpair);
      }
    }
  }
  //end of force density loop
  
  //additional cac flux contributions due to neighbors interacting with neighbors
  //  in the vicinity of this quadrature point
  if (quad_flux_flag) {
    //compute_intersections();
    quad_neigh_flux();
  }
}

/* ---------------------------------------------------------------------- */

double PairCACBuck::pair_interaction(double distancesq, int itype, int jtype) {
  double fpair, forcebuck;
  double r2inv, r;

  r2inv = 1.0 / distancesq;
  r6inv = r2inv*r2inv*r2inv;

  r = sqrt(distancesq);
  rexp = exp(-r*rhoinv[itype][jtype]);
  forcebuck = buck1[itype][jtype] * r*rexp - buck2[itype][jtype] * r6inv;
  fpair = forcebuck*r2inv;
  return fpair;
}
