/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(cac/pb,PairCACPb)

#else

#ifndef LMP_PAIR_PB_CAC_H
#define LMP_PAIR_PB_CAC_H

//#include "asa_user.h"
#include "pair_cac.h"


namespace LAMMPS_NS {

class PairCACPb : public PairCAC {
 public:
	 PairCACPb(class LAMMPS *);
  virtual ~PairCACPb();
 
  
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);

 protected:
             
  int neigh_nodes_per_element;
	
  double **cut, v_sh, rexp, r6inv;
  double **a, **rho, **c;
  double **rhoinv, **buck1, **buck2, **offset;
  double cut_coul, cut_coulsq, alf, cut_buck, e_shift, f_shift;
  double  factor_coul;
  //debug array
  double force_contrib_buck[500], force_contrib_coul[500], fbuckcontrib, fcoulcontrib;

  void allocate();
  void force_densities(int, double, double, double, double, double
	  &fx, double &fy, double &fz);
  
  virtual void settings(int, char **);
  virtual double pair_interaction_q(double, int, int, double, double);
};

}

#endif
#endif
