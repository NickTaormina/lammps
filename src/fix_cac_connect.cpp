#include "fix_cac_connect.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "respa.h"
#include "update.h"
#include <cstdio>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------------------- */
/*
    fix (name) (group-ID) cac/connect
*/
FixConnectCAC::FixConnectCAC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  for (int i = 0; i < narg; i++) { std::cout << "arg: " << arg[i] << std::endl; }
}

/* ---------------------------------------------------------------------- */
int FixConnectCAC::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixConnectCAC::setup(int vflag)
{
  if (!atom->CAC_flag) error->all(FLERR, "fix cac/connect requires a CAC atom style");

  // check for overlapping nodes
  double ****nodal_positions = atom->nodal_positions;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int *nodes_count_list = atom->nodes_per_element_list;
  int nodes_per_element;
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;

  int totalCount = 0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      nodes_per_element = nodes_count_list[element_type[i]];
      for (int ipoly = 0; ipoly < poly_count[i]; ipoly++) {
        for (int n = 0; n < nodes_per_element; n++) { totalCount++; }
      }
    }
  }

  std::cout << "Total number of nodes: " << totalCount << std::endl;
}