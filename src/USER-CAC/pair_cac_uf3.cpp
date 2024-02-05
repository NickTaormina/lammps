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

#include "pair_cac_uf3.h"
#include "asa_user.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#define MAXNEIGHOUT 50
#define MAXNEIGHIN 10
#define EXPAND 10
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCACUF3::PairCACUF3(LAMMPS *lmp) : PairCAC(lmp)
{
  // CAC
  nmax = 0;
  outer_neighflag = 1;
  flux_enable = 1;
  this->lmp = lmp;

  // UF3
  single_enable = 1;    // 1 if single() routine exists
  restartinfo = 0;      // 1 if pair style writes restart info
  maxshort = 10;
  numshort = 0;
  neighshort = nullptr;
  centroidstressflag = CENTROID_AVAIL;
  manybody_flag = 1;
  one_coeff = 0;    // if 1 then allow only one coeff call of form 'pair_coeff * *'
                    // by setting it to 0 we will allow multiple 'pair_coeff' calls
  bsplines_created = 0;
  num_of_elements = 0;
  nbody_flag = 0;
  n2body_pot_files = 0;
  n3body_pot_files = 0;
  tot_pot_files = 0;
  pot_3b = 0;
  coeff_matrix_dim1 = 0;
  coeff_matrix_dim2 = 0;
  coeff_matrix_dim3 = 0;
  coeff_matrix_elements_len = 0;
  setflag_3b = nullptr;
  cut = nullptr;
  cut_3b = nullptr;
  cut_3b_list = nullptr;
  min_cut_3b = nullptr;
  elements = nullptr;
  map = NULL;
  cluster_neighbors = nullptr;
  cluster_neighbor_counts = nullptr;
  flux_max = 0;
  add_ncluster = 0;
  add_cluster_neighbors = nullptr;
  add_cluster_neighbor_counts = nullptr;
  // CAC
}

/* ---------------------------------------------------------------------- */

PairCACUF3::~PairCACUF3()
{
  if (copymode) return;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    if (pot_3b) {
      memory->destroy(setflag_3b);
      memory->destroy(cut_3b);
      memory->destroy(cut_3b_list);
      memory->destroy(min_cut_3b);
      memory->destroy(neighshort);
    }

    // CAC
    memory->destroy(inner_neighbor_coords);
    memory->destroy(inner_neighbor_types);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCACUF3::allocate()
{
  allocated = 1;
  max_nodes_per_element = atom->nodes_per_element;
  int n = atom->ntypes;
  // Contains info about wether UF potential were found for type i and j
  memory->create(setflag, n + 1, n + 1, "pairCAC:setflag");
  for (int i = 1; i <= num_of_elements; i++)
    for (int j = i; j <= num_of_elements; j++) setflag[i][j] = 0;
  map = new int[n + 1];
  // Contains info about 2-body cutoff distance for type i and j
  // cutsq is the global variable
  // Even though we are making cutsq don't manually change the default values
  // Lammps take care of setting the value
  memory->create(cutsq, n + 1, n + 1, "pairCAC:cutsq");
  // cut is specific to this pair style. We will set the values in cut
  memory->create(cut, num_of_elements + 1, num_of_elements + 1, "pairCAC:cut");
  // memory->create(cut, n + 1, n + 1, "pair:cut_lj"); --- old  CAC allocate
  // Contains knot_vect of 2-body potential for type i and j
  n2b_knot.resize(num_of_elements + 1);
  n2b_coeff.resize(num_of_elements + 1);
  UFBS2b.resize(num_of_elements + 1);
  for (int i = 1; i < num_of_elements + 1; i++) {
    n2b_knot[i].resize(num_of_elements + 1);
    n2b_coeff[i].resize(num_of_elements + 1);
    UFBS2b[i].resize(num_of_elements + 1);
  }

  if (pot_3b) {
    // Contains info about wether UF potential were found for type i, j and k
    memory->create(setflag_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1,
                   "pairCAC:setflag_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    memory->create(cut_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1,
                   "pairCAC:cut_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    // for constructing 3-body list
    memory->create(cut_3b_list, num_of_elements + 1, num_of_elements + 1, "pairCAC:cut_3b_list");
    // Contains info about minimum 3-body cutoff distance for type i, j and k
    memory->create(min_cut_3b, num_of_elements + 1, num_of_elements + 1, num_of_elements + 1, 3,
                   "pairCAC:min_cut_3b");

    // setting cut_3b and setflag = 0
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        cut_3b_list[i][j] = 0;
        for (int k = 1; k < num_of_elements + 1; k++) {
          cut_3b[i][j][k] = 0;
          min_cut_3b[i][j][k][0] = 0;
          min_cut_3b[i][j][k][1] = 0;
          min_cut_3b[i][j][k][2] = 0;
        }
      }
    }
    n3b_knot_matrix.resize(num_of_elements + 1);
    UFBS3b.resize(num_of_elements + 1);
    for (int i = 1; i < num_of_elements + 1; i++) {
      n3b_knot_matrix[i].resize(num_of_elements + 1);
      UFBS3b[i].resize(num_of_elements + 1);
      for (int j = 1; j < num_of_elements + 1; j++) {
        n3b_knot_matrix[i][j].resize(num_of_elements + 1);
        UFBS3b[i][j].resize(num_of_elements + 1);
      }
    }
    memory->create(neighshort, maxshort, "pair:neighshort");
  }

  // normal allocate function from lammps, then allocate the CAC arrays
  PairCAC::allocate();
}

void PairCACUF3::uf3_read_pot_file(int itype, int jtype, char *potf_name)
{
  utils::logmesg(lmp, "UF3: {} file should contain UF3 potential for {} {}\n", potf_name, itype,
                 jtype);

  FILE *fp;
  fp = utils::open_potential(potf_name, lmp, nullptr);

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  std::string temp_line = txtfilereader.next_line(1);
  Tokenizer file_header(temp_line);

  if (file_header.count() != 2)
    error->all(FLERR, "UF3: Expected only two words on 1st line of {} but found \n\
            {} word/s",
               potf_name, file_header.count());

  if (file_header.contains("#UF3 POT") == 0)
    error->all(FLERR, "UF3: {} file is not UF3 POT type, 1st line of UF3 POT \n\
            files contain '#UF3 POT'. Found {} in the header",
               potf_name, temp_line);

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp2nd_line(temp_line);

  if (fp2nd_line.count() != 3)
    error->all(FLERR, "UF3: Expected 3 words on 2nd line =>\n\
            nBody leading_trim trailing_trim\n\
            Found {}",
               temp_line);

  std::string nbody_on_file = fp2nd_line.next_string();
  if (utils::strmatch(nbody_on_file, "2B"))
    utils::logmesg(lmp, "UF3: File {} contains 2-body UF3 potential\n", potf_name);
  else
    error->all(FLERR, "UF3: Expected a 2B UF3 file but found {}", nbody_on_file);

  int leading_trim = fp2nd_line.next_int();
  int trailing_trim = fp2nd_line.next_int();
  if (leading_trim != 0)
    error->all(FLERR, "UF3: Current implementation is throughly tested only for\n\
            leading_trim=0\n");
  if (trailing_trim != 3)
    error->all(FLERR, "UF3: Current implementation is throughly tested only for\n\
            trailing_trim=3\n");

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp3rd_line(temp_line);
  if (fp3rd_line.count() != 2)
    error->all(FLERR, "UF3: Expected only 2 numbers on 3rd line =>\n\
            Rij_CUTOFF NUM_OF_KNOTS\n\
            Found {} number/s",
               fp3rd_line.count());

  // cut is used in init_one which is called by pair.cpp at line 267 where the return of init_one is squared abc
  cut[itype][jtype] = fp3rd_line.next_double();
  if (comm->me == 0)
    utils::logmesg(lmp, "UF3 pairpot: Cutoff {} {} {}\n", cut[itype][jtype], itype, jtype);
  cut[jtype][itype] = cut[itype][jtype];

  int num_knots_2b = fp3rd_line.next_int();

  temp_line = txtfilereader.next_line(num_knots_2b);
  ValueTokenizer fp4th_line(temp_line);

  if (fp4th_line.count() != num_knots_2b)
    error->all(FLERR, "UF3: Expected {} numbers on 4th line but found {} numbers", num_knots_2b,
               fp4th_line.count());

  n2b_knot[itype][jtype].resize(num_knots_2b);
  n2b_knot[jtype][itype].resize(num_knots_2b);
  for (int k = 0; k < num_knots_2b; k++) {
    n2b_knot[itype][jtype][k] = fp4th_line.next_double();
    n2b_knot[jtype][itype][k] = n2b_knot[itype][jtype][k];
  }

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp5th_line(temp_line);
  int num_of_coeff_2b = fp5th_line.next_int();

  temp_line = txtfilereader.next_line(num_of_coeff_2b);
  ValueTokenizer fp6th_line(temp_line);

  if (fp6th_line.count() != num_of_coeff_2b)
    error->all(FLERR, "UF3: Expected {} numbers on 6th line but found {} numbers", num_of_coeff_2b,
               fp6th_line.count());

  n2b_coeff[itype][jtype].resize(num_of_coeff_2b);
  n2b_coeff[jtype][itype].resize(num_of_coeff_2b);
  for (int k = 0; k < num_of_coeff_2b; k++) {
    n2b_coeff[itype][jtype][k] = fp6th_line.next_double();
    n2b_coeff[jtype][itype][k] = n2b_coeff[itype][jtype][k];
  }

  if (n2b_knot[itype][jtype].size() != n2b_coeff[itype][jtype].size() + 4) {
    error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
               potf_name);
  }
  setflag[itype][jtype] = 1;
  setflag[jtype][itype] = 1;
}

void PairCACUF3::uf3_read_pot_file(int itype, int jtype, int ktype, char *potf_name)
{
  utils::logmesg(lmp, "UF3: {} file should contain UF3 potential for {} {} {}\n", potf_name, itype,
                 jtype, ktype);

  FILE *fp;
  fp = utils::open_potential(potf_name, lmp, nullptr);

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  std::string temp_line = txtfilereader.next_line(1);
  Tokenizer file_header(temp_line);

  if (file_header.count() != 2)
    error->all(FLERR, "UF3: Expected only two words on 1st line of {} but found \n\
            {} word/s",
               potf_name, file_header.count());

  if (file_header.contains("#UF3 POT") == 0)
    error->all(FLERR, "UF3: {} file is not UF3 POT type, 1st line of UF3 POT \n\
            files contain '#UF3 POT'. Found {} in the header",
               potf_name, temp_line);

  temp_line = txtfilereader.next_line(1);
  ValueTokenizer fp2nd_line(temp_line);

  if (fp2nd_line.count() != 3)
    error->all(FLERR, "UF3: Expected 3 words on 2nd line =>\n\
            nBody leading_trim trailing_trim\n\
            Found {}",
               temp_line);

  std::string nbody_on_file = fp2nd_line.next_string();

  if (utils::strmatch(nbody_on_file, "3B"))
    utils::logmesg(lmp, "UF3: File {} contains 3-body UF3 potential\n", potf_name);
  else
    error->all(FLERR, "UF3: Expected a 3B UF3 file but found {}", nbody_on_file);

  int leading_trim = fp2nd_line.next_int();
  int trailing_trim = fp2nd_line.next_int();
  if (leading_trim != 0)
    error->all(FLERR, "UF3: Current implementation is throughly tested only for\n\
            leading_trim=0\n");
  if (trailing_trim != 3)
    error->all(FLERR, "UF3: Current implementation is throughly tested only for\n\
            trailing_trim=3\n");

  temp_line = txtfilereader.next_line(6);
  ValueTokenizer fp3rd_line(temp_line);

  if (fp3rd_line.count() != 6)
    error->all(FLERR, "UF3: Expected only 6 numbers on 3rd line =>\n\
            Rjk_CUTOFF Rik_CUTOFF Rij_CUTOFF NUM_OF_KNOTS_JK NUM_OF_KNOTS_IK NUM_OF_KNOTS_IJ\n\
            Found {} number/s",
               fp3rd_line.count());

  double cut3b_rjk = fp3rd_line.next_double();
  double cut3b_rij = fp3rd_line.next_double();
  double cut3b_rik = fp3rd_line.next_double();

  if (cut3b_rij != cut3b_rik) {
    error->all(FLERR, "UF3: rij!=rik, Current implementation only works for rij=rik");
  }

  if (2 * cut3b_rik != cut3b_rjk) {
    error->all(FLERR, "UF3: 2rij=2rik!=rik, Current implementation only works \n\
            for 2rij=2rik!=rik");
  }

  cut_3b_list[itype][jtype] = std::max(cut3b_rij, cut_3b_list[itype][jtype]);
  cut_3b_list[itype][ktype] = std::max(cut_3b_list[itype][ktype], cut3b_rik);

  cut_3b[itype][jtype][ktype] = cut3b_rij;
  cut_3b[itype][ktype][jtype] = cut3b_rik;

  int num_knots_3b_jk = fp3rd_line.next_int();
  temp_line = txtfilereader.next_line(num_knots_3b_jk);
  ValueTokenizer fp4th_line(temp_line);

  if (fp4th_line.count() != num_knots_3b_jk)
    error->all(FLERR, "UF3: Expected {} numbers on 4th line but found {} numbers", num_knots_3b_jk,
               fp4th_line.count());

  n3b_knot_matrix[itype][jtype][ktype].resize(3);
  n3b_knot_matrix[itype][ktype][jtype].resize(3);

  n3b_knot_matrix[itype][jtype][ktype][0].resize(num_knots_3b_jk);
  n3b_knot_matrix[itype][ktype][jtype][0].resize(num_knots_3b_jk);

  for (int i = 0; i < num_knots_3b_jk; i++) {
    n3b_knot_matrix[itype][jtype][ktype][0][i] = fp4th_line.next_double();
    n3b_knot_matrix[itype][ktype][jtype][0][i] = n3b_knot_matrix[itype][jtype][ktype][0][i];
  }

  min_cut_3b[itype][jtype][ktype][0] = n3b_knot_matrix[itype][jtype][ktype][0][0];
  min_cut_3b[itype][ktype][jtype][0] = n3b_knot_matrix[itype][ktype][jtype][0][0];
  if (comm->me == 0)
    utils::logmesg(lmp, "UF3: 3b min cutoff {} {}-{}-{}_0={} {}-{}-{}_0={}\n", potf_name, itype,
                   jtype, ktype, min_cut_3b[itype][jtype][ktype][0], itype, ktype, jtype,
                   min_cut_3b[itype][ktype][jtype][0]);

  int num_knots_3b_ik = fp3rd_line.next_int();
  temp_line = txtfilereader.next_line(num_knots_3b_ik);
  ValueTokenizer fp5th_line(temp_line);

  if (fp5th_line.count() != num_knots_3b_ik)
    error->all(FLERR, "UF3: Expected {} numbers on 5th line but found {} numbers", num_knots_3b_ik,
               fp5th_line.count());

  n3b_knot_matrix[itype][jtype][ktype][1].resize(num_knots_3b_ik);
  n3b_knot_matrix[itype][ktype][jtype][2].resize(num_knots_3b_ik);
  for (int i = 0; i < num_knots_3b_ik; i++) {
    n3b_knot_matrix[itype][jtype][ktype][1][i] = fp5th_line.next_double();
    n3b_knot_matrix[itype][ktype][jtype][2][i] = n3b_knot_matrix[itype][jtype][ktype][1][i];
  }

  min_cut_3b[itype][jtype][ktype][1] = n3b_knot_matrix[itype][jtype][ktype][1][0];
  min_cut_3b[itype][ktype][jtype][2] = n3b_knot_matrix[itype][ktype][jtype][2][0];
  if (comm->me == 0)
    utils::logmesg(lmp, "UF3: 3b min cutoff {} {}-{}-{}_1={} {}-{}-{}_2={}\n", potf_name, itype,
                   jtype, ktype, min_cut_3b[itype][jtype][ktype][1], itype, ktype, jtype,
                   min_cut_3b[itype][ktype][jtype][2]);

  int num_knots_3b_ij = fp3rd_line.next_int();
  temp_line = txtfilereader.next_line(num_knots_3b_ij);
  ValueTokenizer fp6th_line(temp_line);

  if (fp6th_line.count() != num_knots_3b_ij)
    error->all(FLERR, "UF3: Expected {} numbers on 6th line but found {} numbers", num_knots_3b_ij,
               fp5th_line.count());

  n3b_knot_matrix[itype][jtype][ktype][2].resize(num_knots_3b_ij);
  n3b_knot_matrix[itype][ktype][jtype][1].resize(num_knots_3b_ij);
  for (int i = 0; i < num_knots_3b_ij; i++) {
    n3b_knot_matrix[itype][jtype][ktype][2][i] = fp6th_line.next_double();
    n3b_knot_matrix[itype][ktype][jtype][1][i] = n3b_knot_matrix[itype][jtype][ktype][2][i];
  }

  min_cut_3b[itype][jtype][ktype][2] = n3b_knot_matrix[itype][jtype][ktype][2][0];
  min_cut_3b[itype][ktype][jtype][1] = n3b_knot_matrix[itype][ktype][jtype][1][0];
  if (comm->me == 0)
    utils::logmesg(lmp, "UF3: 3b min cutoff {} {}-{}-{}_2={} {}-{}-{}_1={}\n", potf_name, itype,
                   jtype, ktype, min_cut_3b[itype][jtype][ktype][2], itype, ktype, jtype,
                   min_cut_3b[itype][ktype][jtype][1]);

  temp_line = txtfilereader.next_line(3);
  ValueTokenizer fp7th_line(temp_line);

  if (fp7th_line.count() != 3)
    error->all(FLERR, "UF3: Expected 3 numbers on 7th line =>\n\
           SHAPE_OF_COEFF_MATRIX[I][J][K] \n\
           found {} numbers",
               fp7th_line.count());

  coeff_matrix_dim1 = fp7th_line.next_int();
  coeff_matrix_dim2 = fp7th_line.next_int();
  coeff_matrix_dim3 = fp7th_line.next_int();

  if (n3b_knot_matrix[itype][jtype][ktype][0].size() != coeff_matrix_dim3 + 3 + 1)
    error->all(FLERR, "UF3: {} has incorrect knot (NUM_OF_KNOTS_JK) and \n\
            coeff (coeff_matrix_dim3) data nknots!=ncoeffs + 3 +1",
               potf_name);

  if (n3b_knot_matrix[itype][jtype][ktype][1].size() != coeff_matrix_dim2 + 3 + 1)
    error->all(FLERR, "UF3: {} has incorrect knot (NUM_OF_KNOTS_IK) and \n\
            coeff (coeff_matrix_dim2) data nknots!=ncoeffs + 3 +1",
               potf_name);

  if (n3b_knot_matrix[itype][jtype][ktype][2].size() != coeff_matrix_dim1 + 3 + 1)
    error->all(FLERR, "UF3: {} has incorrect knot (NUM_OF_KNOTS_IJ) and \n\
            coeff ()coeff_matrix_dim1 data nknots!=ncoeffs + 3 +1",
               potf_name);

  coeff_matrix_elements_len = coeff_matrix_dim3;

  std::string key = std::to_string(itype) + std::to_string(jtype) + std::to_string(ktype);
  n3b_coeff_matrix[key].resize(coeff_matrix_dim1);

  int line_count = 0;
  for (int i = 0; i < coeff_matrix_dim1; i++) {
    n3b_coeff_matrix[key][i].resize(coeff_matrix_dim2);
    for (int j = 0; j < coeff_matrix_dim2; j++) {
      temp_line = txtfilereader.next_line(coeff_matrix_elements_len);
      ValueTokenizer coeff_line(temp_line);
      n3b_coeff_matrix[key][i][j].resize(coeff_matrix_dim3);

      if (coeff_line.count() != coeff_matrix_elements_len)
        error->all(FLERR, "UF3: Expected {} numbers on {}th line but found \n\
                {} numbers",
                   coeff_matrix_elements_len, line_count + 8, coeff_line.count());
      for (int k = 0; k < coeff_matrix_dim3; k++) {
        n3b_coeff_matrix[key][i][j][k] = coeff_line.next_double();
      }
      line_count += 1;
    }
  }

  std::string key2 = std::to_string(itype) + std::to_string(ktype) + std::to_string(jtype);
  n3b_coeff_matrix[key2].resize(coeff_matrix_dim2);
  for (int j = 0; j < coeff_matrix_dim2; j++) {
    n3b_coeff_matrix[key2][j].resize(coeff_matrix_dim1);
    for (int i = 0; i < coeff_matrix_dim1; i++) {
      n3b_coeff_matrix[key2][j][i].resize(coeff_matrix_dim3);
    }
  }

  for (int i = 0; i < coeff_matrix_dim1; i++) {
    for (int j = 0; j < coeff_matrix_dim2; j++) {
      for (int k = 0; k < coeff_matrix_dim3; k++) {
        n3b_coeff_matrix[key2][j][i][k] = n3b_coeff_matrix[key][i][j][k];
      }
    }
  }

  setflag_3b[itype][jtype][ktype] = 1;
  setflag_3b[itype][ktype][jtype] = 1;
}

void PairCACUF3::uf3_read_pot_file(char *potf_name)
{
  if (comm->me == 0) utils::logmesg(lmp, "\nUF3: Opening {} file\n", potf_name);

  FILE *fp;
  fp = utils::open_potential(potf_name, lmp, nullptr);
  // if (fp) error->all(FLERR,"UF3: {} file not found",potf_name);

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  std::string temp_line = txtfilereader.next_line(2);
  Tokenizer fp1st_line(temp_line);

  if (fp1st_line.contains("#UF3 POT") == 0)
    error->all(FLERR, "UF3: {} file is not UF3 POT type, found type {} {} on the file", potf_name,
               fp1st_line.next(), fp1st_line.next());

  if (comm->me == 0)
    utils::logmesg(lmp, "UF3: {} file is of type {} {}\n", potf_name, fp1st_line.next(),
                   fp1st_line.next());

  temp_line = txtfilereader.next_line(1);
  Tokenizer fp2nd_line(temp_line);
  if (fp2nd_line.contains("2B") == 1) {
    temp_line = txtfilereader.next_line(4);
    ValueTokenizer fp3rd_line(temp_line);
    int temp_type1 = fp3rd_line.next_int();
    int temp_type2 = fp3rd_line.next_int();
    if (comm->me == 0)
      utils::logmesg(lmp, "UF3: {} file contains 2-body UF3 potential for {} {}\n", potf_name,
                     temp_type1, temp_type2);

    // cut is used in init_one which is called by pair.cpp at line 267 where the return of init_one is squared
    // test
    cut[temp_type1][temp_type2] = fp3rd_line.next_double();
    if (comm->me == 0) utils::logmesg(lmp, "test");
    if (comm->me == 0) utils::logmesg(lmp, "UF3 pairpot: Cutoff {}\n", cut[temp_type1][temp_type2]);
    cut[temp_type2][temp_type1] = cut[temp_type1][temp_type2];

    int temp_line_len = fp3rd_line.next_int();

    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp4th_line(temp_line);

    n2b_knot[temp_type1][temp_type2].resize(temp_line_len);
    n2b_knot[temp_type2][temp_type1].resize(temp_line_len);
    for (int k = 0; k < temp_line_len; k++) {
      n2b_knot[temp_type1][temp_type2][k] = fp4th_line.next_double();
      n2b_knot[temp_type2][temp_type1][k] = n2b_knot[temp_type1][temp_type2][k];
    }

    temp_line = txtfilereader.next_line(1);
    ValueTokenizer fp5th_line(temp_line);

    temp_line_len = fp5th_line.next_int();

    temp_line = txtfilereader.next_line(temp_line_len);
    // utils::logmesg(lmp,"UF3:11 {}",temp_line); //test
    ValueTokenizer fp6th_line(temp_line);
    // if(comm->me==0) utils::logmesg(lmp,"UF3: {}\n",temp_line_len);
    n2b_coeff[temp_type1][temp_type2].resize(temp_line_len);
    n2b_coeff[temp_type2][temp_type1].resize(temp_line_len);

    for (int k = 0; k < temp_line_len; k++) {
      n2b_coeff[temp_type1][temp_type2][k] = fp6th_line.next_double();
      n2b_coeff[temp_type2][temp_type1][k] = n2b_coeff[temp_type1][temp_type2][k];
      // if(comm->me==0) utils::logmesg(lmp,"UF3: {}\n",n2b_coeff[temp_type1][temp_type2][k]);
    }
    // for(int i=0;i<n2b_coeff[temp_type1][temp_type2].size();i++) if(comm->me==0) utils::logmesg(lmp,"UF3: {}\n",n2b_coeff[temp_type1][temp_type2][i]);
    if (n2b_knot[temp_type1][temp_type2].size() != n2b_coeff[temp_type1][temp_type2].size() + 4) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }
    setflag[temp_type1][temp_type2] = 1;
    setflag[temp_type2][temp_type1] = 1;
  } else if (fp2nd_line.contains("3B") == 1) {
    temp_line = txtfilereader.next_line(9);
    ValueTokenizer fp3rd_line(temp_line);
    int temp_type1 = fp3rd_line.next_int();
    int temp_type2 = fp3rd_line.next_int();
    int temp_type3 = fp3rd_line.next_int();
    if (comm->me == 0)
      utils::logmesg(lmp, "UF3: {} file contains 3-body UF3 potential for {} {} {}\n", potf_name,
                     temp_type1, temp_type2, temp_type3);

    double cut3b_rjk = fp3rd_line.next_double();
    double cut3b_rij = fp3rd_line.next_double();
    // cut_3b[temp_type1][temp_type2] = std::max(cut3b_rij,
    // cut_3b[temp_type1][temp_type2]);
    cut_3b_list[temp_type1][temp_type2] = std::max(cut3b_rij, cut_3b_list[temp_type1][temp_type2]);
    double cut3b_rik = fp3rd_line.next_double();
    if (cut3b_rij != cut3b_rik) {
      error->all(FLERR, "UF3: rij!=rik, Current implementation only works for rij=rik");
    }
    if (2 * cut3b_rik != cut3b_rjk) {
      error->all(FLERR,
                 "UF3: 2rij=2rik!=rik, Current implementation only works for 2rij=2rik!=rik");
    }
    // cut_3b[temp_type1][temp_type3] = std::max(cut_3b[temp_type1][temp_type3],cut3b_rik);
    cut_3b_list[temp_type1][temp_type3] = std::max(cut_3b_list[temp_type1][temp_type3], cut3b_rik);

    cut_3b[temp_type1][temp_type2][temp_type3] = cut3b_rij;
    cut_3b[temp_type1][temp_type3][temp_type2] = cut3b_rik;

    int temp_line_len = fp3rd_line.next_int();
    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp4th_line(temp_line);

    n3b_knot_matrix[temp_type1][temp_type2][temp_type3].resize(3);
    n3b_knot_matrix[temp_type1][temp_type3][temp_type2].resize(3);

    n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0].resize(temp_line_len);
    n3b_knot_matrix[temp_type1][temp_type3][temp_type2][0].resize(temp_line_len);
    for (int i = 0; i < temp_line_len; i++) {
      n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0][i] = fp4th_line.next_double();
      n3b_knot_matrix[temp_type1][temp_type3][temp_type2][0][i] =
          n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0][i];
    }

    min_cut_3b[temp_type1][temp_type2][temp_type3][0] =
        n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0][0];
    min_cut_3b[temp_type1][temp_type3][temp_type2][0] =
        n3b_knot_matrix[temp_type1][temp_type3][temp_type2][0][0];
    if (comm->me == 0)
      utils::logmesg(lmp, "UF3: 3b min cutoff {} {}-{}-{}_0={} {}-{}-{}_0={}\n", potf_name,
                     temp_type1, temp_type2, temp_type3,
                     min_cut_3b[temp_type1][temp_type2][temp_type3][0], temp_type1, temp_type3,
                     temp_type2, min_cut_3b[temp_type1][temp_type3][temp_type2][0]);

    temp_line_len = fp3rd_line.next_int();
    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp5th_line(temp_line);
    n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1].resize(temp_line_len);
    n3b_knot_matrix[temp_type1][temp_type3][temp_type2][2].resize(temp_line_len);
    for (int i = 0; i < temp_line_len; i++) {
      n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1][i] = fp5th_line.next_double();
      n3b_knot_matrix[temp_type1][temp_type3][temp_type2][2][i] =
          n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1][i];
    }

    min_cut_3b[temp_type1][temp_type2][temp_type3][1] =
        n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1][0];
    min_cut_3b[temp_type1][temp_type3][temp_type2][2] =
        n3b_knot_matrix[temp_type1][temp_type3][temp_type2][2][0];
    if (comm->me == 0)
      utils::logmesg(lmp, "UF3: 3b min cutoff {} {}-{}-{}_1={} {}-{}-{}_2={}\n", potf_name,
                     temp_type1, temp_type2, temp_type3,
                     min_cut_3b[temp_type1][temp_type2][temp_type3][1], temp_type1, temp_type3,
                     temp_type2, min_cut_3b[temp_type1][temp_type3][temp_type2][2]);

    temp_line_len = fp3rd_line.next_int();
    temp_line = txtfilereader.next_line(temp_line_len);
    ValueTokenizer fp6th_line(temp_line);
    n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2].resize(temp_line_len);
    n3b_knot_matrix[temp_type1][temp_type3][temp_type2][1].resize(temp_line_len);
    for (int i = 0; i < temp_line_len; i++) {
      n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2][i] = fp6th_line.next_double();
      n3b_knot_matrix[temp_type1][temp_type3][temp_type2][1][i] =
          n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2][i];
    }

    min_cut_3b[temp_type1][temp_type2][temp_type3][2] =
        n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2][0];
    min_cut_3b[temp_type1][temp_type3][temp_type2][1] =
        n3b_knot_matrix[temp_type1][temp_type3][temp_type2][1][0];
    if (comm->me == 0)
      utils::logmesg(lmp, "UF3: 3b min cutoff {} {}-{}-{}_2={} {}-{}-{}_1={}\n", potf_name,
                     temp_type1, temp_type2, temp_type3,
                     min_cut_3b[temp_type1][temp_type2][temp_type3][2], temp_type1, temp_type3,
                     temp_type2, min_cut_3b[temp_type1][temp_type3][temp_type2][2]);

    temp_line = txtfilereader.next_line(3);
    ValueTokenizer fp7th_line(temp_line);

    coeff_matrix_dim1 = fp7th_line.next_int();
    coeff_matrix_dim2 = fp7th_line.next_int();
    coeff_matrix_dim3 = fp7th_line.next_int();
    if (n3b_knot_matrix[temp_type1][temp_type2][temp_type3][0].size() !=
        coeff_matrix_dim3 + 3 + 1) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }
    if (n3b_knot_matrix[temp_type1][temp_type2][temp_type3][1].size() !=
        coeff_matrix_dim2 + 3 + 1) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }
    if (n3b_knot_matrix[temp_type1][temp_type2][temp_type3][2].size() !=
        coeff_matrix_dim1 + 3 + 1) {
      error->all(FLERR, "UF3: {} has incorrect knot and coeff data nknots!=ncoeffs + 3 +1",
                 potf_name);
    }

    coeff_matrix_elements_len = coeff_matrix_dim3;

    std::string key =
        std::to_string(temp_type1) + std::to_string(temp_type2) + std::to_string(temp_type3);
    n3b_coeff_matrix[key].resize(coeff_matrix_dim1);
    for (int i = 0; i < coeff_matrix_dim1; i++) {
      n3b_coeff_matrix[key][i].resize(coeff_matrix_dim2);
      for (int j = 0; j < coeff_matrix_dim2; j++) {
        temp_line = txtfilereader.next_line(coeff_matrix_elements_len);
        ValueTokenizer coeff_line(temp_line);
        n3b_coeff_matrix[key][i][j].resize(coeff_matrix_dim3);
        for (int k = 0; k < coeff_matrix_dim3; k++) {
          n3b_coeff_matrix[key][i][j][k] = coeff_line.next_double();
        }
      }
    }

    key = std::to_string(temp_type1) + std::to_string(temp_type3) + std::to_string(temp_type2);
    n3b_coeff_matrix[key] =
        n3b_coeff_matrix[std::to_string(temp_type1) + std::to_string(temp_type2) +
                         std::to_string(temp_type3)];
    setflag_3b[temp_type1][temp_type2][temp_type3] = 1;
    setflag_3b[temp_type1][temp_type3][temp_type2] = 1;
  } else
    error->all(
        FLERR,
        "UF3: {} file does not contain right words indicating whether it is 2 or 3 body potential",
        potf_name);
}

/* ----------------------------------------------------------------------
global settings. pair_style call in script
narg is literally number of args that are passed in. we want 1 or 2
char**arg is the arguments themselves
uf3
------------------------------------------------------------------------- */
void PairCACUF3::settings(int narg, char **arg)
{

  if (narg != 2) { error->all(FLERR, "Illegal pair_style cac/uf3 command"); }

  nbody_flag =
      utils::numeric(FLERR, arg[0], true, lmp);    // how many body interactions to consider
  num_of_elements = utils::numeric(FLERR, arg[1], true,
                                   lmp);    // number of elements in the system. atoms ->ntypes;

  if (num_of_elements !=
      atom->ntypes) {    // why is this even here if you can just take from lammps
    error->all(FLERR, "Number of elements in the system does not match the number of atom types");
  }

  if (nbody_flag == 2) {
    pot_3b = false;
    n2body_pot_files = num_of_elements * (num_of_elements + 1) / 2;
    tot_pot_files = n2body_pot_files;
  } else if (nbody_flag == 3) {
    pot_3b = true;
    n2body_pot_files = num_of_elements * (num_of_elements + 1) / 2;
    n3body_pot_files = num_of_elements * (num_of_elements * (num_of_elements + 1) / 2);
    tot_pot_files = n2body_pot_files + n3body_pot_files;
  } else
    error->all(FLERR, "UF3: UF3 not yet implemented for {}-body", nbody_flag);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs. pair_coeff.
  narg is literally number of args that are passed in. we want 3 or 4
  uf3
 ------------------------------------------------------------------------- */

void PairCACUF3::coeff(int narg, char **arg)
{
  if (narg != 3 && narg != 4)
    error->all(FLERR, "UF3: WARNING!! It seems that you are using the \n\
             older style of specifying UF3 POT files. This style of listing \n\
             all the potential files on a single line will be depcrecated in \n\
             the next version of ML-UF3");
  if (!allocated) allocate();

  if (narg == 3 || narg == 4) {
    int ilo, ihi, jlo, jhi, klo, khi;
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
    if (narg == 4) utils::bounds(FLERR, arg[2], 1, atom->ntypes, klo, khi, error);

    if (narg == 3) {
      if (narg == 3) {
        for (int i = ilo; i <= ihi; i++) {
          for (int j = MAX(jlo, i); j <= jhi; j++) {
            if (comm->me == 0) utils::logmesg(lmp, "\nUF3: Opening {} file\n", arg[2]);
            uf3_read_pot_file(i, j, arg[2]);
          }
        }
      }
    }

    if (narg == 4) {
      for (int i = ilo; i <= ihi; i++) {
        for (int j = jlo; j <= jhi; j++) {
          for (int k = MAX(klo, jlo); k <= khi; k++) {
            if (comm->me == 0) utils::logmesg(lmp, "\nUF3: Opening {} file\n", arg[3]);
            uf3_read_pot_file(i, j, k, arg[3]);
          }
        }
      }
    }
  } else {
    if (narg != tot_pot_files + 2)
      error->all(FLERR, "UF3: Invalid number of argument in pair coeff; \n\
              Number of potential files provided is not correct");

    error->warning(FLERR, "\nUF3: WARNING!! It seems that you are using the \n\
            older style of specifying UF3 POT files. This style of listing \n\
            all the potential files on a single line will be depcrecated in \n\
            the next version of ML-UF3");

    // open UF3 potential file on all proc
    for (int i = 2; i < narg; i++) { uf3_read_pot_file(arg[i]); }
    if (!bsplines_created) create_bsplines();
  }
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 Basic lammps
 uf3
 ------------------------------------------------------------------------- */

double PairCACUF3::init_one(int i, int j)
{
  if (!bsplines_created) create_bsplines();
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

// in normal uf3, netwon pair is forced on here
void PairCACUF3::init_style()
{
  PairCAC::init_style();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  atom->outer_neigh_flag = 1;
}

void PairCACUF3::create_bsplines()
{
  bsplines_created = 1;
  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = 1; j < num_of_elements + 1; j++) {
      if (setflag[i][j] != 1)
        error->all(FLERR, "UF3: Not all 2-body UF potentials are set, \n\
                missing potential file for {}-{} interaction",
                   i, j);
    }
  }
  if (pot_3b) {
    for (int i = 1; i < num_of_elements + 1; i++) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = 1; k < num_of_elements + 1; k++) {
          if (setflag_3b[i][j][k] != 1)
            error->all(FLERR, "UF3: Not all 3-body UF potentials are set, \n\
                    missing potential file for {}-{}-{} interaction",
                       i, j, k);
        }
      }
    }
  }

  for (int i = 1; i < num_of_elements + 1; i++) {
    for (int j = i; j < num_of_elements + 1; j++) {
      UFBS2b[i][j] = uf3_pair_bspline(lmp, n2b_knot[i][j], n2b_coeff[i][j]);
      UFBS2b[j][i] = UFBS2b[i][j];
    }
    if (pot_3b) {
      for (int j = 1; j < num_of_elements + 1; j++) {
        for (int k = j; k < num_of_elements + 1; k++) {
          std::string key = std::to_string(i) + std::to_string(j) + std::to_string(k);
          UFBS3b[i][j][k] =
              uf3_triplet_bspline(lmp, n3b_knot_matrix[i][j][k], n3b_coeff_matrix[key]);
          std::string key2 = std::to_string(i) + std::to_string(k) + std::to_string(j);
          UFBS3b[i][k][j] =
              uf3_triplet_bspline(lmp, n3b_knot_matrix[i][k][j], n3b_coeff_matrix[key2]);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
/*
  CAC-specific code below
*/
void PairCACUF3::force_densities(int iii, double s, double t, double w, double coefficients,
                                 double &force_densityx, double &force_densityy,
                                 double &force_densityz)
{

  double delx, dely, delz;
  double cutshortsq = cutmax * cutmax;
  int timestep = update->ntimestep;
  double *special_lj = force->special_lj;
  double fpair, flux_interaction[3];
  int *type = atom->type;
  double distancesq;
  double scan_position[3], scan_position2[3];
  int current_type = poly_counter;
  double rcut;
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;
  int origin_type = type_array[poly_counter];
  int listtype;
  int scan_type, scan_type2, scan_type3;
  int element_index;
  int poly_index, short_scan, short_scan2, short_scan3, short_scan_index;
  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int neigh_max_outer = outer_quad_lists_counts[pqi];
  int neigh_max_add, all_neigh;
  int itype, jtype, ktype;
  double energy_contribution;
  double rsq, rsq1, rsq2, rsq3;
  double rij, rik, rkj, del_rji[3], del_rki[3], del_rkj[3];
  double delr1[3], delr2[3], ndelr1[3], delr3[3], fij[3], fik[3], fjk[3];
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  double inner_scan_position[3];
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  int **outer_quad_indices = outer_quad_lists_index[pqi];
  double cut_add = atom->cut_add;
  energy_contribution = 0;

  //printf("pqi is: %d\n", pqi);

  // allocate arrays that store neighbor information around just this quadrature point
  if (neigh_max_inner > local_inner_max) {
    memory->grow(cluster_neighbor_counts, neigh_max_inner + EXPAND + 1,
                 "pair_cac_uf3:cluster_neighbor_counts");
    cluster_neighbors =
        (int **) memory->srealloc(cluster_neighbors, sizeof(int *) * (neigh_max_inner + EXPAND + 1),
                                  "pair_cac_uf3:cluster_neighbors");
    // initialize sizing of cluster neighbors using neigh_max_inner
    for (int cinit = 0; cinit < neigh_max_inner + EXPAND + 1; cinit++) {
      if (cinit >= local_inner_max + 1 || local_inner_max == 0) cluster_neighbors[cinit] = NULL;
      memory->grow(cluster_neighbors[cinit], neigh_max_inner + neigh_max_outer + EXPAND,
                   "pair_cac_uf3:inner_neighbor_types");
    }
  }
  allocate_quad_memory();
  // set virtual neighbor types, etc.
  init_quad_arrays();
  // interpolate virtual atom coordinates from shape functions corresponding to unit cells
  interpolation(iii, s, t, w);

  // initialize cluster counts
  for (int cinit = 0; cinit < local_inner_max + 1; cinit++) { cluster_neighbor_counts[cinit] = 0; }

  int counter = 0;
  // two body contribution
  for (int l = 0; l < neigh_max_inner; l++) {

    scan_type = inner_neighbor_types[l];
    scan_position[0] = inner_neighbor_coords[l][0];
    scan_position[1] = inner_neighbor_coords[l][1];
    scan_position[2] = inner_neighbor_coords[l][2];
    delx = current_position[0] - scan_position[0];
    dely = current_position[1] - scan_position[1];
    delz = current_position[2] - scan_position[2];
    distancesq = delx * delx + dely * dely + delz * delz;

    if (distancesq <= cutsq[origin_type][scan_type]) {
      rij = sqrt(distancesq);
      twobody(rij, fpair, 1, scan_type, quad_eflag, energy_contribution);
      // quadrature_energy += energy_contribution / 2;
      force_densityx += delx * fpair;
      force_densityy += dely * fpair;
      force_densityz += delz * fpair;
      // if 3 body potential is used, store this pair in the shortlist
      if (pot_3b) {
        if (rij <= cut_3b_list[origin_type][scan_type]) {
          if (counter >= neigh_max_inner) {
            error->all(FLERR, "3-body counter exceeded");
          } else {
            cluster_neighbors[0][cluster_neighbor_counts[0]++] = l;
            counter++;
          }
        }
      }
      if (atom->CAC_virial) {}
      // cac flux contribution due to current quadrature point and neighbor pair interactions
    }
  }

  // make list of cluster neighbors in the short cutoff
  // construct the list of cluster neighbors using the short cutoff
  for (int l = 0; l < cluster_neighbor_counts[0]; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];
    for (int j = 0; j < neigh_max_inner; j++) {
      if (short_scan == j) continue;
      scan_position2[0] = inner_neighbor_coords[j][0];
      scan_position2[1] = inner_neighbor_coords[j][1];
      scan_position2[2] = inner_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
      if (distancesq < cut_3b_list[origin_type][scan_type])
        cluster_neighbors[l + 1][cluster_neighbor_counts[l + 1]++] = j;
    }
    for (int j = 0; j < neigh_max_outer; j++) {
      scan_position2[0] = outer_neighbor_coords[j][0];
      scan_position2[1] = outer_neighbor_coords[j][1];
      scan_position2[2] = outer_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
      if (distancesq < cut_3b_list[origin_type][scan_type])
        cluster_neighbors[l + 1][cluster_neighbor_counts[l + 1]++] = j + neigh_max_inner;
    }
  }

  int loop_limit = cluster_neighbor_counts[0] - 1;

  // ith three body contributions
  // loop over neighbors (j)
  for (int l = 0; l < loop_limit; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_type = inner_neighbor_types[short_scan];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];
    delr1[0] = delx = scan_position[0] - current_position[0];
    delr1[1] = dely = scan_position[1] - current_position[1];
    delr1[2] = delz = scan_position[2] - current_position[2];

    rsq1 = sqrt(delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2]);    // rij
    if (rsq1 >= cut_3b_list[origin_type][scan_type]) continue;
    if (quad_flux_flag) flux_interaction[0] = flux_interaction[1] = flux_interaction[2] = 0;
    for (int k = l + 1; k < cluster_neighbor_counts[0]; k++) {
      short_scan2 = cluster_neighbors[0][k];
      if (short_scan == short_scan2) continue;
      scan_type2 = inner_neighbor_types[short_scan2];
      scan_position2[0] = inner_neighbor_coords[short_scan2][0];
      scan_position2[1] = inner_neighbor_coords[short_scan2][1];
      scan_position2[2] = inner_neighbor_coords[short_scan2][2];

      delr2[0] = scan_position2[0] - current_position[0];
      delr2[1] = scan_position2[1] - current_position[1];
      delr2[2] = scan_position2[2] - current_position[2];
      rsq2 = sqrt(delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2]);    // rik
      if (rsq2 >= cut_3b[origin_type][scan_type2][scan_type2]) continue;

      // now check if everything is in the cutoff
      if ((rsq1 <= cut_3b[origin_type][scan_type][scan_type]) &&
          (rsq2 <= cut_3b[origin_type][scan_type][scan_type2]) &&
          (rsq1 >= min_cut_3b[origin_type][scan_type][scan_type][0]) &&
          (rsq2 >= min_cut_3b[origin_type][scan_type][scan_type2][1])) {
        del_rkj[0] = scan_position2[0] - scan_position[0];
        del_rkj[1] = scan_position2[1] - scan_position[1];
        del_rkj[2] = scan_position2[2] - scan_position[2];
        rkj = sqrt(
            ((del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2])));
        if (rkj >= min_cut_3b[origin_type][scan_type][scan_type2][2]) {
          double distances[] = {rsq1, rsq2, rkj};
          double *del[] = {delr1, delr2, del_rkj};
          int type[] = {origin_type, scan_type, scan_type2};
          threebody(distances, del, fij, fik, fjk, type, quad_eflag, energy_contribution);
          // quadrature_energy += energy_contribution / 3;

          // substracting it first because itll get added back later?
          force_densityx -= fij[0] + fik[0];
          force_densityy -= fij[1] + fik[1];
          force_densityz -= fij[2] + fik[2];

          if (atom->CAC_virial) {}
        }
      }
    }
  }
  // jk three body contribution to i
  for (int l = 0; l < loop_limit; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_type = inner_neighbor_types[short_scan];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];
    delr1[0] = current_position[0] - scan_position[0];
    delr1[1] = current_position[1] - scan_position[1];
    delr1[2] = current_position[2] - scan_position[2];

    rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];    // rij
    if (rsq1 >= cut_3b_list[scan_type][origin_type]) continue;
    if (quad_flux_flag) flux_interaction[0] = flux_interaction[1] = flux_interaction[2] = 0;

    for (int k = 0; k < cluster_neighbor_counts[l + 1]; k++) {
      short_scan2 = cluster_neighbors[l + 1][k];
      if (short_scan2 >= neigh_max_inner) {
        short_scan2 -= neigh_max_inner;
        scan_type2 = outer_neighbor_types[short_scan2];
        scan_position2[0] = outer_neighbor_coords[short_scan2][0];
        scan_position2[1] = outer_neighbor_coords[short_scan2][1];
        scan_position2[2] = outer_neighbor_coords[short_scan2][2];
      } else {
        scan_type2 = inner_neighbor_types[short_scan2];
        scan_position2[0] = inner_neighbor_coords[short_scan2][0];
        scan_position2[1] = inner_neighbor_coords[short_scan2][1];
        scan_position2[2] = inner_neighbor_coords[short_scan2][2];
      }

      // ik param = scan-scan2
      // ijk param = scan-origin-scan2

      delr2[0] = scan_position2[0] - scan_position[0];
      delr2[1] = scan_position2[1] - scan_position[1];
      delr2[2] = scan_position2[2] - scan_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];    // rik
      if (rsq2 >= cut_3b[scan_type][scan_type2][scan_type2]) continue;
      del_rkj[0] = current_position[0] - scan_position2[0];
      del_rkj[1] = current_position[1] - scan_position2[1];
      del_rkj[2] = current_position[2] - scan_position2[2];
      rkj =
          sqrt(((del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2])));

      double distances[] = {rsq1, rsq2, rkj};
      double *del[] = {delr1, delr2, del_rkj};
      int type[] = {scan_type, origin_type, scan_type2};
      threebody(distances, del, fij, fik, fjk, type, quad_eflag, energy_contribution);
      // quadrature_energy += energy_contribution / 3;

      force_densityx += fij[0];
      force_densityy += fij[1];
      force_densityz += fij[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCACUF3::twobody(double distancesq, double &fforce, int itype, int jtype, int eflag,
                         double &eng)
{
  // two body
  double *pair_eval = UFBS2b[itype][jtype].eval(distancesq);
  fforce = -1 * pair_eval[1] / (distancesq);

  if (eflag) {
    // eng = pair_eval[1]; // PLACEHOLDER
  }
}

// computes the UF3 three body potential.
void PairCACUF3::threebody(double *r, double **del, double *fij, double *fik, double *fjk,
                           int *types, int eflag, double &eng)
{
  // three body placeholder
  double *triangle_eval = UFBS3b[types[0]][types[1]][types[2]].eval(r[0], r[1], r[2]);

  // x-axis
  fij[0] = *(triangle_eval + 1) * (del[0][0] / r[0]);
  fik[0] = *(triangle_eval + 2) * (del[1][0] / r[1]);
  fjk[0] = *(triangle_eval + 3) * (del[2][0] / r[2]);

  // y-axis
  fij[1] = *(triangle_eval + 1) * (del[0][1] / r[0]);
  fik[1] = *(triangle_eval + 2) * (del[1][1] / r[1]);
  fjk[1] = *(triangle_eval + 3) * (del[2][1] / r[2]);

  // z-axis
  fij[2] = *(triangle_eval + 1) * (del[0][2] / r[0]);
  fik[2] = *(triangle_eval + 2) * (del[1][2] / r[1]);
  fjk[2] = *(triangle_eval + 3) * (del[2][2] / r[2]);
}
