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

PairStyle(cac/uf3, PairCACUF3)

#else

#ifndef LMP_PAIR_UF3_CAC_H
#define LMP_PAIR_UF3_CAC_H

#include "pair_cac.h"
#include "uf3_pair_bspline.h"
#include "uf3_triplet_bspline.h"
#include <unordered_map>
#include <string>
#include "tokenizer.h"
#include "text_file_reader.h"

namespace LAMMPS_NS
{

    class PairCACUF3 : public PairCAC
    {
    public:
        PairCACUF3(class LAMMPS *);
        virtual ~PairCACUF3();

        void coeff(int, char **);
        virtual void init_style();
        virtual double init_one(int, int);

    protected:
        // uf3
        LAMMPS *lmp;
        void uf3_read_pot_file(char *potf_name);
        void uf3_read_pot_file(int i, int j, char *potf_name);
        void uf3_read_pot_file(int i, int j, int k, char *potf_name);
        int num_of_elements, nbody_flag, n2body_pot_files, n3body_pot_files, tot_pot_files;
        int bsplines_created;
        int coeff_matrix_dim1, coeff_matrix_dim2, coeff_matrix_dim3, coeff_matrix_elements_len;
        bool pot_3b;
        int ***setflag_3b;
        double **cut, ***cut_3b, **cut_3b_list, ****min_cut_3b;
        void create_bsplines();
        std::vector<std::vector<std::vector<double>>> n2b_knot, n2b_coeff;
        std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> n3b_knot_matrix;
        std::unordered_map<std::string, std::vector<std::vector<std::vector<double>>>> n3b_coeff_matrix;
        std::vector<std::vector<uf3_pair_bspline>> UFBS2b;
        std::vector<std::vector<std::vector<uf3_triplet_bspline>>> UFBS3b;
        // short neighbor list variables

        // neighbor arrays
        char **elements; // names of unique elements
        int *map;        // mapping from atom types to elements
        int **
            cluster_neighbors; // stores neighbors of neighbors within cutshort for the given quadrature point
        int *
            cluster_neighbor_counts; // number of neighbors withing cutshort for each neighbor within cutshort
        int flux_max;                // array storage maximum for additional cluster neighbor array for flux calculation
        int add_ncluster;            // number of additional sites to store neighbors around for the flux calculation
        int **
            add_cluster_neighbors; // stores neighbors of neighbors for flux calculation around the quadrature point
        int *
            add_cluster_neighbor_counts; // stores neighbors of neighbors counts for flux calculation around the quadrature point

        // previous uf3-cac
        int *neighshort, maxshort;                      // short neighbor list array for 3body interaction
        int numshort, j, jtype, jnum, jj, k, kk, ktype; // short neighbor list variables
        double rij, rik, rjk;

        // cac
        void allocate();
        void force_densities(int, double, double, double, double, double &fx, double &fy, double &fz);
        virtual void settings(int, char **);
        void twobody(double distancesq, double &fforce, int itype, int jtype, int eflag,
                     double &eng);
        void threebody(double *, double **, double *fij, double *fik, double *fjk, int *, int eflag,
                       double &eng);

        // virtual void quad_neigh_flux(); TODO
    };

} // namespace LAMMPS_NS

#endif
#endif

    /* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unexpected argument in pair cac/buck invocation; only accepts cutoff and the 'one' keyword

Self-explanatory.  Check the input script. See the documentation for the proper syntax.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
