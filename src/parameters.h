#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "interface.h"
#include "misc.h"
#include <stdio.h>  // fprintf
#include <stdlib.h> // malloc
#include <errno.h>  // errno
#include <string.h> // strtok

typedef struct io_s io_t;
typedef struct parameter_s {
    unsigned int num_species;           /* Number of species */
    unsigned int num_reactions;         /* Number of reactions */
    unsigned int num_groups;            /* Number of functional groups */
    int *z;                             /* valence of species */
    int *group_valence;                 /* valence of functional group */
    double *mw;                         /* molecular weight */
    double *D;                          /* diffusion coeff. */
    double *hydration;                  /* hydration number */
    int **mu;                           /* forward reaction coeff. */
    int **nu;                           /* backward reaction coeff. */
    int **group_nu;                     /* group coeff. */
    double *k_f, *k_b;                  /* reaction speeds */
    double T, eta, rho;                 /* temperature, viscosity, density */
    double closest_approach;            /* PDH closest approach parameter */
    double F, R;                        /* Faraday const., id. gas const. */
    double epsilon;                     /* medium permittivity */
    char **name;                        /* species names */
    double c0;                          /* concentration normalization */
    double D0;                          /* diffusion coeff. normalization */
    double mu0;                         /* viscosity normalization */
    double L0;                          /* length normalization */
    double p0;                          /* pressure normalization */
    double u0;                          /* velocity normalization */
    double t0;                          /* time normalization */
    double phi0;                        /* potential normalization */
    double debye;                       /* system debye number */
    double kappa;                       /* membrane permeability */
    char *filename;                     /* system definition filename */
    double u;                           /* convective velocity, dimensionless */
    struct param_impedance_s {
        double f;
        double bias;
        double modulus;
    } impedance;
} parameter_t;

typedef struct param_norm_s {
    double c0;
    double D0;
} param_norm_t;


void param_set(parameter_t *p, char *filename, param_norm_t n);
int param_load(parameter_t *p, io_t *io);
void param_destroy(parameter_t *p);
double param_get_debye_length(parameter_t *p, double *c);
#endif
