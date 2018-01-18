#ifndef EQUI_H
#define EQUI_H
#include <petscsnes.h>
#include "parameters.h"

#define EQUI_MAX_ITER 100
#define EQUI_TOLF 1E-10
#define EQUI_TOLX 1E-6
typedef struct equi_component_s {
    unsigned int num_acid_reacs;
    unsigned int num_base_reacs;
    double *Ka;
    double *Kb;
    double c;
} equi_component_t;

typedef struct equi_s {
    unsigned int num_components;
    equi_component_t **component;
    unsigned int *tags;
} equi_t;

/* Selfmade solver */
int equi_create(equi_t **equi);
int equi_parse(equi_t *equi, parameter_t *p);
void equi_destroy(equi_t *equi);

int equi_is_acid_reaction(parameter_t *p, unsigned int j);
int equi_add_reaction(equi_t *e, unsigned int reaction, unsigned int component,
        parameter_t *p);
double equi_charge_residual(equi_t *e, double pH);
double equi_solve_ph(equi_t *e);
double equi_bifurcate(double (*func)(), double a, double b, equi_t *e, int n);
void equi_get_initial(equi_t *e, double *c);

/* Petsc solving */
PetscErrorCode equi_solve(PetscScalar *initial, parameter_t *p);
PetscErrorCode equi_function(SNES snes, Vec x, Vec r, void *context);
PetscErrorCode equi_jacobian(SNES snes, Vec x, Mat jac, Mat B,
        void *context);

#endif
