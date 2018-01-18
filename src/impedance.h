#ifndef IMPEDANCE_H
#define IMPEDANCE_H

#include "app.h"
#include "mesh.h"
#include "parameters.h"
#include "simulation.h"
#include <petscts.h>
#include <math.h>
#define TAU 2*3.14159265359

void impedance_set_boundaries(app_t *app, PetscReal time);
PetscErrorCode impedance_sampler(TS ts, PetscInt steps, PetscReal time, Vec u,
                                 void *ctx);
PetscErrorCode impedance_setter(TS ts);
PetscErrorCode impedance_setup(app_t *app, double *initial, unsigned int n);
PetscErrorCode impedance_run(app_t *app);

#endif
