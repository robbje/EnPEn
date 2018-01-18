#ifndef SIMULATION_H
#define SIMULATION_H
#include "app.h"
#include "iface_petsc.h"
#include "interface.h"
#include "modules.h"
#include "equi.h"
#include <petscsys.h>
#include <petscts.h>
#include <stdio.h>

typedef struct sim_iter_s {
    struct sim_iter_s *next;
    double value;
} sim_iter_t;

PetscErrorCode snes_solve(app_t *app);
PetscErrorCode sweep_zone_adaptive(app_t *app, double *p, double start, double end);
PetscErrorCode sweep_zones(app_t *app);
PetscErrorCode sweep_boundary_chargeall(app_t *app, int bc, double *initial,
        double *p, double start, double end, double step);
PetscErrorCode sweep_boundary_with_charge(app_t *app, int bc, double *initial,
        double *p, double start, double end, double step, zone_t *z);
PetscErrorCode sweep_flow(app_t *app, double *p, double start, double end, double step);
PetscErrorCode sweep_boundary_chargeall_flow(app_t *app, int bc, double *initial,
    double *p, double start, double end, double step, double *u, double u_start,
    double u_end, double u_step);
PetscErrorCode time_current_control(app_t *app, double *initial, unsigned int n);
PetscErrorCode time_ivc(app_t *app, double *initial, unsigned int n);
PetscErrorCode tsmon(TS ts, PetscInt steps, PetscReal time, Vec u, void *context);
#endif
