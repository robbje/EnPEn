#ifndef IFACE_PETSC_H
#define IFACE_PETSC_H
#include "app.h"
#include <petscsys.h>
#include <petscts.h>
#include <stdarg.h>

typedef struct model_s {
    unsigned int num_add_vars;
    PetscErrorCode (*init)(app_t *ctx);
    PetscErrorCode (*finish)(app_t *ctx);
    Mat J;
    Vec x,x1,x2,r;
    SNES snes;
    TS ts;
    PetscInt bs;
    /* User needs to implement (some of) these functions: */
    PetscErrorCode (*fnc_spec)(volume_t *, PetscScalar *, PetscScalar *,
            unsigned int, app_t *);
    PetscErrorCode (*fnc_pois)(volume_t *, PetscScalar *, PetscScalar *,
            app_t *);
    PetscErrorCode (*fnc_time)(volume_t *, PetscScalar *, PetscScalar *,
            PetscScalar *, unsigned int, app_t *);
    PetscErrorCode (*jac_spec_local)(volume_t *, PetscScalar *, PetscScalar *,
            unsigned int, app_t *);
    PetscErrorCode (*jac_pois_local)(volume_t *, PetscScalar *, PetscScalar *,
            app_t *);
    PetscErrorCode (*jac_time)(volume_t *, PetscScalar *, PetscScalar *,
            PetscScalar *, unsigned int i, app_t *);
    PetscErrorCode (*jac_spec_neighbour)(volume_t *, PetscScalar *, PetscScalar *,
            unsigned int, unsigned int, app_t *);
    PetscErrorCode (*jac_pois_neighbour)(volume_t *, PetscScalar *, PetscScalar *,
            unsigned int, app_t *);
    PetscErrorCode (*fnc_flux)(face_t *, app_t *);
} model_t;

int petsc_init(int *, char ***);
int petsc_load(app_t *);
int petsc_finish(app_t *);
int petsc_fopen(const char *, const char *, FILE **);
int petsc_fprintf(FILE *, const char *, ...);
int petsc_getline(FILE *, char **, size_t *);
int petsc_fclose(FILE *);

/* SNES Interface */
PetscErrorCode cb_func_snes(SNES snes, Vec x, Vec r, void *context);
PetscErrorCode cb_jac_snes(SNES snes, Vec x, Mat jac, Mat B,
        void *context);
PetscErrorCode cb_vol_jac(volume_t *v, PetscScalar *in, Mat jac, app_t *ctx,
        PetscInt bs);

/* TS Interface */
PetscErrorCode cb_func_ts(TS ts, PetscReal t, Vec x, Vec xdot, Vec f,
        void *context);
PetscErrorCode cb_jac_ts(TS ts, PetscReal t, Vec x, Vec xdot, PetscReal a,
        Mat jac, Mat B, void *context);
PetscErrorCode cb_vol_jac_ts(volume_t *v, PetscScalar *in, PetscScalar *din,
        Mat jac, app_t *ctx, PetscInt bs, PetscReal a);

/* Other stuff */
PetscErrorCode petsc_check_solution(model_t *model, int *is_positive);
PetscErrorCode copy_solution_from_mesh(app_t *ctx);
PetscErrorCode copy_solution_to_mesh(app_t *ctx);
#endif
