#ifndef NPP_1D_H
#define NPP_1D_H
#include <petscts.h>
#include "../app.h"
#include "../iface_petsc.h"

PetscErrorCode npp_1D_reac_ts_init(app_t *ctx);
PetscErrorCode npp_1D_reac_ts_spec(volume_t *v, PetscScalar *in, PetscScalar *dst,
        unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_pois(volume_t *v, PetscScalar *in, PetscScalar *dst,
        app_t *ctx);
PetscErrorCode npp_1D_reac_ts_time(volume_t *v, PetscScalar *in, PetscScalar *din,
        PetscScalar *dst, unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_jac_spec_local(volume_t *v, PetscScalar *in, 
        PetscScalar *dst, unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_jac_pois_local(volume_t *v, PetscScalar *in,
        PetscScalar *dst, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_jac_time(volume_t *v, PetscScalar *in,
        PetscScalar *din, PetscScalar *dst, unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_jac_spec_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int i, unsigned int f, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_jac_pois_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int f, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_create_matrix(app_t *ctx);
PetscErrorCode npp_1D_reac_ts_flux(face_t *f, app_t *ctx);
PetscErrorCode npp_1D_reac_ts_destroy(app_t *ctx);

#endif
