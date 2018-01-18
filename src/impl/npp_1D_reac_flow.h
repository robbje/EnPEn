#ifndef NPP_1D_H
#define NPP_1D_H
#include <petscts.h>
#include "../app.h"
#include "../iface_petsc.h"

PetscErrorCode npp_1D_reac_flow_init(app_t *ctx);
PetscErrorCode npp_1D_reac_flow_spec(volume_t *v, PetscScalar *in, PetscScalar *dst,
        unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_reac_flow_pois(volume_t *v, PetscScalar *in, PetscScalar *dst,
        app_t *ctx);
PetscErrorCode npp_1D_reac_flow_jac_spec_local(volume_t *v, PetscScalar *in, 
        PetscScalar *dst, unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_reac_flow_jac_pois_local(volume_t *v, PetscScalar *in,
        PetscScalar *dst, app_t *ctx);
PetscErrorCode npp_1D_reac_flow_jac_spec_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int i, unsigned int f, app_t *ctx);
PetscErrorCode npp_1D_reac_flow_jac_pois_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int f, app_t *ctx);
PetscErrorCode npp_1D_reac_flow_create_matrix(app_t *ctx);
PetscErrorCode npp_1D_reac_flow_flux(face_t *f, app_t *ctx);
PetscErrorCode npp_1D_reac_flow_destroy(app_t *ctx);

#endif
