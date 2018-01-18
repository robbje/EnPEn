#ifndef NPP_1D_ACT_H
#define NPP_1D_ACT_H
#include <petscts.h>
#include "../app.h"
#include "../iface_petsc.h"

PetscScalar npp_1D_act_PDH(PetscScalar Ic, PetscInt z);
PetscScalar npp_1D_act_PDH_diff(PetscScalar Ic, PetscInt z);
PetscErrorCode npp_1D_act_init(app_t *ctx);
PetscErrorCode npp_1D_act_spec(volume_t *v, PetscScalar *in, PetscScalar *dst,
        unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_act_pois(volume_t *v, PetscScalar *in, PetscScalar *dst,
        app_t *ctx);
PetscErrorCode npp_1D_act_jac_spec_local(volume_t *v, PetscScalar *in, 
        PetscScalar *dst, unsigned int i, app_t *ctx);
PetscErrorCode npp_1D_act_jac_pois_local(volume_t *v, PetscScalar *in,
        PetscScalar *dst, app_t *ctx);
PetscErrorCode npp_1D_act_jac_spec_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int i, unsigned int f, app_t *ctx);
PetscErrorCode npp_1D_act_jac_pois_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int f, app_t *ctx);
PetscErrorCode npp_1D_act_create_matrix(app_t *ctx);
PetscErrorCode npp_1D_act_flux(face_t *f, app_t *ctx);
PetscErrorCode npp_1D_act_destroy(app_t *ctx);

#endif
