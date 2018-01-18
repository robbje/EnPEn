#include "npp_1D_act.h"

model_t NPP_1D_ACT = {.init = &npp_1D_act_init,
                  .finish = &npp_1D_act_destroy,
                  .fnc_spec = &npp_1D_act_spec,
                  .fnc_pois = &npp_1D_act_pois,
                  .fnc_time = NULL,
                  .jac_spec_local = &npp_1D_act_jac_spec_local,
                  .jac_pois_local = &npp_1D_act_jac_pois_local,
                  .jac_time = NULL,
                  .jac_spec_neighbour = &npp_1D_act_jac_spec_neighbour,
                  .jac_pois_neighbour = &npp_1D_act_jac_pois_neighbour,
                  .fnc_flux = &npp_1D_act_flux,
                  .num_add_vars = 1}; /* electric potential */

PetscScalar rho = 14.9;
PetscScalar A = 0.39;
PetscScalar S = 18.5/997.0;
PetscScalar npp_1D_act_PDH(PetscScalar Ic, PetscInt z)
{
    PetscScalar Ix = S*Ic;
    PetscScalar gamma = \
    exp( \
        -sqrt(1/S) * A * ( \
               (2*z*z/rho) * log(1+rho*sqrt(Ix)) +
               (z*z*sqrt(Ix) - 2*pow(Ix, 3.0/2.0)) /
               (1+rho*sqrt(Ix))
               ) \
       );
    return gamma;
}

PetscScalar npp_1D_act_PDH_diff(PetscScalar Ic, PetscInt z)
{
    PetscScalar Ix = S*Ic;
    /* The exponential term of the differentiation of PDH */
    PetscScalar exp = \
    (-sqrt(1/S)*A)*(
        ((S*pow(z,4))/(4*sqrt(Ix)) - 1.5*S*pow(Ix,2.0)/pow(Ix,3.0/2.0))/
        (rho*sqrt(Ix)+1) +
        (S*pow(z,4)/(2*(rho*sqrt(Ix)+1)*sqrt(Ix))) +
        (S*rho*pow(z,2)*(2*pow(Ix,3.0/2.0)-pow(z,2)*sqrt(Ix)))/
        (4*pow(rho*sqrt(Ix)+1, 2.0)*sqrt(Ix))
    );
    return npp_1D_act_PDH(Ic, z)*exp;
}

PetscErrorCode npp_1D_act_spec(volume_t *v, PetscScalar *in, PetscScalar *dst,
        unsigned int i, app_t *ctx)
{
    model_t *model = (model_t *) ctx->solverdata;
    unsigned int vi = v->index * model->bs;
    parameter_t *p = ctx->p;

    PetscScalar c = in[vi+i];
    PetscScalar phi = in[vi+p->num_species];

    /* Flux term */
    PetscScalar flux = 0;
    unsigned int f;
    for(f = 0; f < v->num_faces; f++) {
        PetscScalar A = v->face[f]->area;
        PetscScalar ksi = v->distance[f];
        volume_t *nv = v->nvol[f];

        PetscScalar c_n, phi_n;
        if(nv->type == VOLTYPE_BULK) {
            c_n = nv->state[i];
            phi_n = nv->state[p->num_species];
        } else {
            unsigned int ni = nv->index * model->bs;
            c_n = in[ni+i];
            phi_n = in[ni+p->num_species];
        }

        PetscScalar alpha = v->alpha;
        if(nv->type == VOLTYPE_MEMBRANE)
            alpha = nv->alpha;

        PetscScalar c_f = 0.5 * (c + c_n);
        flux += A/ksi*alpha*p->D[i] * (c_n-c);
        flux += A/ksi*alpha*p->D[i] * (p->z[i]*c_f*(phi_n-phi));
    }
    *dst = 0 - flux;
    return 0;
}

PetscErrorCode npp_1D_act_pois(volume_t *v, PetscScalar *in, PetscScalar *dst,
        app_t *ctx)
{
    model_t *model = (model_t *) ctx->solverdata;
    unsigned int vi = v->index * model->bs;
    parameter_t *p = ctx->p;

    PetscScalar phi = in[vi+p->num_species];

    /* Source term */
    PetscScalar rho = v->charge, source;
    unsigned int i;
    for(i = 0; i < p->num_species; i++) rho += in[vi+i] * p->z[i];
    source = v->volume * rho * p->debye;

    /* Flux term */
    PetscScalar flux = 0;
    unsigned int f;
    for(f = 0; f < v->num_faces; f++) {
        PetscScalar A = v->face[f]->area;
        PetscScalar ksi = v->distance[f];

        volume_t *nv = v->nvol[f];
        PetscScalar phi_n;
        if(nv->type == VOLTYPE_BULK) {
            phi_n = nv->state[p->num_species];
        } else {
            unsigned int ni = nv->index * model->bs;
            phi_n = in[ni+p->num_species];
        }
        flux += A / ksi * (phi_n - phi);
    }
    *dst = flux + source;
    return 0;
}

PetscErrorCode npp_1D_act_jac_spec_local(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int i, app_t *ctx)
{
    model_t *model = (model_t *) ctx->solverdata;
    unsigned int vi = v->index * model->bs;
    parameter_t *p = ctx->p;

    PetscScalar c = in[vi+i];
    PetscScalar phi = in[vi+p->num_species];

    PetscScalar dc = 0;
    PetscScalar dphi = 0;

    unsigned int f;
    for(f = 0; f < v->num_faces; f++) {
        PetscScalar A = v->face[f]->area;
        PetscScalar ksi = v->distance[f];

        volume_t *nv = v->nvol[f];
        unsigned int ni = nv->index * model->bs;

        PetscScalar alpha = v->alpha;
        if(nv->type == VOLTYPE_MEMBRANE)
            alpha = nv->alpha;

        PetscScalar c_n, phi_n, c_f;
        if(nv->type == VOLTYPE_BULK) {
            c_n = nv->state[i];
            phi_n = nv->state[p->num_species];
        } else {
            c_n = in[ni+i];
            phi_n = in[ni+p->num_species];
        }
        c_f = 0.5 * (c + c_n);

        PetscScalar factor = A*alpha*p->D[i]/ksi;
        dc += factor - factor * p->z[i] * 0.5 * (phi_n-phi);
        dphi += factor * p->z[i] * c_f;
    }
    dst[i] = dc;
    dst[p->num_species] = dphi;
    return 0;
}

PetscErrorCode npp_1D_act_jac_pois_local(volume_t *v, PetscScalar *in,
        PetscScalar *dst, app_t *ctx)
{
    parameter_t *p = ctx->p;

    PetscScalar dphi = 0;

    unsigned int f;
    for(f = 0; f < v->num_faces; f++) {
        PetscScalar A = v->face[f]->area;
        PetscScalar ksi = v->distance[f];
        dphi += -A/ksi;
    }
    dst[p->num_species] = dphi;

    PetscScalar dc = 0;
    unsigned int i;
    for(i = 0; i < p->num_species; i++) {
        dc = v->volume * p->debye * p->z[i];
        dst[i] = dc;
    }
    return 0;
}

PetscErrorCode npp_1D_act_jac_spec_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int i, unsigned int f, app_t *ctx)
{
    model_t *model = (model_t *) ctx->solverdata;
    unsigned int vi = v->index * model->bs;
    parameter_t *p = ctx->p;

    volume_t *nv = v->nvol[f];
    unsigned int ni = nv->index * model->bs;

    PetscScalar A = v->face[f]->area;
    PetscScalar ksi = v->distance[f];
    PetscScalar alpha = v->alpha;
    if(nv->type == VOLTYPE_MEMBRANE)
        alpha = nv->alpha;
    PetscScalar factor = A*alpha*p->D[i]/ksi;

    PetscScalar c, phi, c_n, phi_n, c_f;
    c = in[vi+i];
    phi = in[vi+p->num_species];
    c_n = in[ni+i];
    phi_n = in[ni+p->num_species];
    c_f  = 0.5 * (c + c_n);

    PetscScalar dc = -factor * (1 + p->z[i] * 0.5 * (phi_n-phi));
    PetscScalar dphi = -factor * p->z[i] * c_f;

    dst[i] = dc;
    dst[p->num_species] = dphi;
    return 0;
}

PetscErrorCode npp_1D_act_jac_pois_neighbour(volume_t *v, PetscScalar *in,
        PetscScalar *dst, unsigned int f, app_t *ctx)
{
    parameter_t *p = ctx->p;
    dst[p->num_species] = v->face[f]->area/v->distance[f];
    return 0;
}

PetscErrorCode npp_1D_act_flux(face_t *face, app_t *ctx)
{
    volume_t *v1, *v2;
    v1 = face->nvol[0];
    v2 = face->nvol[1];
    if(!v1||!v2) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
            "A face is not properly connected");

    PetscScalar alpha = v1->alpha;
    if(v2->alpha < v1->alpha)
        alpha = v2->alpha;

    PetscScalar ksi = node_distance(&(v1->center), &(v2->center));
    if(v1->center.x < v2->center.x) ksi *= -1;
    PetscScalar dphi, dc, c_f, cd = 0;
    dphi = (v2->state[ctx->p->num_species] - v1->state[ctx->p->num_species]);

    unsigned int i;
    for(i = 0; i < ctx->p->num_species; i++) {
        dc = v2->state[i] - v1->state[i];
        c_f = 0.5*(v2->state[i] + v1->state[i]);
        face->flux[i] = alpha*ctx->p->D[i]/ksi * (dc + ctx->p->z[i]*c_f*dphi);
        cd += ctx->p->z[i]*face->flux[i];
    }
    face->flux[ctx->p->num_species] = cd;
    return 0;
}

PetscErrorCode npp_1D_act_init(app_t *ctx)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;
    model->bs = ctx->p->num_species + model->num_add_vars;
    unsigned int N = ctx->mesh->num_volumes * model->bs;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &model->x); CHKERRQ(ierr);
    ierr = VecSetBlockSize(model->x, model->bs); CHKERRQ(ierr);
    ierr = VecSetSizes(model->x, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetType(model->x, VECMPI); CHKERRQ(ierr);

    ierr = VecDuplicate(model->x, &model->x1); CHKERRQ(ierr);
    ierr = VecDuplicate(model->x, &model->x2); CHKERRQ(ierr);
    ierr = VecDuplicate(model->x, &model->r); CHKERRQ(ierr);

    /* Create Jacobian */
    ierr = npp_1D_act_create_matrix(ctx); CHKERRQ(ierr);

    ierr = SNESCreate(PETSC_COMM_WORLD, &model->snes);
    ierr = SNESSetFromOptions(model->snes);
    /* Set function and jacobian */ 
    ierr = SNESSetFunction(model->snes, model->r, cb_func_snes, ctx);
    CHKERRQ(ierr);
    ierr = SNESSetJacobian(model->snes, model->J, model->J, cb_jac_snes, ctx);
    CHKERRQ(ierr);
    return 0;
}

PetscErrorCode npp_1D_act_create_matrix(app_t *ctx)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;

    unsigned int n_charged, i;
    for(i = 0; i < ctx->p->num_species; i++)
        n_charged += abs(ctx->p->z[i]);
    PetscInt N = ctx->mesh->num_volumes * model->bs;

    PetscInt rstart, rend;
    ierr = VecGetOwnershipRange(model->x, &rstart, &rend); CHKERRQ(ierr);
    PetscInt nr = rend-rstart;

    PetscInt *d_nnz, *o_nnz;
    ierr = PetscMalloc(nr * sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
    ierr = PetscMalloc(nr * sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);
    memset(d_nnz, 0, nr * sizeof(PetscInt));
    memset(o_nnz, 0, nr * sizeof(PetscInt));

    unsigned int vi, j;
    for(vi = rstart/model->bs; vi < rend/model->bs; vi++) {
        i = vi * model->bs - rstart;
        for(j = 0; j < ctx->p->num_species; j++) {
            d_nnz[i+j] = 3;
        }
        d_nnz[i+ctx->p->num_species] = 3;
    }

    ierr = MatCreateBAIJ(PETSC_COMM_WORLD, model->bs,
            PETSC_DECIDE, PETSC_DECIDE, N, N,
            0, d_nnz, 0, o_nnz, &(model->J)); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "%i Variables\n", N);
    return 0;
}

PetscErrorCode npp_1D_act_destroy(app_t *ctx)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;

    ierr = SNESDestroy(&model->snes); CHKERRQ(ierr);
    ierr = TSDestroy(&model->ts); CHKERRQ(ierr);
    ierr = VecDestroy(&model->x); CHKERRQ(ierr);
    ierr = VecDestroy(&model->x1); CHKERRQ(ierr);
    ierr = VecDestroy(&model->x2); CHKERRQ(ierr);
    ierr = VecDestroy(&model->r); CHKERRQ(ierr);
    ierr = MatDestroy(&model->J); CHKERRQ(ierr);
    return 0;
}
