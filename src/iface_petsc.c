#include "interface.h"
#include "iface_petsc.h"
PETSC_EXTERN PetscErrorCode (*PetscVFPrintf)(FILE*,const char[],va_list);

io_t IO_PETSC = {
    .fopen = petsc_fopen,
    .fclose = petsc_fclose,
    .fprintf = petsc_fprintf,
    .getline = petsc_getline
};

const interface_t IFACE_PETSC = {
    .init = petsc_init,
    .load = petsc_load,
    .finish = petsc_finish,
};

PetscErrorCode cb_func_ts(TS ts, PetscReal t, Vec x, Vec xdot, Vec f,
        void *context)
{
    PetscErrorCode ierr;
    app_t *ctx = (app_t *) context;
    model_t *model = (model_t *) ctx->solverdata;

    VecScatter sctr;
    Vec xglobal,xdotglobal;
    ierr = VecScatterCreateToAll(x, &sctr, &xglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);

    ierr = VecScatterCreateToAll(xdot, &sctr, &xdotglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, xdot, xdotglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, xdot, xdotglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);

    PetscInt rstart, rend, bs;
    ierr = VecGetOwnershipRange(x, &rstart, &rend); CHKERRQ(ierr);
    ierr = VecGetBlockSize(x, &bs); CHKERRQ(ierr);

    PetscScalar *in, *din, *out, *dst;
    ierr = VecGetArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecGetArray(xdotglobal, &din); CHKERRQ(ierr);
    ierr = VecGetArray(f, &out); CHKERRQ(ierr);

    unsigned int vi;
    for(vi = rstart; vi < rend; vi+=bs) {
        dst = &out[vi - rstart];
        volume_t *v = ctx->mesh->volumes[vi/bs];
        unsigned int i;
        for(i = 0; i < ctx->p->num_species; i++) {
            ierr = model->fnc_spec(v, in, &dst[i], i, ctx); CHKERRQ(ierr);
            PetscScalar t_term;
            ierr = model->fnc_time(v, in, din, &t_term, i, ctx); CHKERRQ(ierr);
            dst[i] += t_term;
        }
        ierr = model->fnc_pois(v, in, &dst[ctx->p->num_species], ctx);
        CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecRestoreArray(xdotglobal, &din); CHKERRQ(ierr);
    ierr = VecRestoreArray(f, &out); CHKERRQ(ierr);

    ierr = VecDestroy(&xglobal); CHKERRQ(ierr);
    ierr = VecDestroy(&xdotglobal); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode cb_func_snes(SNES snes, Vec x, Vec f, void *context)
{
    PetscErrorCode ierr;
    app_t *ctx = (app_t *) context;
    model_t *model = (model_t *) ctx->solverdata;

    VecScatter sctr;
    Vec xglobal;
    ierr = VecScatterCreateToAll(x, &sctr, &xglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);

    PetscInt rstart, rend, bs;
    ierr = VecGetOwnershipRange(x, &rstart, &rend); CHKERRQ(ierr);
    ierr = VecGetBlockSize(x, &bs); CHKERRQ(ierr);
    
    PetscScalar *in, *out, *dst;
    ierr = VecGetArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecGetArray(f, &out); CHKERRQ(ierr);

    unsigned int vi;
    for(vi = rstart; vi < rend; vi += bs) {
        dst = &out[vi - rstart];
        unsigned int i;
        volume_t *v = ctx->mesh->volumes[vi/bs];
        for(i = 0; i < ctx->p->num_species; i++) {
            ierr = model->fnc_spec(v, in, &dst[i], i, ctx); CHKERRQ(ierr);
        }
        ierr = model->fnc_pois(v, in, &dst[ctx->p->num_species], ctx);
        CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecRestoreArray(f, &out); CHKERRQ(ierr);

    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);
    ierr = VecDestroy(&xglobal); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode cb_jac_ts(TS ts, PetscReal t, Vec x, Vec xdot, PetscReal a,
        Mat jac, Mat B, void *context)
{
    PetscErrorCode ierr;
    app_t *ctx = (app_t *) context;

    VecScatter sctr;
    Vec xglobal,xdotglobal;
    ierr = VecScatterCreateToAll(x, &sctr, &xglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);

    ierr = VecScatterCreateToAll(xdot, &sctr, &xdotglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, xdot, xdotglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, xdot, xdotglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);

    PetscInt rstart, rend, bs;
    ierr = VecGetOwnershipRange(x, &rstart, &rend); CHKERRQ(ierr);
    ierr = VecGetBlockSize(x, &bs); CHKERRQ(ierr);

    PetscScalar *in, *din;
    ierr = VecGetArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecGetArray(xdotglobal, &din); CHKERRQ(ierr);

    unsigned int vi;
    ierr = MatZeroEntries(jac); CHKERRQ(ierr);
    for(vi = rstart; vi < rend; vi+=bs) {
        ierr = cb_vol_jac_ts(ctx->mesh->volumes[vi/bs], in, din, jac, ctx, bs, a);
        CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = VecRestoreArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecRestoreArray(xdotglobal, &din); CHKERRQ(ierr);

    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);
    ierr = VecDestroy(&xglobal); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode cb_jac_snes(SNES snes, Vec x, Mat jac, Mat B,
        void *context)
{
    PetscErrorCode ierr;
    app_t *ctx = (app_t *) context;

    VecScatter sctr;
    Vec xglobal;
    ierr = VecScatterCreateToAll(x, &sctr, &xglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, x, xglobal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);

    PetscInt rstart, rend, bs;
    ierr = VecGetOwnershipRange(x, &rstart, &rend); CHKERRQ(ierr);
    ierr = VecGetBlockSize(x, &bs); CHKERRQ(ierr);
    
    PetscScalar *in;
    ierr = VecGetArray(xglobal, &in); CHKERRQ(ierr);

    unsigned int vi;
    ierr = MatZeroEntries(jac); CHKERRQ(ierr);
    for(vi = rstart; vi < rend; vi+=bs) {
        ierr = cb_vol_jac(ctx->mesh->volumes[vi/bs], in, jac, ctx, bs);
        CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = VecRestoreArray(xglobal, &in); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);
    ierr = VecDestroy(&xglobal); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode cb_vol_jac_ts(volume_t *v, PetscScalar *in, PetscScalar *din,
        Mat jac, app_t *ctx, PetscInt bs, PetscReal a)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;
    
    /* This will hold the values for the block for volume v */
    PetscScalar *values;
    ierr = PetscMalloc(bs * bs * sizeof(PetscScalar), &values);
    CHKERRQ(ierr);
    memset(values, 0, bs * bs * sizeof(PetscScalar));

    /* These will hold the indices to set */
    PetscInt *idx, *cidx;
    ierr = PetscMalloc(bs * sizeof(PetscInt), &idx); CHKERRQ(ierr);
    ierr = PetscMalloc(bs * sizeof(PetscInt), &cidx); CHKERRQ(ierr);

    /* Local volume */
    unsigned int i;
    unsigned int vi = v->index * bs;
    for(i = 0; i < ctx->p->num_species; i++) {
        idx[i] = vi+i;
        ierr = model->jac_spec_local(v, in, &values[i*bs], i, ctx);
        CHKERRQ(ierr);
        PetscScalar t_term;
        ierr = model->jac_time(v, in, din, &t_term, i, ctx);
        values[i*bs+i] += a * t_term;
        CHKERRQ(ierr);
    }
    idx[ctx->p->num_species] = vi+ctx->p->num_species;
    ierr = model->jac_pois_local(v, in, &values[bs*ctx->p->num_species], ctx);
    CHKERRQ(ierr);

    ierr = MatSetValues(jac, bs, idx, bs, idx, values, INSERT_VALUES);
    CHKERRQ(ierr);

    /* All neighbouring volumes */
    unsigned int f, ni;
    for(f = 0; f < v->num_faces; f++) {
        volume_t *nv = v->nvol[f];
        /* Ignore boundary conditions */
        if(nv->type == VOLTYPE_BULK ||
           nv->type == VOLTYPE_ELECTRODE) continue;
        memset(values, 0, bs * bs * sizeof(PetscScalar));

        ni = nv->index * bs;
        for(i = 0; i < ctx->p->num_species; i++) {
            idx[i] = vi+i;
            cidx[i] = ni+i;
            ierr = model->jac_spec_neighbour(v, in, &values[i*bs], i, f, ctx);
            CHKERRQ(ierr);
        }

        idx[ctx->p->num_species] = vi + ctx->p->num_species;
        cidx[ctx->p->num_species] = ni + ctx->p->num_species;
        ierr = model->jac_pois_neighbour(v, in, &values[i*bs],
                f, ctx); CHKERRQ(ierr);

        ierr = MatSetValues(jac, bs, idx, bs, cidx, values, INSERT_VALUES);
        CHKERRQ(ierr);
    }

    ierr = PetscFree(values); CHKERRQ(ierr);
    ierr = PetscFree(idx); CHKERRQ(ierr);
    ierr = PetscFree(cidx); CHKERRQ(ierr); // Make these static?
    return 0;
}

PetscErrorCode cb_vol_jac(volume_t *v, PetscScalar *in, Mat jac, app_t *ctx,
        PetscInt bs)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;
    
    /* This will hold the values for the block for volume v */
    PetscScalar *values;
    ierr = PetscMalloc(bs * bs * sizeof(PetscScalar), &values);
    CHKERRQ(ierr);
    memset(values, 0, bs * bs * sizeof(PetscScalar));

    /* These will hold the indices to set */
    PetscInt *idx, *cidx;
    ierr = PetscMalloc(bs * sizeof(PetscInt), &idx); CHKERRQ(ierr);
    ierr = PetscMalloc(bs * sizeof(PetscInt), &cidx); CHKERRQ(ierr);

    /* Local volume */
    unsigned int i;
    unsigned int vi = v->index * bs;
    for(i = 0; i < ctx->p->num_species; i++) {
        idx[i] = vi+i;
        ierr = model->jac_spec_local(v, in, &values[i*bs], i, ctx);
        CHKERRQ(ierr);
    }
    idx[ctx->p->num_species] = vi+ctx->p->num_species;
    ierr = model->jac_pois_local(v, in, &values[bs*ctx->p->num_species], ctx);
    CHKERRQ(ierr);

    ierr = MatSetValues(jac, bs, idx, bs, idx, values, INSERT_VALUES);
    CHKERRQ(ierr);

    /* All neighbouring volumes */
    unsigned int f, ni;
    for(f = 0; f < v->num_faces; f++) {
        volume_t *nv = v->nvol[f];
        /* Ignore boundary conditions */
        if(nv->type == VOLTYPE_BULK ||
           nv->type == VOLTYPE_DOWNSTREAM ||
           nv->type == VOLTYPE_ELECTRODE) continue;
        memset(values, 0, bs * bs * sizeof(PetscScalar));

        ni = nv->index * bs;
        for(i = 0; i < ctx->p->num_species; i++) {
            idx[i] = vi+i;
            cidx[i] = ni+i;
            ierr = model->jac_spec_neighbour(v, in, &values[i*bs], i, f, ctx);
            CHKERRQ(ierr);
        }

        idx[ctx->p->num_species] = vi + ctx->p->num_species;
        cidx[ctx->p->num_species] = ni + ctx->p->num_species;
        ierr = model->jac_pois_neighbour(v, in, &values[i*bs],
                f, ctx); CHKERRQ(ierr);

        ierr = MatSetValues(jac, bs, idx, bs, cidx, values, INSERT_VALUES);
        CHKERRQ(ierr);
    }

    ierr = PetscFree(values); CHKERRQ(ierr);
    ierr = PetscFree(idx); CHKERRQ(ierr);
    ierr = PetscFree(cidx); CHKERRQ(ierr); // Make these static?
    return 0;
}

PetscErrorCode petsc_check_solution(model_t *model, int *is_positive)
{
    PetscErrorCode ierr;
    PetscScalar *x;
    PetscInt rstart, rend;
    ierr = VecGetOwnershipRange(model->x, &rstart, &rend); CHKERRQ(ierr);
    ierr = VecGetArray(model->x, &x); CHKERRQ(ierr);
    PetscInt i;
    *is_positive = 1;
    for(i = rstart; i < rend; i++) {
        if((i+1) % model->bs > 0 && x[i] < 0) {
            *is_positive ^= *is_positive;
            break;
        }
    }
    ierr = VecRestoreArray(model->x, &x); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode copy_solution_from_mesh(app_t *ctx)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;

    PetscInt rstart, rend, bs;
    ierr = VecGetOwnershipRange(model->x, &rstart, &rend); CHKERRQ(ierr);
    ierr = VecGetBlockSize(model->x, &bs); CHKERRQ(ierr);

    unsigned int vi, i;
    for(vi = rstart/bs; vi < rend/bs; vi++) {
        double *state = ctx->mesh->volumes[vi]->state;
        for(i = 0; i < ctx->p->num_species; i++) {
            ierr = VecSetValue(model->x, vi*bs+i, state[i], INSERT_VALUES);
            CHKERRQ(ierr);
        }
        ierr = VecSetValue(model->x, vi*bs+ctx->p->num_species,
                state[ctx->p->num_species], INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(model->x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(model->x); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode copy_solution_to_mesh(app_t *ctx)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;

    PetscInt bs;
    ierr = VecGetBlockSize(model->x, &bs); CHKERRQ(ierr);

    VecScatter sctr;
    Vec xglobal;
    ierr = VecScatterCreateToAll(model->x, &sctr, &xglobal); CHKERRQ(ierr);
    ierr = VecScatterBegin(sctr, model->x, xglobal, INSERT_VALUES,
            SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(sctr, model->x, xglobal, INSERT_VALUES,
            SCATTER_FORWARD); CHKERRQ(ierr);

    PetscScalar *x;
    ierr = VecGetArray(xglobal, &x); CHKERRQ(ierr);

    unsigned int vi;
    for(vi = 0; vi < ctx->mesh->num_volumes; vi++) {
        double *dst = ctx->mesh->volumes[vi]->state;
        memcpy(dst, &x[vi*bs], bs*sizeof(double));
    }

    ierr = VecRestoreArray(xglobal, &x); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&sctr); CHKERRQ(ierr);
    ierr = VecDestroy(&xglobal); CHKERRQ(ierr);
    return 0;
}

int petsc_init(int *argc, char ***argv)
{
    PetscErrorCode ierr;
    ierr = PetscInitialize(argc, argv, NULL, NULL); CHKERRQ(ierr);
    return 0;
}

int petsc_load(app_t *ctx)
{
    int err;
    model_t *model = (model_t *) ctx->solverdata;

    err = app_init(ctx, &IO_PETSC, model->num_add_vars); if(err) return err;
    err = model->init(ctx); if(err) return err;
    return 0;
}

int petsc_finish(app_t *ctx)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) ctx->solverdata;
    ierr = model->finish(ctx); CHKERRQ(ierr);
    app_destroy(ctx);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}

int petsc_fopen(const char *filename, const char *mode, FILE **fp)
{
    PetscErrorCode ierr;
    ierr = PetscFOpen(PETSC_COMM_WORLD, filename, mode, fp); CHKERRQ(ierr);
    return 0;
}

int petsc_fprintf(FILE *fp, const char *fmt, ...)
{
    PetscErrorCode ierr;
    PetscInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(rank == 0) {
        va_list arg;
        va_start(arg, fmt);
        ierr = PetscVFPrintf(fp, fmt, arg); CHKERRQ(ierr);
        va_end(arg);
    }
    return 0;
}

int petsc_getline(FILE *fp, char **dst, size_t *len)
{
    PetscErrorCode ierr;
    ierr = PetscSynchronizedFGets(PETSC_COMM_WORLD, fp, *len, *dst);
        CHKERRQ(ierr);
    return 0;
}

int petsc_fclose(FILE *fp)
{
    PetscErrorCode ierr;
    ierr = PetscFClose(PETSC_COMM_WORLD, fp); CHKERRQ(ierr);
    return 0;
}
