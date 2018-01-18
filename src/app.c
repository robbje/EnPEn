#include "app.h"
#include <stdio.h>

int app_init(app_t *app, io_t *io, unsigned int num_add_vars)
{
    int err;
    /* Load parameters */
    err = param_load(app->p, io); if(err) return err;
    err = io->fclose(io->fh); if(err) return err;

    app->nv = app->p->num_species + num_add_vars;
    err = asprintf(&(app->fname), "outdata/sweep.json");
    if(err == -1) return err;

    /* Load mesh */
    err = mesh_parse(app->mesh, io);
    if(err) {
        return err;
    }
    err = io->fclose(io->fh); if(err) return err;

    /* Initialize state space */
    err = mesh_init_states(app->mesh, app->nv);
    if(err) return err;

    /* Normalization and constants */
    app->p->debye = pow(app->p->F, 2.0) * app->p->c0 *
        pow(app->mesh->def->mesh_unit, 2.0) /
        (app->p->R * app->p->T * app->p->epsilon);

    app->p->mu0 = pow(app->mesh->def->mesh_unit, 2.0) * app->p->c0 *
                    app->p->R * app->p->T / app->p->D0;
    app->p->L0 = app->mesh->def->mesh_unit;
    app->p->phi0 = app->p->R * app->p->T / app->p->F;
    app->p->p0 = app->p->c0 * app->p->R * app->p->T;
    app->p->u0 = app->p->D0/app->p->L0;
    app->p->t0 = pow(app->p->L0, 2.0)/app->p->D0;

    unsigned int i, j;
    for(i = 0; i < app->p->num_species; i++) {
        app->p->D[i] /= app->p->D0;
    }

    for(i = 0; i < app->p->num_reactions; i++) {
        int mu_sum = 0, nu_sum = 0;
        for(j = 0; j < app->p->num_species; j++) {
            mu_sum += app->p->mu[j][i];
            nu_sum += app->p->nu[j][i];
        }
        app->p->k_f[i] *= pow(app->mesh->def->mesh_unit, 2.0) / app->p->D0;
        app->p->k_b[i] *= pow(app->mesh->def->mesh_unit, 2.0) / app->p->D0;
    }
    return 0;
}

void app_destroy(app_t *app)
{
    mesh_destroy(app->mesh);
    param_destroy(app->p);
}

int app_dump_metadata(app_t *app, io_t *io)
{
    int err;
    err = io->fprintf(io->fh, "\"meta\":{"); if(err) return err;
    err = io->fprintf(io->fh, "\"nv\":%i,", app->nv); if(err) return err;
    err = io->fprintf(io->fh, "\"species\":["); if (err) return err;
    unsigned int i;
    for(i = 0; i < app->p->num_species; i++) {
        err = io->fprintf(io->fh, "\"%s\",", app->p->name[i]); if(err) return err;
    }
    err = io->fprintf(io->fh, "\"phi\"],\"diffcoeff\":["); if(err) return err;
    for(i = 0; i < app->p->num_species; i++) {
        err = io->fprintf(io->fh, "%.12f,",app->p->D[i]*app->p->D0);
        if(err) return err;
    }
    err = io->fprintf(io->fh, "0],\"charge\":["); if(err) return err;
    for(i = 0; i < app->p->num_species; i++) {
        err = io->fprintf(io->fh, "%i,", app->p->z[i]); if(err) return err;
    }
    err = io->fprintf(io->fh, "0]"); if(err) return err;
    volumelist_t *head = app->mesh->boundaries[FACETYPE_BOUNDARY_LEFT];
    if(head && head->volume && head->volume->state) {
        err = io->fprintf(io->fh, ",\"left_boundary\":["); if(err) return err;
        for(i = 0; i < app->p->num_species; i++) {
            err = io->fprintf(io->fh, "%.2f,", head->volume->state[i]); if(err) return err;
        }
        err = io->fprintf(io->fh, "%.2f]", head->volume->state[app->p->num_species]); if(err) return err;
    }
    head = app->mesh->boundaries[FACETYPE_BOUNDARY_RIGHT];
    if(head && head->volume && head->volume->state) {
        err = io->fprintf(io->fh, ",\"right_boundary\":["); if(err) return err;
        for(i = 0; i < app->p->num_species; i++) {
            err = io->fprintf(io->fh, "%.2f,", head->volume->state[i]); if(err) return err;
        }
        err = io->fprintf(io->fh, "%.2f]", head->volume->state[app->p->num_species]); if(err) return err;
    }
    err = io->fprintf(io->fh, "},\n"); if(err) return err;
    return 0;
}

int app_set_name(app_t *app, const char *fname)
{
    free(app->fname);
    int err = asprintf(&(app->fname), "outdata/%s.json", fname);
    if(err == -1) {
        perror("asprintf");
        return err;
    }
    return 0;
}
