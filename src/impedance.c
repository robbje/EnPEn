#include "impedance.h"
#define SAMPLES_PER_PERIOD 50
#define NUM_PERIODS 15
#define NUM_IGNORE_PERIODS 10
#define NUM_LOCS 3

void impedance_set_boundaries(app_t *app, PetscReal time)
{
    parameter_t *p = app->p;

    PetscReal V = 0.5 * (p->impedance.bias +
        p->impedance.modulus * sin(p->impedance.f*TAU*time*p->t0));

    boundary_modify_state(app->mesh, FACETYPE_BOUNDARY_RIGHT, V, p->num_species);
    boundary_modify_state(app->mesh, FACETYPE_BOUNDARY_LEFT, -V, p->num_species);
}

PetscReal impedance_get_field(volume_t *v, unsigned int n)
{
    face_t *f = v->face[0];
    volume_t *nv = f->nvol[0] == v?f->nvol[1]:f->nvol[0];
    PetscReal delta = v->distance[0];
    delta += nv->face[0] == f ? nv->distance[0] : nv->distance[1];
    PetscReal E = (nv->state[n] - v->state[n])/delta;
    return (nv->center.x > v->center.x)? E : -E;
}

PetscReal impedance_last_sample = 0;
PetscReal E_last = 0.0;
PetscErrorCode impedance_sampler(TS ts, PetscInt steps, PetscReal time, Vec u,
                                 void *ctx)
{
    if (time == 0) return 0;
    PetscErrorCode ierr;
    app_t *app = (app_t *) ctx;
    model_t *model = (model_t *) app->solverdata;
    parameter_t *p = app->p;
    PetscReal T = 1/p->impedance.f;
    PetscReal dt = (1.0/app->p->impedance.f)/SAMPLES_PER_PERIOD/app->p->t0;
    PetscReal V = (p->impedance.bias +
        p->impedance.modulus * sin(p->impedance.f*TAU*(time-dt)*p->t0));

    // Record last 5 periods
    if(time*p->t0 < NUM_IGNORE_PERIODS*T-dt) return 0;

    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    unsigned int i;
    for(i = 0; i < app->mesh->num_faces; i++) {
        ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
    }
    io_t *io = &IO_PETSC;
    
    PetscReal j = -1.0*app->mesh->volumes[0]->face[1]->flux[p->num_species];
    PetscReal E = impedance_get_field(app->mesh->volumes[0], p->num_species);
    E *= p->phi0/p->L0;
    PetscReal j_D = p->epsilon / (dt*p->t0) * (E - E_last);
    E_last = E;

    io->fprintf(io->fh, "%g, %g, %g", time*p->t0, V*p->phi0, j*p->F+j_D);
    //for(i = 0; i < app->mesh->num_volumes; i++) {
    //    PetscReal eta = p->R * p->T * log(app->mesh->volumes[i]->state[0]);
    //    eta += p->F * app->mesh->volumes[i]->state[p->num_species] * p->phi0;
    //    io->fprintf(io->fh, ",%g", eta);
    //}

    io->fprintf(io->fh, "\n");

    return 0;
}

PetscErrorCode impedance_setter(TS ts)
{
    PetscErrorCode ierr;
    PetscReal time;
    app_t *app;
    ierr = TSGetTime(ts, &time); CHKERRQ(ierr);
    ierr = TSGetApplicationContext(ts, &app); CHKERRQ(ierr);
    impedance_set_boundaries(app, time);
    return 0;
}

PetscErrorCode impedance_setup(app_t *app, double *initial, unsigned int n)
{
    PetscErrorCode ierr;
    //model_t *model = (model_t *) app->solverdata;
    mesh_set_state(app->mesh, initial, n);
    boundary_set_state(app->mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(app->mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    impedance_set_boundaries(app, 0.0);
    copy_solution_from_mesh(app);

    printf("Bias=%g ", app->p->impedance.bias);
    ierr = sweep_zones(app); CHKERRQ(ierr);
    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    //ierr = app_set_name(app, "ts/initial"); CHKERRQ(ierr);
    //ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    //ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    //ierr = mesh_dump(app->mesh, &IO_PETSC, 0.0, model->bs); CHKERRQ(ierr);
    //ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode impedance_run(app_t *app)
{
    PetscErrorCode ierr;
    PetscReal dt, t_final;
    t_final = (1.0/app->p->impedance.f)*NUM_PERIODS/app->p->t0;
    dt = (1.0/app->p->impedance.f)/SAMPLES_PER_PERIOD/app->p->t0;
    model_t *model = (model_t *) app->solverdata;

    ierr = TSSetInitialTimeStep(model->ts, 0.0, dt); CHKERRQ(ierr);
    ierr = TSSetSolution(model->ts, model->x); CHKERRQ(ierr);
    ierr = TSSetDuration(model->ts, 2*t_final/dt, t_final+dt); CHKERRQ(ierr);
    ierr = TSMonitorSet(model->ts, impedance_sampler, (void *) app, NULL);
    CHKERRQ(ierr);
    ierr = TSSetApplicationContext(model->ts, (void *) app); CHKERRQ(ierr);
    ierr = TSSetPreStep(model->ts, impedance_setter); CHKERRQ(ierr);
    ierr = TSSetUp(model->ts); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    io_t *io = &IO_PETSC;
    io->fprintf(io->fh, "\"imp\":{\"f\":%g,\"bias\":%g,\"modulus\":%g},",
                app->p->impedance.f,
                app->p->impedance.bias,
                app->p->impedance.modulus);
    io->fprintf(io->fh, "\"x\":[");
    int i;
    for(i = 0; i < app->mesh->num_volumes; i++) {
        if(i == app->mesh->num_volumes-1) {
            io->fprintf(io->fh, "%g", app->mesh->volumes[i]->center.x);
        } else {
            io->fprintf(io->fh, "%g,", app->mesh->volumes[i]->center.x);
        }
    }
    io->fprintf(io->fh, "]}\n");
    io->fprintf(io->fh, "t, V, i\n");
    ierr = TSSolve(model->ts, model->x); CHKERRQ(ierr);
    TSConvergedReason reason;
    ierr = TSGetConvergedReason(model->ts, &reason); CHKERRQ(ierr);
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    if(reason != TS_CONVERGED_TIME) return -1;
    return 0;
}
