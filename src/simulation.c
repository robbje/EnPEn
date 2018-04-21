#include "simulation.h"

PetscErrorCode snes_solve(app_t *app)
{
    PetscErrorCode ierr;
    model_t *model = (model_t *) app->solverdata;
    ierr = copy_solution_from_mesh(app); CHKERRQ(ierr);
    ierr = SNESSolve(model->snes, NULL, model->x); CHKERRQ(ierr);
    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    /* Calculate fluxes */
    unsigned int i;
    for(i = 0; i < app->mesh->num_faces; i++) {
        ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
    }
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    ierr = mesh_dump(app->mesh, &IO_PETSC, 0, app->nv); CHKERRQ(ierr);
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode sweep_zone_adaptive(app_t *app, double *p, double start, double end)
{
    PetscErrorCode ierr;
    SNESConvergedReason reason;
    //PetscReal fnorm;
    model_t *model = (model_t *) app->solverdata;
    int is_positive;
    double old_p = *p;
    double step = (end-start)/10000.0;

    Vec tmp;
    ierr = VecDuplicate(model->x, &tmp); CHKERRQ(ierr);
    ierr = VecCopy(model->x, tmp); CHKERRQ(ierr);

    *p = start;
    PetscPrintf(PETSC_COMM_WORLD, "[");
    while(1) {
        /* Solve */
        mesh_apply_zones(app->mesh);
        ierr = SNESSolve(model->snes, NULL, model->x); CHKERRQ(ierr);
        ierr = SNESGetConvergedReason(model->snes, &reason); CHKERRQ(ierr);
        ierr = petsc_check_solution(model, &is_positive); CHKERRQ(ierr);
        if(is_positive && reason >= 0) {
            /* Accept step */
            old_p = *p;
            ierr = VecCopy(model->x, tmp); CHKERRQ(ierr);
            /* Choose new step */
            step *= 1.3;
            PetscPrintf(PETSC_COMM_WORLD, ".");
        } else {
            /* Revert step */
            *p = old_p;
            ierr = VecCopy(tmp, model->x); CHKERRQ(ierr);
            /* Choose new step */
            step *= 0.8;
            PetscPrintf(PETSC_COMM_WORLD, "<");
        }
        if(fabs(*p) >= fabs(end))
            break;
        if(fabs(step) <= 1e-09) {
            SETERRQ1(PETSC_COMM_WORLD, -1, "Stepsize too small: %g]\n", step);
        }

        /* Update step */
        if(fabs(*p+step) > fabs(end)) {
            *p = end;
        } else {
            *p += step;
        }
    }
    PetscReal res;
    ierr = SNESComputeFunction(model->snes, model->x, model->r); CHKERRQ(ierr);
    ierr = VecNorm(model->r, NORM_2, &res);
    PetscPrintf(PETSC_COMM_WORLD, "] %g\n", res);
    return 0;
}

PetscErrorCode sweep_zones(app_t *app)
{
    PetscErrorCode ierr;
    SNESConvergedReason reason;
    model_t *model = (model_t *) app->solverdata;
    int is_positive;
    double step = 1e-5; // Chosen by dice roll

    Vec tmp;
    ierr = VecDuplicate(model->x, &tmp); CHKERRQ(ierr);
    ierr = VecCopy(model->x, tmp); CHKERRQ(ierr);

    zonelist_t *zonelist = app->mesh->zones;
    if(zonelist == NULL) {
        PetscPrintf(PETSC_COMM_WORLD, "No zones found. Not sweeping.\n");
        ierr = SNESSolve(model->snes, NULL, model->x); CHKERRQ(ierr);
        ierr = SNESGetConvergedReason(model->snes, &reason); CHKERRQ(ierr);
        ierr = petsc_check_solution(model, &is_positive); CHKERRQ(ierr);
        PetscReal res;
        ierr = SNESComputeFunction(model->snes, model->x, model->r); CHKERRQ(ierr);
        ierr = VecNorm(model->r, NORM_2, &res);
        PetscPrintf(PETSC_COMM_WORLD, "[%g]: %s\n", res, SNESConvergedReasons[reason]);
        goto sweep_zones_out;
    }

    while(zonelist) {
        zonelist->zone->charge = 0.0;
        zonelist->zone->old_charge = 0.0;
        zonelist = zonelist->next;
    }

    PetscPrintf(PETSC_COMM_WORLD, "[");
    int done = 0;
    while(!done) {
        /* Solve */
        mesh_apply_zones(app->mesh);
        ierr = SNESSolve(model->snes, NULL, model->x); CHKERRQ(ierr);
        ierr = SNESGetConvergedReason(model->snes, &reason); CHKERRQ(ierr);
        ierr = petsc_check_solution(model, &is_positive); CHKERRQ(ierr);
        if(is_positive && reason >= 0) {
            /* Accept step */
            zonelist = app->mesh->zones;
            while(zonelist) {
                zonelist->zone->old_charge = zonelist->zone->charge;
                zonelist = zonelist->next;
            }
            ierr = VecCopy(model->x, tmp); CHKERRQ(ierr);
            /* Choose new step */
            step *= 1.5;
            PetscPrintf(PETSC_COMM_WORLD, ".");
        } else {
            /* Revert step */
            zonelist = app->mesh->zones;
            while(zonelist) {
                zonelist->zone->charge = zonelist->zone->old_charge;
                zonelist = zonelist->next;
            }
            ierr = VecCopy(tmp, model->x); CHKERRQ(ierr);
            /* Choose new step */
            step *= 0.8;
            if(is_positive) {
                PetscPrintf(PETSC_COMM_WORLD, "<");
            } else {
                PetscPrintf(PETSC_COMM_WORLD, "-");
            }
        }
        zonelist = app->mesh->zones;
        if(fabs(step) <= 1e-09) {
            SETERRQ1(PETSC_COMM_WORLD, -1, "Stepsize too small: %g]\n", step);
        }

        done = 1;
        /* Apply steps */
        zonelist = app->mesh->zones;
        while(zonelist) {
            if(fabs(zonelist->zone->charge) < fabs(zonelist->zone->end_charge)) {
                done = 0;
                if(zonelist->zone->charge < zonelist->zone->end_charge) {
                    zonelist->zone->charge += step;
                } else {
                    zonelist->zone->charge -= step;
                }
                if(fabs(zonelist->zone->charge) > fabs(zonelist->zone->end_charge)) {
                    zonelist->zone->charge = zonelist->zone->end_charge;
                }
            }
            zonelist = zonelist->next;
        }
    }
    PetscReal res;
    ierr = SNESComputeFunction(model->snes, model->x, model->r); CHKERRQ(ierr);
    ierr = VecNorm(model->r, NORM_2, &res);
    PetscPrintf(PETSC_COMM_WORLD, "] %g Reason: %s\n", res, SNESConvergedReasons[reason]);
sweep_zones_out:
    ierr = VecDestroy(&tmp); if(ierr) return ierr;
    return 0;
}

PetscErrorCode sweep_boundary_with_charge(app_t *app, int bc, double *initial,
        double *p, double start, double end, double step, zone_t *z)
{
    PetscErrorCode ierr;
    Vec init;
    model_t *model = (model_t *) app->solverdata;
    ierr = VecDuplicate(model->x, &init);
    ierr = VecCopy(model->x, init); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    for(*p = start; *p <= end; *p += step) {
        ierr = VecCopy(init, model->x); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "V=%g ", *p);
        boundary_set_state(app->mesh, bc, initial, app->nv);
        ierr = sweep_zone_adaptive(app, &(z->charge), 0, -1); CHKERRQ(ierr);
        ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
        /* Calculate fluxes */
        unsigned int i;
        for(i = 0; i < app->mesh->num_faces; i++) {
            ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
        }
        ierr = mesh_dump(app->mesh, &IO_PETSC, *p, app->nv); CHKERRQ(ierr);
    }
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    ierr = VecDestroy(&init); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode sweep_boundary_chargeall(app_t *app, int bc, double *initial,
    double *p, double start, double end, double step)
{
    PetscErrorCode ierr;
    Vec init;
    model_t *model = (model_t *) app->solverdata;
    ierr = VecDuplicate(model->x, &init); CHKERRQ(ierr);
    ierr = VecCopy(model->x, init); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    for(*p = start; *p <= end; *p += step) {
        ierr = VecCopy(init, model->x); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "V=%g ", *p);
        boundary_set_state(app->mesh, bc, initial, app->nv);
        ierr = sweep_zones(app); CHKERRQ(ierr);
        ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
        /* Calculate fluxes */
        unsigned int i;
        for(i = 0; i < app->mesh->num_faces; i++) {
            ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
        }
        ierr = mesh_dump(app->mesh, &IO_PETSC, *p, app->nv); CHKERRQ(ierr);
    }
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    ierr = VecDestroy(&init); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode sweep_flow(app_t *app, double *p, double start, double end, double step)
{
    PetscErrorCode ierr;
    Vec init;
    SNESConvergedReason reason;
    int is_positive;
    model_t *model = (model_t *) app->solverdata;
    ierr = VecDuplicate(model->x, &init); CHKERRQ(ierr);
    ierr = VecCopy(model->x, init); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "%g [", start);
    for(*p = start; fabs(*p) <= fabs(end); *p += step) {
        ierr = SNESSolve(model->snes, NULL, model->x); CHKERRQ(ierr);
        ierr = SNESGetConvergedReason(model->snes, &reason); CHKERRQ(ierr);
        ierr = petsc_check_solution(model, &is_positive); CHKERRQ(ierr);
        if(is_positive && reason >= 0) {
            PetscPrintf(PETSC_COMM_WORLD, ".");
        } else {
            break;
        }
        if((int)(*p/step) % 10 == 0) {
            ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
            /* Calculate fluxes */
            unsigned int i;
            for(i = 0; i < app->mesh->num_faces; i++) {
                ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
            }
            ierr = mesh_dump(app->mesh, &IO_PETSC, (*p/step), app->nv); CHKERRQ(ierr);
        }
    }
    PetscPrintf(PETSC_COMM_WORLD, "] %g\n", *p-step);
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    ierr = VecDestroy(&init); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode sweep_boundary_chargeall_flow(app_t *app, int bc, double *initial,
    double *p, double start, double end, double step, double *u, double u_start,
    double u_end, double u_step)
{
    PetscErrorCode ierr;
    SNESConvergedReason reason;
    int is_positive, i;
    Vec init;
    model_t *model = (model_t *) app->solverdata;
    ierr = VecDuplicate(model->x, &init); CHKERRQ(ierr);
    ierr = VecCopy(model->x, init); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    for(*p = start; fabs(*p) <= fabs(end); *p += step) {
        *u = u_start;
        ierr = VecCopy(init, model->x); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "V=%g ", *p);
        boundary_set_state(app->mesh, bc, initial, app->nv);
        ierr = sweep_zones(app); CHKERRQ(ierr);

        PetscPrintf(PETSC_COMM_WORLD, "\tu [");
        for (*u = u_start, i = 0; fabs(*u) <= fabs(u_end); *u += u_step, i++) {
            ierr = SNESSolve(model->snes, NULL, model->x); CHKERRQ(ierr);
            ierr = SNESGetConvergedReason(model->snes, &reason); CHKERRQ(ierr);
            ierr = petsc_check_solution(model, &is_positive); CHKERRQ(ierr);
            if(is_positive && reason >= 0) {
                if(i % 10 == 0)
                    PetscPrintf(PETSC_COMM_WORLD, ".");
            } else {
                PetscPrintf(PETSC_COMM_WORLD, "x");
                break;
            }
        }
        PetscPrintf(PETSC_COMM_WORLD, "] %g\n", *u-u_step);
        ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
        /* Calculate fluxes */
        unsigned int i;
        for(i = 0; i < app->mesh->num_faces; i++) {
            ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
        }
        ierr = mesh_dump(app->mesh, &IO_PETSC, *p, app->nv); CHKERRQ(ierr);
    }
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    ierr = VecDestroy(&init); CHKERRQ(ierr);
    return 0;
}

PetscScalar last_R = 0.0;
PetscScalar last_dR = 0.0;
PetscScalar maxrel = 0.0;
PetscErrorCode ts_controller(TS ts)
{
    /* Here, we abuse p->u as the target current */
    PetscErrorCode ierr;
    app_t *app;
    ierr = TSGetApplicationContext(ts, &app); CHKERRQ(ierr);
    model_t *model = (model_t *) app->solverdata;

    unsigned int i;
    for(i = 0; i < app->mesh->num_faces; i++) {
        ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
    }

    PetscReal time, dt;
    ierr = TSGetTime(ts, &time); CHKERRQ(ierr);
    ierr = TSGetTimeStep(ts, &dt); CHKERRQ(ierr);
    parameter_t *p = app->p;
    PetscReal current = app->mesh->faces[0]->flux[p->num_species];
    volumelist_t *head = app->mesh->boundaries[FACETYPE_BOUNDARY_RIGHT];

    PetscScalar err = (p->u - current);
    PetscScalar old_bias = head->volume->state[app->p->num_species];
    if(old_bias == 0) old_bias = -0.0725;

    if(current < 1e-6) current = 1e-6;

    PetscScalar cur_R = fabs(old_bias/current);
    if(last_R == 0) last_R = cur_R;
    PetscScalar cur_dR = (cur_R - last_R)/(dt * p->t0);
    PetscScalar cur_ddR = (last_dR - cur_dR)/(dt * p->t0);
    last_dR = cur_dR;
    PetscScalar new_bias = old_bias - 1.4*cur_R * err - 0.01 * (dt * p->t0) * (cur_dR) * current;

    if(time*p->t0 > 1.0) {
        if(-0.01 * cur_ddR * current > 0.4)
            new_bias += 0.0123 * cur_ddR * current;
    }

    if(time*p->t0 > 0.2 && fabs(err/p->u) > maxrel) maxrel = fabs(err/p->u);


    PetscPrintf(PETSC_COMM_WORLD, "[%0.4fs] RELDIFF %g, new bias %g, MAXRELDIFF: %g\n", time*p->t0, err/p->u, new_bias, maxrel);
    boundary_modify_state(app->mesh, FACETYPE_BOUNDARY_RIGHT, new_bias, app->p->num_species);
    old_bias = new_bias;
    return 0;
};

PetscErrorCode ts_setter(TS ts)
{
    PetscErrorCode ierr;
    PetscReal time;
    app_t *app;
    ierr = TSGetTime(ts, &time); CHKERRQ(ierr);
    ierr = TSGetApplicationContext(ts, &app); CHKERRQ(ierr);
    PetscReal final_bias = 15.0;
    PetscReal tau = 10/app->p->t0; // t0 = 0.1

    PetscReal bias = final_bias*(1-exp(-time/tau));
    PetscPrintf(PETSC_COMM_WORLD, "%gs: %g\n", time*app->p->t0, bias);
    
    boundary_modify_state(app->mesh, FACETYPE_BOUNDARY_RIGHT, bias, app->p->num_species);
    return 0;
}

PetscErrorCode time_ivc(app_t *app, double *initial, unsigned int n)
{
    PetscErrorCode ierr;

    mesh_set_state(app->mesh, initial, n);
    boundary_set_state(app->mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(app->mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    copy_solution_from_mesh(app);
    printf("Initializing ");
    ierr = sweep_zones(app); CHKERRQ(ierr);
    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    boundary_modify_state(app->mesh, FACETYPE_BOUNDARY_LEFT, 0.0, app->p->num_species);

    model_t *model = (model_t *) app->solverdata;
    double final_time = 80.0/app->p->t0; /* In seconds */
    /* In seconds: If you encounter negative concentrations, lower this value.
     * For reactive systems I found 0.1/app->p->t0 to be a good value.
     */
    double dt = 0.5/app->p->t0;
    TSType t;
    ierr = TSGetType(model->ts, &t);
    PetscPrintf(PETSC_COMM_WORLD, "%s\n", t);
    ierr = TSSetTime(model->ts, 0.0); CHKERRQ(ierr);
    ierr = TSSetTimeStep(model->ts, dt); CHKERRQ(ierr);
    ierr = TSSetSolution(model->ts, model->x); CHKERRQ(ierr);
    ierr = TSSetMaxSteps(model->ts, final_time/dt); CHKERRQ(ierr);
    ierr = TSSetMaxTime(model->ts, final_time); CHKERRQ(ierr);
    ierr = TSMonitorSet(model->ts, tsmon, (void *) app, NULL); CHKERRQ(ierr);
    ierr = TSSetApplicationContext(model->ts, (void *) app); CHKERRQ(ierr);
    ierr = TSSetPreStep(model->ts, ts_setter); CHKERRQ(ierr);
    ierr = TSSetUp(model->ts); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    ierr = TSSolve(model->ts, model->x); CHKERRQ(ierr);
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode time_current_control(app_t *app, double *initial, unsigned int n)
{
    PetscErrorCode ierr;
    mesh_set_state(app->mesh, initial, n);
    boundary_set_state(app->mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(app->mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    copy_solution_from_mesh(app);
    ierr = sweep_zones(app); CHKERRQ(ierr);
    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    boundary_modify_state(app->mesh, FACETYPE_BOUNDARY_LEFT, 0.0, app->p->num_species);
    model_t *model = (model_t *) app->solverdata;
    double final_time = 30.0/app->p->t0; /* In seconds */
    double dt = 0.1/app->p->t0;
    TSType t;
    ierr = TSGetType(model->ts, &t);
    PetscPrintf(PETSC_COMM_WORLD, "%s\n", t);
    ierr = TSSetTime(model->ts, 0.0); CHKERRQ(ierr);
    ierr = TSSetTimeStep(model->ts, dt); CHKERRQ(ierr);
    ierr = TSSetSolution(model->ts, model->x); CHKERRQ(ierr);
    ierr = TSSetMaxSteps(model->ts, final_time/dt); CHKERRQ(ierr);
    ierr = TSSetMaxTime(model->ts, final_time); CHKERRQ(ierr);
    ierr = TSMonitorSet(model->ts, tsmon, (void *) app, NULL); CHKERRQ(ierr);
    ierr = TSSetApplicationContext(model->ts, (void *) app); CHKERRQ(ierr);
    ierr = TSSetPreStep(model->ts, ts_controller); CHKERRQ(ierr);
    ierr = TSSetUp(model->ts); CHKERRQ(ierr);
    ierr = mesh_open_dumpfile(&IO_PETSC, app->fname); CHKERRQ(ierr);
    ierr = app_dump_metadata(app, &IO_PETSC); CHKERRQ(ierr);
    ierr = TSSolve(model->ts, model->x); CHKERRQ(ierr);
    ierr = mesh_close_dumpfile(&IO_PETSC); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode tsmon(TS ts, PetscInt steps, PetscReal time, Vec u, void *context)
{
    PetscErrorCode ierr;
    app_t *app = (app_t *) context;
    model_t *model = (model_t *) app->solverdata;
    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    unsigned int i;
    for(i = 0; i < app->mesh->num_faces; i++) {
        ierr = model->fnc_flux(app->mesh->faces[i], app); CHKERRQ(ierr);
    }
    ierr = mesh_dump(app->mesh, &IO_PETSC, time*app->p->t0, app->nv); CHKERRQ(ierr);
    return 0;
}
