#include "app.h"
#include "iface_petsc.h"
#include "modules.h"
#include "parameters.h"
#include "mesh.h"
#include "simulation.h"
#include "equi.h"

int main(int argc, char **argv)
{
    int err;
    interface_t *iface = &IFACE_PETSC;
    err = iface->init(&argc, &argv); if(err) return err;

    /* Set up simulation */
    parameter_t p;
    mesh_t mesh;
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC};
    param_norm_t norm = {.c0=1000,.D0=1.0e-9}; /* Molar scale (c0=1000) */
    param_set(&p, "sysdata/nacl.def", norm);

    meshdef_t def = KHAIR;
    mesh_set(&mesh, &def);

    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D_REAC.num_add_vars;
    double initial[3] = {1e-3, 1e-3, 0.5};
    mesh_set_state(&mesh, &initial[0], n);
    //mesh.volumes[mesh.num_volumes-1]->state[0] -= 1e-4;
    //mesh.volumes[mesh.num_volumes-1]->state[1] += 1e-4;
    double lp = -1.0;
    double rp = 1.0;
    //unsigned int i;
    //for(i = 0; i < mesh.num_volumes; i++) {
    //    double x = mesh.volumes[i]->center.x;
    //    double potential = lp+0.25 + (rp - lp - 0.5)/100.0 * x;
    //    mesh.volumes[i]->state[p.num_species] = potential;
    //}

    boundary_set_type(&mesh, FACETYPE_BOUNDARY_LEFT, VOLTYPE_ELECTRODE);
    boundary_set_type(&mesh, FACETYPE_BOUNDARY_RIGHT, VOLTYPE_ELECTRODE);
    initial[2] = lp;
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &initial[0], n);
    initial[2] = rp;
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);
    copy_solution_from_mesh(&app);
  
    /* Run experiments */
    err = app_set_name(&app, "khair"); if(err) return err;
    //err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
    //        &initial[p.num_species], 1.0, 15, 0.5); if(err) return err;
    err = snes_solve(&app); if(err) return err;
    model_t *model = (model_t *) app.solverdata;
    PetscReal res;
    SNESConvergedReason reason;
    SNESComputeFunction(model->snes, model->x, model->r);
    SNESGetConvergedReason(model->snes, &reason);
    VecNorm(model->r, NORM_2, &res);
    PetscPrintf(PETSC_COMM_WORLD, "[%s] %g\n", SNESConvergedReasons[reason], res);

    /* Finish and exit */
    return iface->finish(&app);
}
