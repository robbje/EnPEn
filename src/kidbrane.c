#include "app.h"
#include "iface_petsc.h"
#include "modules.h"
#include "parameters.h"
#include "mesh.h"
#include "simulation.h"
#include "equi.h"

PetscScalar u = 0.0001;    /*  Mean flow velocity dimless*/

int main(int argc, char **argv)
{
    int err;
    interface_t *iface = &IFACE_PETSC;

    err = iface->init(&argc, &argv); if(err) return err;

    /* Set up simulation */
    parameter_t p;
    mesh_t mesh;
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC_FLOW};
    param_norm_t norm = {.c0=1000,.D0=1.0e-9}; /* Molar scale (c0=1000) */
    param_set(&p, "sysdata/kidbrane.def", norm);


    meshdef_t def = TRILAYER;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t gel = TRILAYER_MEMB;
    gel.alpha = 0.2;
    gel.end_charge = -0.0040; // aem.charge*norm.c0 = [mol/m^3]
    gel.epsilon = 999999;       // Let's use this for MWC
    err = mesh_add_zone(&mesh, &gel); if(err) return err;
    
    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double blood[6] = {54.0541, 1e-7, 1e-7, 7e-4, 23*7e-4, 0.0};
    double urine[6] = {54.0541, 1e-7, 1e-7, 1e-04, 23e-04, 0.0};

    mesh_set_state(&mesh, &blood[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &blood[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &urine[0], n);

    copy_solution_from_mesh(&app);
    /* Run experiments */
    err = app_set_name(&app, "kidbrane/test"); if(err) return err;
    //err = sweep_boundary_chargeall_flow(&app, FACETYPE_BOUNDARY_LEFT, initial,
    //        &initial[p.num_species], 0, 40, 0.5, &u, 0.0, 0.01, 0.0001);
    PetscPrintf(PETSC_COMM_WORLD, "Calculate initial condition for ");
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_LEFT, &blood[0],
            &blood[p.num_species], 0, 0, 0.5); if(err) return err;
    //PetscPrintf(PETSC_COMM_WORLD, "Ramping flow:\n");
    //err = sweep_flow(&app, &u, 0.0, -0.01, -0.00001); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
