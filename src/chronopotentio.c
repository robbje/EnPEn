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
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC_TS};
    param_norm_t norm = {.c0=1000,.D0=1.0e-9}; /* Molar scale (c0=1000) */
    param_set(&p, "sysdata/h2o_nacl.def", norm);

    meshdef_t def = TRILAYER;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t aem = TRILAYER_MEMB;
    aem.alpha = 0.1;
    aem.end_charge = 1; // aem.charge*norm.c0 = [mol/m^3]
    aem.epsilon = 1;
    err = mesh_add_zone(&mesh, &aem); if(err) return err;
    
    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double initial[6] = {54, 1e-7, 1e-7, 1e-3, 1e-3, 0.0};

    p.impedance.bias = 0;
    p.impedance.modulus = 1;
    p.impedance.f = 1000;

    mesh_set_state(&mesh, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);
    copy_solution_from_mesh(&app);
  
    PetscScalar V_final = 10;
    /* Find initial solution */
    printf("Calculating final solution:\n");
    initial[n-1] = V_final;
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);
    err = app_set_name(&app, "ts/final"); if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 0, 0.5); if(err) return err;

    intial[n-1] = 0;
    printf("Calculating initial solution:\n");
    err = app_set_name(&app, "ts/initial"); if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 0, 0.5); if(err) return err;
    copy_solution_to_mesh(&app);

    /* Run experiments */
    initial[n-1] = V_final;
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);

    printf("Starting simulation...\n");
    copy_solution_from_mesh(&app);
    time_ivc(&app, 1e-4, 5);
    /* Finish and exit */
    return iface->finish(&app);
}
