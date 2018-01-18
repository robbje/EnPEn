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
    param_set(&p, "sysdata/joeri.def", norm);

    meshdef_t def = JOERI;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t cem = JOERI_CEM;
    cem.alpha = 0.1;
    cem.end_charge = -1.0; // aem.charge*norm.c0 = [mol/m^3]
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

    /* Set up a membrane */
    zone_t cem2 = JOERI_AEM;
    cem2.alpha = 0.1;
    cem2.end_charge = -1.0; // aem.charge*norm.c0 = [mol/m^3]
    cem2.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem2); if(err) return err;
    
    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;

    /* Find initial conditions */
    //double *initial = alloca(app.nv*sizeof(double));
    //equi_t *e;
    //err = equi_create(&e); if(err) return err;
    //err = equi_parse(e, &p); if(err) return err;
    //e->component[0]->c = 0.1;
    //equi_get_initial(e, initial);
    //equi_destroy(e);
    //double initial2[9] = {54.0541,      // H2O
    //                      0.000031765,  // H+
    //                      4.34e-10,     // OH-
    //                      0.0068994,    // Cu2+
    //                      0.0069138,    // SO42-
    //                      0.0030687,    // CuSO4
    //                      0.000011409,  // HSO4-
    //                      9.36e-6-9.15565e-7-1.26906e-12,  // CuHSO4
    //                      0.0};         // el. potential
    double initial2[3] = {1e-2, 1e-2, 0.0};
    double *initial = &initial2[0];

    //err = equi_solve(initial, &p); if(err<0) return err;
    //mesh_set_state(&mesh, initial, n);
    //boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    //boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    //copy_solution_from_mesh(&app);
  
    /* Run experiments */
    err = app_set_name(&app, "joeri/doubleCEM_ts"); if(err) return err;
    //err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_LEFT, initial,
    //        &initial[p.num_species], 0, 10, 0.5); if(err) return err;
    err = time_ivc(&app, initial, n); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
