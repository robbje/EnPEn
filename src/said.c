#include "app.h"
#include "iface_petsc.h"
#include "modules.h"
#include "parameters.h"
#include "mesh.h"
#include "simulation.h"
#include "equi.h"

double P_sodium(face_t *face, double *initial)
{
    return initial[1]/initial[0] * face->flux[0] / (2*face->flux[1]);
}

int main(int argc, char **argv)
{
    int err;
    interface_t *iface = &IFACE_PETSC;

    err = iface->init(&argc, &argv); if(err) return err;

    /* Set up simulation */
    parameter_t p;
    mesh_t mesh;
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC};
    param_norm_t norm = {.c0=1000,.D0=1e-9}; /* Molar scale (c0=1000) */
    //param_set(&p, "sysdata/h2o_said.def", norm);
    param_set(&p, "sysdata/said.def", norm);

    meshdef_t def = SAID;
    mesh_set(&mesh, &def);

    /* Load application */
    err = iface->load(&app); if(err) return err;
    unsigned int n = p.num_species + NPP_1D.num_add_vars;

    /* Initial conditions */
    //double initial2[7] = {54.0541, 1e-7, 1e-7, 0.005, 0.0025, 0.010, 0};
    double initial2[4] = {1e-3, 1e-3, 3e-3, 0};
    double *initial = &initial2[0];

    /* Settings */
    /* Try three systems */
    zone_t cem = SAID_CEM;
    cem.alpha = 0.1;
    cem.end_charge = -1.0;
    cem.epsilon = 1.0;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

    //zone_t aem = SAID_AEM;
    //aem.alpha = 0.1;
    //aem.end_charge = 0.1;
    //aem.epsilon = 1.0;
    //err = mesh_add_zone(&mesh, &aem); if(err) return err;

    //zone_t cem2 = SAID_CEM2;
    //cem2.alpha = 0.1;
    //cem2.end_charge = -0.1;
    //cem2.epsilon = 1.0;
    //err = mesh_add_zone(&mesh, &cem2); if(err) return err;

    //zone_t aem2 = SAID_AEM2;
    //aem2.alpha = 0.1;
    //aem2.end_charge = 0.1;
    //aem2.epsilon = 1.0;
    //err = mesh_add_zone(&mesh, &aem2); if(err) return err;

    //zone_t cem3 = SAID_CEM3;
    //cem3.alpha = 0.1;
    //cem3.end_charge = -0.1;
    //cem3.epsilon = 1.0;
    //err = mesh_add_zone(&mesh, &cem3); if(err) return err;

    mesh_set_state(&mesh, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    copy_solution_from_mesh(&app);

    err = app_set_name(&app, "said");
    if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 10, 0.5); if(err) return err;

    /* Finish and exit */
    return iface->finish(&app);
}
