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
    param_set(&p, "sysdata/nacl.def", norm);
    meshdef_t def = MAARTEN;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t aem = MAARTEN_MEMB;
    aem.alpha = 0.1;
    aem.end_charge = -1; // aem.charge*norm.c0 = [mol/m^3]
    aem.epsilon = 1;
    err = mesh_add_zone(&mesh, &aem); if(err) return err;

    /* Load application */
    err = iface->load(&app); if(err) return err;

    /* Setup boundary condition */
    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double initial[3] = {1e-3, 1e-3, 0.0};
    mesh_set_state(&mesh, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);

    /* Adjust diffusion coefficient */
    unsigned int i;
    for(i = 0; i < mesh.num_volumes; i++) {
        volume_t *v = mesh.volumes[i];
        if (v->center.x < 100.0) {
            v->epsilon = exp((100.0 - v->center.x)/10.0);
        }
        if (v->center.x > 110.0) {
            v->epsilon = exp((v->center.x-110)/10.0);
        }
    }
    copy_solution_from_mesh(&app);
  
    /* Run experiments */
    err = app_set_name(&app, "test"); if(err) return err;
    err = time_ivc(&app, &initial[0], n); if(err) return err;
    //err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
    //        &initial[p.num_species], 0, 50, 0.5); if(err) return err;

    /* Finish and exit */
    return iface->finish(&app);
}
