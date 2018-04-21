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
    param_norm_t norm = {.c0=1000,.D0=1e-9}; /* Molar scale (c0=1000) */
    param_set(&p, "sysdata/h2o_real.def", norm);

    meshdef_t def = BPM;
    mesh_set(&mesh, &def);

    /* Set up bpm */
    zone_t aem = BPM_AEM;
    aem.alpha = 0.1;
    aem.end_charge = 0.5; // aem.charge*norm.c0 = [mol/m^3]
    aem.epsilon = 1;
    err = mesh_add_zone(&mesh, &aem); if(err) return err;

    zone_t cem = BPM_CEM;
    cem.alpha = 0.1;
    cem.end_charge = -0.5;
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;
    
    /* Load application */
    err = iface->load(&app); if(err) return err;
    unsigned int n = p.num_species + NPP_1D.num_add_vars;

    /* Find initial conditions */
    double *initial = alloca(app.nv*sizeof(double));
    equi_t *e;
    err = equi_create(&e); if(err) return err;
    err = equi_parse(e, &p); if(err) return err;
    e->component[0]->c = 50/130.1;
    err = equi_solve(initial, &p); if(err<0) return err;
    equi_get_initial(e, initial);
    equi_destroy(e);

    mesh_set_state(&mesh, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    copy_solution_from_mesh(&app);
  
    /* Run experiments */
    err = app_set_name(&app, "outdata/bpm");
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 10, 1); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
