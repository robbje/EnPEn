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
    //param_set(&p, "sysdata/h2o_h2co3.def", norm);
    //param_set(&p, "sysdata/nacl.def", norm);
    //param_set(&p, "sysdata/h2o_nacl.def", norm);

    meshdef_t def = TRILAYER;
    //meshdef_t def = MAARTEN;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t cem = TRILAYER_MEMB;
    //zone_t cem = MAARTEN_MEMB;
    cem.alpha = 0.1;
    cem.end_charge = +1.0; // aem.charge*norm.c0 = [mol/m^3]
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double initial[3] = {1e-3,1e-3,0};
    //double initial[6] = {54.0, 1e-7, 1e-7, 1e-3, 1e-3, 0.0};
    //double initial[7] = {54.0, 1e-7, 1e-7, 1e-3, 1e-3, 1e-3};
    //double initial[9] = {54, 1e-7, 1e-7, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0.0};
    //equi_t *e;
    //err = equi_create(&e); if(err) return err;
    //err = equi_parse(e, &p); if(err) return err;
    //e->component[0]->c = 5e-3;
    //equi_get_initial(e, initial);
    //equi_destroy(e);
    //unsigned int i;
    //for(i = 0; i < p.num_species; i++)
    //    PetscPrintf(PETSC_COMM_WORLD, "[%s] %g\n", p.name[i], initial[i]);
    PetscReal lambda = sqrt(p.epsilon*p.R*p.T/(2.0*initial[0]*p.c0*pow(p.F, 2.0)));
    PetscPrintf(PETSC_COMM_WORLD, "Debye length: %gm\n", lambda);
    mesh_set_state(&mesh, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);
    /* Adjust diffusion coefficient */
    //unsigned int i;
    //for(i = 0; i < mesh.num_volumes; i++) {
    //    volume_t *v = mesh.volumes[i];
    //    if (v->center.x < 100.0) {
    //        v->epsilon = exp((100.0 - v->center.x)/10.0);
    //    }
    //    if (v->center.x > 110.0) {
    //        v->epsilon = exp((v->center.x-110)/10.0);
    //    }
    //}
    copy_solution_from_mesh(&app);
  
    /* Run experiments */
    err = app_set_name(&app, "test");
    if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 20, 0.5); if(err) return err;

    //err = time_ivc(&app, initial, n); if(err) return err;
    //err = time_current_control(&app, initial, n); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
