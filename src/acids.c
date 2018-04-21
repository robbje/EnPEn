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
    //param_set(&p, "sysdata/h2o_ita.def", norm);
    param_set(&p, "sysdata/h2o_oxal.def", norm);

    meshdef_t def = TRILAYER;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t cem = TRILAYER_MEMB;
    cem.alpha = 0.1;
    cem.end_charge = 1.0; // aem.charge*norm.c0 = [mol/m^3]
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double initial[9] = {54, 1e-7, 1e-7, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0.0};
    equi_t *e;
    err = equi_create(&e); if(err) return err;
    err = equi_parse(e, &p); if(err) return err;
    PetscBool set;
    err = PetscOptionsGetReal(NULL, NULL, "-acid_conc", &(e->component[0]->c), &set);
    if(err||!set) {
        printf("Acid concentration not set\n");
        return -1;
    }
    err = PetscOptionsGetReal(NULL, NULL, "-naoh_conc", &(e->component[1]->c), &set);
    if(err||!set) {
        printf("NaOH concentration not set\n");
        return -1;
    }
    PetscPrintf(PETSC_COMM_WORLD, "[%g|%g]", e->component[0]->c, e->component[1]->c);
    equi_get_initial(e, initial);
    equi_destroy(e);
    mesh_set_state(&mesh, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);
    copy_solution_from_mesh(&app);
  
    /* Run experiments */
    char name[256];
    err = PetscOptionsGetString(NULL, NULL, "-name", name, 256, &set);
    if(err||!set) {
        printf("Name of application not set. use -name\n");
        iface->finish(&app);
        return err;
    }
    err = app_set_name(&app, name); if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0.00, 25.002, 1); if(err) return err;

    //err = time_ivc(&app, initial, n); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
