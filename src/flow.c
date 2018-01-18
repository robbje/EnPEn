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
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC_FLOW};
    param_norm_t norm = {.c0=1000,.D0=1.0e-9}; /* Molar scale (c0=1000) */
    param_set(&p, "sysdata/h2o_ita.def", norm);
    PetscBool set;
    err = PetscOptionsGetReal(NULL, NULL, "-velocity", &(p.u), &set);
    if(err||!set) {
        return -1;
    }

    meshdef_t def = TRILAYER;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t aem = TRILAYER_MEMB;
    aem.alpha = 0.1;
    aem.end_charge = -0.1; // aem.charge*norm.c0 = [mol/m^3]
    err = PetscOptionsGetReal(NULL, NULL, "-charge", &(aem.end_charge), &set);
    aem.epsilon = 1.0;
    err = mesh_add_zone(&mesh, &aem); if(err) return err;

    /* Load application */
    err = iface->load(&app); if(err) return err;

    /* Set system */
    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double initial[9] = {54.0, 1e-7, 1e-7, 1e-4, 1e-4, 3e-4, 1e-5, 1e-5, 0.0};
    equi_t *e;
    err = equi_create(&e); if(err) return err;
    err = equi_parse(e, &p); if(err) return err;
    e->component[0]->c = 1e-3;
    e->component[1]->c = 1e-4;
    err = PetscOptionsGetReal(NULL, NULL, "-naoh_conc", &(e->component[1]->c), &set);
    if(err||!set) {
        return -1;
    }
    equi_get_initial(e, initial);
    equi_destroy(e);

    printf("Pe = %g\n", 1e-6*p.u/norm.D0);
    printf("pH = %g\n", -1.0*log(initial[1])/log(10));
    PetscPrintf(PETSC_COMM_WORLD, "ITA=%g NaOH=%g\n", e->component[0]->c, e->component[1]->c);
    if(err||!set) {
        return -1;
    }

    mesh_set_state(&mesh, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    boundary_set_type(&mesh, FACETYPE_BOUNDARY_RIGHT, VOLTYPE_DOWNSTREAM);
    copy_solution_from_mesh(&app);

    /* Run experiments */
    char name[256];
    err = PetscOptionsGetString(NULL, NULL, "-name", name, 256, &set);
    if(err||!set) {
        return -1;
    }
    err = app_set_name(&app, name); if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 0, 1); if(err) return err;

    //PetscPrintf(PETSC_COMM_WORLD, "Ramping flow:\n");
    //err = sweep_flow(&app, &u, 0.0, 0.01, 0.0001); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
