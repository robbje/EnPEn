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
    param_set(&p, "sysdata/memristor.def", norm);

    meshdef_t def = TRILAYER;
    mesh_set(&mesh, &def);

    /* Set up a membrane with same charge as our P- concentration in equilibrium */
    zone_t cem = TRILAYER_MEMB;
    cem.alpha = 0.1;
    cem.end_charge = -0.1; // aem.charge*norm.c0 = [mol/m^3]
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    /* Initial concentrations in electrolyte solution, 0 PH/P- */
    //double initial[8] = {54.0, 1e-7, 1e-7, 0.0, 0.0, 1e-3, 1e-3, 0.0};

    /* Initial concentrations in polyelectrolyte */
    double initial_polyelectrolyte[8] = {54.0, 2e-7, 1e-7, 1e-7, 1e-7, 1e-3, 1e-3, 0.0};
    double *initial = &initial_polyelectrolyte[0];
    equi_t *e;
    err = equi_create(&e); if(err) return err;
    err = equi_parse(e, &p); if(err) return err;
    e->component[0]->c = 1e-4;
    equi_get_initial(e, initial_polyelectrolyte);
    equi_destroy(e);

    double charge = 0;
    for(int i = 0; i < 7; i++) {
        printf("%g\n", initial_polyelectrolyte[i]);
        charge += p.z[i]*initial_polyelectrolyte[i];
    }
    printf("pH = %g\n", -log(initial_polyelectrolyte[1])/log(10));
    printf("charge = %g\n", charge);
    printf("[P-] = %g\n", initial_polyelectrolyte[3]);
    PetscReal lambda = sqrt(p.epsilon*p.R*p.T/(2.0*initial[0]*p.c0*pow(p.F, 2.0)));
    PetscPrintf(PETSC_COMM_WORLD, "Debye length: %gm\n", lambda);


    mesh_set_state(&mesh, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, &initial[0], n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, &initial[0], n);
    copy_solution_from_mesh(&app);

    //cem.end_charge = -1.0*initial_polyelectrolyte[3];

    /* Solve initial solution */
    //err = sweep_zones(&app); if(err) return err;

    /* Perform ye olde switcharoo: Switch "artificial membrane" with "polyelectrolyte"
     * Make zone charge-ineffective */
    //mesh_apply_zones(&mesh);

    /* Replace zeroed charge zone with our equilibrium PE */
    //mesh_set_state_range(&mesh, 100.0, 200.0, initial_polyelectrolyte, n);

    /* Run experiments */
    err = app_set_name(&app, "test"); if(err) return err;
    err = sweep_boundary_chargeall(&app, FACETYPE_BOUNDARY_RIGHT, initial,
            &initial[p.num_species], 0, 1.5, 1.0); if(err) return err;

    //err = time_ivc(&app, initial_polyelectrolyte, n); if(err) return err;
    /* Finish and exit */
    return iface->finish(&app);
}
