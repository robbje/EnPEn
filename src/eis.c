#include "app.h"
#include "iface_petsc.h"
#include "modules.h"
#include "parameters.h"
#include "mesh.h"
#include "simulation.h"
#include "impedance.h"
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
    //param_set(&p, "sysdata/nacl.def", norm);
    //param_set(&p, "sysdata/h2o_weaksalt.def", norm);
    //param_set(&p, "sysdata/said.def", norm);
    //param_set(&p, "sysdata/cacl.def", norm);
    param_set(&p, "sysdata/cesar.def", norm);

    meshdef_t def = SAID;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t cem = SAID_CEM;
    cem.alpha = 0.1;
    cem.end_charge = -1.0; // aem.charge*norm.c0 = [mol/m^3]
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

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

    zone_t aem = SAID_AEM;
    aem.alpha = 0.1;
    //aem.end_charge = 0.1;
    aem.epsilon = 1;
    PetscBool set;
    err = PetscOptionsGetReal(NULL, NULL, "-charge", &(aem.end_charge), &set);
    if(err||!set) {
        return -1;
    }
    err = mesh_add_zone(&mesh, &aem); if(err) return err;
    
    /* Load application */
    err = iface->load(&app); if(err) return err;

    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    //double initial[7] = {54, 1e-7, 1e-7, 1e-3, 1e-3, 1e-3, 0.0};
    //double initial[3] = {1e-3, 1e-3, 0.0};
    double initial[5] = {1e-3, 1e-3, 1e-3, 5e-3, 0.0};
    //equi_t *e;
    //err = equi_create(&e); if(err) return err;
    //err = equi_parse(e, &p); if(err) return err;
    //e->component[0]->c = 1e-3;
    //equi_get_initial(e, initial);
    //equi_destroy(e);

    err = PetscOptionsGetReal(NULL, NULL, "-eis-bias", &(p.impedance.bias), &set);
    if(err||!set) {
        return -1;
    }
    err = PetscOptionsGetReal(NULL, NULL, "-eis-modulus", &(p.impedance.modulus), &set);
    if(err||!set) {
        return err;
    }
    err = PetscOptionsGetReal(NULL, NULL, "-eis-frequency", &(p.impedance.f), &set);
    if(err||!set) {
        return err;
    }

    err = impedance_setup(&app, &initial[0], n); if(err) return err;

    char name[256];
    err = PetscOptionsGetString(NULL, NULL, "-eis-name", name, 256, &set);
    if(err||!set) {
        printf("Name of application not set. use -eis-name\n");
        iface->finish(&app);
        return err;
    }
    err = app_set_name(&app, name); if(err) return err;
    err = impedance_run(&app); if(err) return err;

    /* Finish and exit */
    return iface->finish(&app);
}
