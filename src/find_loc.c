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

    meshdef_t def = TRILAYER;
    mesh_set(&mesh, &def);

    /* Set up a membrane */
    zone_t cem = TRILAYER_MEMB;
    cem.alpha = 1;
    cem.end_charge = -100.0; // aem.charge*norm.c0 = [mol/m^3]
    cem.epsilon = 1;
    err = mesh_add_zone(&mesh, &cem); if(err) return err;

    //zone_t aem = SAID_AEM;
    //aem.alpha = 0.1;
    //aem.end_charge = 0.9;
    //aem.epsilon = 1;
    //err = mesh_add_zone(&mesh, &aem); if(err) return err;

    /* Load application */
    err = iface->load(&app); if(err) return err;

    int i;
    for(i = 1; i < argc; i++) {
        double loc = atof(argv[i]);
        int j;
        for(j = 0; j < mesh.num_volumes; j++) {
            if(mesh.volumes[j]->center.x > loc) {
                printf("%g ~= %g: %i\n",
                    loc,
                    mesh.volumes[j]->center.x,
                    j);
                    break;
            }
        }
    }


    return iface->finish(&app);
}
