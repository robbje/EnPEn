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
    parameter_t p;
    mesh_t mesh;
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC_FLOW};
    param_norm_t norm = {.c0=1000,.D0=1.0e-9}; /* Molar scale (c0=1000) */
    param_set(&p, "sysdata/nacl.def", norm);
    meshdef_t def = JACTEST;
    mesh_set(&mesh, &def);
    //zone_t memb = JAC_MEMB;
    //memb.alpha = 0.1;
    //memb.charge = 1.0;
    //memb.end_charge = 1.0;
    //memb.epsilon = 0.1;
    //err = mesh_add_zone(&mesh, &memb); if(err) return err;
    err = iface->load(&app); if(err) return err;
    unsigned int n = p.num_species + NPP_1D.num_add_vars;
    double initial[3] = {1e-3, 1e-3, 0.0};
    p.u = 0.0;
    mesh_set_state(&mesh, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_LEFT, initial, n);
    boundary_set_state(&mesh, FACETYPE_BOUNDARY_RIGHT, initial, n);
    boundary_set_type(&mesh, FACETYPE_BOUNDARY_RIGHT, VOLTYPE_DOWNSTREAM);
    copy_solution_from_mesh(&app);
    model_t *model = (model_t *) app.solverdata;
    err = SNESSolve(model->snes, NULL, model->x); CHKERRQ(err);
    return iface->finish(&app);
}
