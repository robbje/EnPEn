#include "parameters.h"
#include "equi.h"
#include "modules.h"
#include "iface_petsc.h"
#include "app.h"

PetscScalar u = 0;

void produce_line(double *c, unsigned int n)
{
    //int i;
    //for(i = 0; i < n; i++) {
    //    printf("%g, ", c[i]);
    //}
    printf("%g, %g\n", -c[3],  -log(c[1])/log(10));
}

int main(int argc, char **argv)
{
    int err;
    parameter_t p;
    mesh_t mesh;
    interface_t *iface = &IFACE_PETSC;
    err = iface->init(&argc, &argv); if(err) return err;
    app_t app = {.mesh=&mesh, .p=&p, .solverdata = (void *) &NPP_1D_REAC};
    meshdef_t def = JACTEST;
    mesh_set(&mesh, &def);

    param_norm_t norm = {.c0=1000,.D0=1.0e-9};
    param_set(&p, "sysdata/h2o_weakacid.def", norm);
    err = iface->load(&app);
    //int i;
    //for(i = 0; i < app.nv-1; i++) printf("[%s], ", p.name[i]); printf("pH\n");
    equi_t *e;
    equi_create(&e);
    equi_parse(e, &p);
    double *initial = malloc((app.nv-1)*sizeof(double));

    e->component[0]->c = 0.01;
    e->component[1]->c = 0.00001;
    equi_get_initial(e, initial);

    while(-log(initial[1])/log(10) <= 12.0) {
        equi_get_initial(e, initial);
        printf("%g, ", e->component[1]->c);produce_line(initial, app.nv-1);
        e->component[1]->c += 0.00001;
    }

    free(initial);

    return iface->finish(&app);
}
