#ifndef APP_H
#define APP_H

#include "interface.h"
#include "mesh.h"
#include "parameters.h"

typedef struct mesh_s mesh_t;
typedef struct parameter_s parameter_t;
typedef struct model_s model_t;

typedef struct app_s {
    parameter_t *p;
    mesh_t *mesh;
    void *solverdata;
    unsigned int nv;
    char *fname;
} app_t;

int app_init(app_t *app, io_t *io, unsigned int num_vars);
int app_dump_metadata(app_t *app, io_t *io);
int app_set_name(app_t *app, const char *fname);
void app_destroy(app_t *app);

#endif
