#include <Python.h>
#include <structmember.h>
#include "app.h"
#include "iface_petsc.h"
#include "modules.h"
#include "parameters.h"
#include "mesh.h"
#include "simulation.h"
#include "equi.h"

#define PRIMED_MESH 1
#define PRIMED_PARAMS 2

typedef struct {
    PyObject_HEAD
    parameter_t p;
    mesh_t mesh;
    app_t app;
    meshdef_t meshdef;
    interface_t *iface;
    int initialized;
    int primed;
    double *initial;
    zone_t zone[2];
} fvm_PetscObject;
typedef fvm_PetscObject petsc;

/* Class methods definitions */
static PyObject *petsc_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int petsc_init_(petsc *self, PyObject *args, PyObject *kwds);
static petsc_dealloc(petsc *self);

static PyObject *petsc_set_system(petsc *self, PyObject *args);
static PyObject *petsc_set_mesh(petsc *self, PyObject *args);
static PyObject *petsc_add_zone(petsc *self, PyObject *args);
static PyObject *petsc_set_initial(petsc *self, PyObject *args);
static PyObject *petsc_simulate_voltage(petsc *self, PyObject *args);
static PyObject *petsc_load_app(petsc *self);
