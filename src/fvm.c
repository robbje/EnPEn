#include "fvm.h"

/* Module stuff */
static char module_docstring[] = "FVM implementation by R. Femmer";

/* Methods def */
static PyMethodDef fvm_methods[] = {
    {"set_system", (PyCFunction)petsc_set_system, METH_VARARGS, "Set system and normalization constants"},
    {"set_mesh", (PyCFunction)petsc_set_mesh, METH_VARARGS, "Set mesh"},
    {"add_zone", (PyCFunction)petsc_add_zone, METH_VARARGS, "Add zone to the existing mesh"},
    {"set_composition", (PyCFunction)petsc_set_initial, METH_VARARGS, "Set the solution composition"},
    {"simulate", (PyCFunction)petsc_simulate_voltage, METH_VARARGS, "Simulate given potential drop"},
    {"load", (PyCFunction)petsc_load_app, METH_NOARGS, "Load mesh and parameters"},
    {NULL}  /* Sentinel */
};

/* Type definition */
static PyTypeObject fvm_PetscType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "fvm.petsc",               /*tp_name*/
    sizeof(fvm_PetscObject),   /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)petsc_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Petsc Interface",         /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    fvm_methods,               /* tp_methods */
    0,                         /* tp_members */ 
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)petsc_init_,     /* tp_init */
    0,                         /* tp_alloc */
    petsc_new,                 /* tp_new */
};

/* Methods implementation */
static PyObject *petsc_set_system(petsc *self, PyObject *args)
{
    char *fname;
    double c0, D0;
    if(!PyArg_ParseTuple(args, "sdd:", &fname, &c0, &D0))
        return 0;
    param_destroy(&(self->p));
    param_norm_t norm = {.c0=c0,.D0=D0};
    param_set(&(self->p), fname, norm);
    self->primed |= PRIMED_PARAMS;
    printf("[+] set_system: %s, %g, %g\n", fname, c0, D0);
    return Py_BuildValue("");
}

static PyObject *petsc_set_mesh(petsc *self, PyObject *args)
{
    meshdef_t *def;
    def = PyMem_Malloc(sizeof(meshdef_t));
    def->dim = 1;
    def->size.a.y = -1.0;
    def->size.a.z = -1.0;
    def->size.b.y = -1.0;
    def->size.b.z = -1.0;
    if(!PyArg_ParseTuple(args, "sddd:", &(def->filename), &(def->mesh_unit),
                                        &(def->size.a.x), &(def->size.b.x))) {
        return NULL;
    }
    mesh_destroy(&(self->mesh));
    mesh_set(&(self->mesh), def);
    mesh_add_zone(&(self->mesh), &(self->zone[0]));
    mesh_add_zone(&(self->mesh), &(self->zone[1]));
    self->primed |= PRIMED_MESH;
    printf("[+] set_mesh: %s\n", def->filename);
    return Py_BuildValue("");
}

static PyObject *petsc_add_zone(petsc *self, PyObject *args)
{
    unsigned int zi;
    double x0, x1, sigma;
    if(!PyArg_ParseTuple(args, "Iddd:", &zi, &x0, &x1, &sigma))
        return NULL;
    if(zi != 0 && zi != 1)
        return PyErr_Format(PyExc_ValueError, "Invalid zone index: %i\n", zi);
    self->zone[zi].end_charge = sigma;
    self->zone[zi].type = VOLTYPE_MEMBRANE;
    self->zone[zi].alpha = sigma==0?1:0.2;
    self->zone[zi].box.a.x = x0<x1?x0:x1;
    self->zone[zi].box.a.y = 0;
    self->zone[zi].box.a.z = 0;
    self->zone[zi].box.b.x = x0<x1?x1:x0;
    self->zone[zi].box.b.y = 0;
    self->zone[zi].box.b.z = 0;
    printf("[+] %s added (%g): %g -> %g\n", sigma>0?"AEM":"CEM", sigma, self->zone[zi].box.a.x, self->zone[zi].box.b.x);
    return Py_BuildValue("");
}

static PyObject *petsc_set_initial(petsc *self, PyObject *args)
{
    PyObject *list;
    if(!self->initialized) {
        PyErr_SetString(PyExc_StandardError, "Object not initialized: Call load() first!");
        return NULL;
    }
    if(!PyArg_ParseTuple(args, "O!", &PyList_Type, &list))
        return NULL;
    unsigned int len = PyList_Size(list);
    if(len != self->p.num_species) {
        return PyErr_Format(PyExc_ValueError, "List should be of size %u not %u", self->p.num_species, len);
    }
    self->initial = PyMem_Realloc(self->initial, (len+1) * sizeof(double));
    unsigned int i;
    for(i = 0; i < len; i++) {
        self->initial[i] = PyFloat_AsDouble(PyList_GetItem(list, i));
        if(self->initial[i] < 0.0) {
            PyErr_SetString(PyExc_ValueError, "Not a valid concentration");
            return NULL;
        }
    }
    self->initial[len] = 0.0;
    if(equi_solve(self->initial, &(self->p)) < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Could not solve equilibrium. Check your initial conditions");
        return NULL;
    }
    mesh_set_state(&(self->mesh), self->initial, len+1);
    boundary_set_state(&(self->mesh), FACETYPE_BOUNDARY_LEFT, self->initial, len+1);
    boundary_set_state(&(self->mesh), FACETYPE_BOUNDARY_RIGHT, self->initial, len+1);
    copy_solution_from_mesh(&(self->app));
    return Py_BuildValue("");
}

static PyObject *petsc_simulate_voltage(petsc *self, PyObject *args)
{
    if(!(self->initialized & (PRIMED_PARAMS|PRIMED_MESH))) {
        PyErr_SetString(PyExc_StandardError, "Not loaded");
        return NULL;
    }
    if(!self->initial) {
        PyErr_SetString(PyExc_StandardError, "No composition set: use set_composition first");
        return NULL;
    }
    double V;
    char *name;
    if(!PyArg_ParseTuple(args, "sd:", &name, &V)) return NULL;
    app_set_name(&(self->app), name);
    if(sweep_boundary_chargeall(&(self->app), FACETYPE_BOUNDARY_RIGHT, self->initial,
        &(self->initial[self->p.num_species]), V, V, 1.0)) {
        PyErr_SetString(PyExc_StandardError, "Something went wrong");
        return NULL;
    }
    return Py_BuildValue("");
}

static PyObject *petsc_load_app(petsc *self)
{
    if(self->initialized) {
        PyErr_SetString(PyExc_StandardError, "Already loaded!");
        return NULL;
    }
    if(!self->primed) {
        PyErr_SetString(PyExc_StandardError, "No system and mesh set: use set_system and set_mesh!");
        return NULL;
    }
    if(self->primed < 3 && self->primed & PRIMED_MESH) {
        PyErr_SetString(PyExc_StandardError, "No system set: use set_system!");
        return NULL;
    }
    if(self->primed < 3 && self->primed & PRIMED_PARAMS) {
        PyErr_SetString(PyExc_StandardError, "No mesh set: use set_mesh!");
        return NULL;
    }
    if(self->iface->load(&(self->app))) {
        PyErr_SetString(PyExc_StandardError, "Could not load application");
        return NULL;
    }
    model_t *model = (model_t *) self->app.solverdata;
    SNESSetTolerances(model->snes, 1e-8, 1e-30, 1e-15, 50, 1000);
    self->initialized = 1;
    self->primed = 0;
    printf("[+] App loaded\n");
    return Py_BuildValue("");
}

static PyObject *petsc_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    petsc *self;
    self = (petsc *)type->tp_alloc(type, 0);
    self->iface = &IFACE_PETSC;
    self->app.mesh = &(self->mesh);
    self->app.p = &(self->p);
    self->app.solverdata = (void *) &NPP_1D_REAC;
    self->initialized = 0;
    self->primed = 0;
    self->initial = NULL;
    self->zone[0].charge = 0;
    self->zone[1].charge = 0;
    self->zone[0].alpha = 1;
    self->zone[1].alpha = 1;
    /* TODO: parse the tuple *args into argc and argv to pass to iface->init() */
    if(self->iface->init(NULL, NULL)) {
        Py_DECREF(self);
        return NULL;
    }
    return (PyObject *) self;
}

static int petsc_init_(petsc *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

static petsc_dealloc(petsc *self)
{
    model_t *model;
    model = (model_t *) self->app.solverdata;
    model->finish(&(self->app));
    app_destroy(&(self->app));
    self->ob_type->tp_free((PyObject *) self);
}

/* Modules init */
PyMODINIT_FUNC initfvm(void)
{
    PyObject *m;

    if(PyType_Ready(&fvm_PetscType) < 0)
        return;
    m = Py_InitModule3("fvm", fvm_methods, module_docstring);
    Py_INCREF(&fvm_PetscType);
    PyModule_AddObject(m, "Petsc", (PyObject *) &fvm_PetscType);
    return;
}
