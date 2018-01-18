/* MESHES */
/* filename, Lref, dimension, x0,y0,z0, x1,y1,y2 */
#define TRILAYER MESH_DEF("meshdata/trilayer-1D.msh", 1e-6, 1,\
        0, -1, -1, 300.0, -1, -1)
#define TRILAYER_MEMB {.box = MESH_BOX(100, 0, 0, 200.0, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define BPM MESH_DEF("meshdata/bpm-1D.msh", 1e-6, 1,\
        0, -1, -1, 400, -1, -1)
#define BPM_CEM {.box = MESH_BOX(100, 0, 0, 200, 0, 0),\
        .type=VOLTYPE_MEMBRANE}
#define BPM_AEM {.box = MESH_BOX(200.010, 0, 0, 300.000, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define BPMW MESH_DEF("meshdata/bpm-wide-1D.msh", 1e-6, 1,\
        0, -1, -1, 600, -1, -1)
#define BPMW_CEM {.box = MESH_BOX(200, 0, 0, 300, 0, 0),\
        .type=VOLTYPE_MEMBRANE}
#define BPMW_AEM {.box = MESH_BOX(300.002, 0, 0, 400.002, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define SAID MESH_DEF("meshdata/said-1D.msh", 1e-6, 1,\
        0, -1, -1, 300.0, -1, -1)
#define SAID_CEM {.box = MESH_BOX(100, 0, 0, 200, 0, 0),\
        .type=VOLTYPE_MEMBRANE}
#define SAID_AEM {.box = MESH_BOX(200.000, 0, 0, 200.1, 0, 0),\
        .type=VOLTYPE_MEMBRANE}
#define SAID_CEM2 {.box = MESH_BOX(200.400, 0, 0, 200.800, 0, 0),\
        .type=VOLTYPE_MEMBRANE}
#define SAID_AEM2 {.box = MESH_BOX(200.800, 0, 0, 201.200, 0, 0),\
        .type=VOLTYPE_MEMBRANE}
#define SAID_CEM3 {.box = MESH_BOX(201.200, 0, 0, 201.600, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define ELLIPSOCELL MESH_DEF("meshdata/ellipso-1D.msh", 1e-6, 1,\
        0, -1, -1, 1000, -1, -1)

#define ELLIPSOCELL_COATING {.box = MESH_BOX(0.0, 0, 0, 0.65, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define JOERI MESH_DEF("meshdata/joeri-1D.msh", 1e-6, 1,\
        0, -1, -1, 2500.0, -1, -1)
#define JOERI_CEM {.box = MESH_BOX(50.0, 0, 0, 250.0, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define JOERI_AEM {.box = MESH_BOX(2250.0, 0, 0, 2450.0, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define MAARTEN MESH_DEF("meshdata/maarten-1D.msh", 1e-5, 1,\
        0,-1,-1,210,-1,-1)
#define MAARTEN_MEMB {.box = MESH_BOX(100, 0, 0, 110, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define KHAIR MESH_DEF("meshdata/khair-1D.msh", 1e-6, 1,\
        0,-1,-1,100,-1,-1)

/* MESHES FOR TESTS - DO NOT MODIFY */
#define DBL MESH_DEF("meshdata/test-dbl-1D.msh", 1e-6, 1,\
        0, -1, -1, 1, -1, -1)

#define YARIV MESH_DEF("meshdata/test-ivc-1D.msh", 1e-6, 1,\
        0, -1, -1, 1.11, -1, -1)
#define YARIV_MEMB {.box = MESH_BOX(1, 0, 0, 1.1, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

#define JACTEST MESH_DEF("meshdata/jactest-1D.msh", 1e-6, 1,\
        0, -1, -1, 10, -1, -1)
#define JAC_MEMB {.box = MESH_BOX(4, 0, 0, 6, 0, 0),\
        .type=VOLTYPE_MEMBRANE}

/* MODELS */
#include "impl/npp_1D.h"
#include "impl/npp_1D_act.h"
#include "impl/npp_1D_reac.h"
#include "impl/npp_1D_reac_flow.h"
extern model_t NPP_1D;
extern model_t NPP_1D_REAC;
extern model_t NPP_1D_REAC_FLOW;
extern model_t NPP_1D_REAC_TS;
extern model_t NPP_1D_ACT;
/* SOLVERS */
extern interface_t IFACE_PETSC;

/* IO */
extern io_t IO_PETSC;
extern io_t IO_DEFAULT;
