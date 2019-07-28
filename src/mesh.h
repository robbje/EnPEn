#ifndef MESH_H
#define MESH_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "interface.h"
#include "gmsh.h"

#define MESH_EXTENSION_LENGTH 1
#define MESH_BOX(x0,y0,z0,x1,y1,z1) {.a={.x=(x0),.y=(y0),.z=(z0)},\
        .b={.x=(x1),.y=(y1),.z=(z1)}}
#define MESH_DEF(fname,unit,dimension,x0,y0,z0,x1,y1,z1) {.filename=(fname),\
        .mesh_unit=(unit),\
        .dim=(dimension),\
        .size=MESH_BOX((x0),(y0),(z0),(x1),(y1),(z1))}

typedef struct io_s io_t;
typedef enum elem_type {ELEM_LINE, ELEM_TRIANGLE, ELEM_NODE} elem_type_t;
typedef enum face_type {FACETYPE_BOUNDARY_LEFT = 0,
                        FACETYPE_BOUNDARY_RIGHT = 1,
                        FACETYPE_BOUNDARY_TOP = 2,
                        FACETYPE_BOUNDARY_BOTTOM = 3,
                        FACETYPE_INTERNAL} face_type_t;
typedef enum vol_type {VOLTYPE_INTERNAL,
                       VOLTYPE_MEMBRANE,
                       VOLTYPE_BULK,
                       VOLTYPE_DOWNSTREAM,
                       VOLTYPE_ELECTRODE} vol_type_t;
typedef struct node_s {
    double x;
    double y;
    double z;
} node_t;

typedef struct box_s {
    struct node_s a;
    struct node_s b;
} box_t;

typedef struct face_s {
    face_type_t type;
    struct node_s center;
    double area;                /* The dimensionless area of this face */
    double *flux;
    struct volume_s *nvol[2];   /* Associated volumes */
} face_t;

typedef struct volume_s {
    unsigned int num_faces;
    unsigned int index;
    vol_type_t type;
    double *state;          /* State of this volume element */
    int *bcond_type;        /* Type of boundary condition per state */
    struct node_s center;   /* Geometrical center of volume element */
    double volume;          /* The dimensionless volume of this element */
    double charge;          /* The dimensionless background charge */
    double alpha;           /* The membrane diffusion factor */
    double emp;             /* Electromechanical pressure */
    double hsp;             /* Hydrostatic pressure */
    double epsilon;         /* Free volume ration */
    struct face_s **face;   /* Associated faces */
    struct volume_s **nvol; /* Neighbouring volumes */
    double *distance;       /* Neighbour-associated distance */
} volume_t;

typedef struct volumelist_s {
    struct volume_s *volume;
    struct volumelist_s *next;
} volumelist_t;

typedef struct zone_s {
    struct box_s box;
    int type;
    double alpha;
    double charge; 
    double end_charge;
    double old_charge;
    double epsilon;
} zone_t;

typedef struct zonelist_s {
    struct zone_s *zone;
    struct zonelist_s *next;
} zonelist_t;

typedef struct meshdef_s {
    char *filename;
    double mesh_unit;
    struct box_s size;
    unsigned int dim;
} meshdef_t;

typedef struct mesh_s {
    struct meshdef_s *def;
    unsigned int num_volumes;
    struct volume_s **volumes;
    unsigned int num_faces;
    struct face_s **faces;
    struct volumelist_s *boundaries[4];
    struct zonelist_s *zones;
} mesh_t;

/**** mesh.c */
/* Mesh User Interface */
void mesh_set(mesh_t *mesh, meshdef_t *def);
int mesh_parse(mesh_t *mesh, io_t *io);
int mesh_add_zone(mesh_t *mesh, zone_t *zone);
void mesh_set_state(mesh_t *mesh, double *state, unsigned int len);
void mesh_set_state_range(mesh_t *mesh, double a, double b, double *state, unsigned int len);
void boundary_set_state(mesh_t *mesh, int boundary, double *state,
        unsigned int len);
void boundary_set_type(mesh_t *mesh, int boundary, int type);
void boundary_modify_state(mesh_t *mesh, int boundary, double state, unsigned int index);
void mesh_apply_zones(mesh_t *mesh);

int mesh_open_dumpfile(io_t *io, const char *filename);
int mesh_dump(mesh_t *mesh, io_t *io, double indexval, unsigned int nv);
int mesh_close_dumpfile(io_t *io);

/**** mesh_1D.c */
int mesh_parse_1D(mesh_t *mesh, io_t *io);
int mesh_parse_elem_1D(mesh_t *mesh, io_t *io, node_t *nodes,
        face_t **facelist);
int mesh_add_vol_1D(mesh_t *mesh, node_t *nodes, face_t **facelist, 
        const unsigned int a, const unsigned int b);

/**** mesh_common.c */
/* IO */
int mesh_read_nodes(mesh_t *mesh, io_t *io, node_t **nodes, 
        unsigned int *num_nodes);
int mesh_read_node(io_t *io, node_t *node);
int mesh_dump_position(node_t *node, io_t *io);
int mesh_dump_vector(double *vec, size_t len, io_t *io);
int mesh_dump_volume(volume_t *v, io_t *io, unsigned int nv);
int mesh_dump_flux(face_t *f, io_t *io, unsigned int nv);
int mesh_find_volidx_for_node(mesh_t *mesh, node_t *n);
/* Geometry */
double node_distance(node_t *a, node_t *b);
double triangle_area(node_t *a, node_t *b, node_t *c);
int node_in_box(node_t *a, box_t *s);
int node_get_boundary(node_t *a, box_t *s);
/* Linked lists */
/* Zonelist */
int zonelist_add(zonelist_t *zonelist, zone_t *zone);
int zonelist_new(zonelist_t **zonelist);
void zonelist_destroy(zonelist_t *head);
/* Volumelist */
int volumelist_new(volumelist_t **list);
int volumelist_add(volumelist_t **list, volume_t *v);
void volumelist_free(volumelist_t *head);
/* Face */
int face_new(face_t **face, const double area);
void face_free(face_t *face);
/* Volume */
int volume_new(volume_t **vol,
        const unsigned int num_faces,
        const double volume,
        const unsigned int index);
void volume_free(volume_t *vol);
/* Mesh */
void mesh_destroy(mesh_t *mesh);
void mesh_apply_zone_to_vol(mesh_t *mesh, volume_t *volume);
int mesh_connect_volumes(mesh_t *mesh);
int mesh_create_boundaries(mesh_t *mesh);
int mesh_init_states(mesh_t *mesh, unsigned int num_states);
/**** mesh_debug.c */
void mesh_print_volumes(mesh_t *mesh);
void volume_print(volume_t *vol);
#endif
