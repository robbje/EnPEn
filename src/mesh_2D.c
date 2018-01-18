#include "mesh.h"

int mesh_parse_2D(mesh_t *mesh, io_t *io)
{
    int err;
    node_t *nodes;
    unsigned int num_nodes;
    err = mesh_read_nodes(mesh, io, &nodes, &num_nodes); if(err) return err;
    err = io_read_until(io, GMSH_MAGIC_ELEMS); if(err) return err;

    return 0;
}

int mesh_parse_elem_2D(mesh_t *mesh, io_t *io, node_t *nodes, face_t **faces)
{

    return 0;
}

int mesh_add_vol_2D(mesh_t *mesh, node_t *nodes, face_t **faces)
{
    return 0;
}
