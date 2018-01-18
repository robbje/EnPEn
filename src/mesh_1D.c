#include "mesh.h"

int mesh_parse_1D(mesh_t *mesh, io_t *io)
{
    int err;
    /* Read nodes */
    node_t *nodes;
    unsigned int num_nodes;
    err = mesh_read_nodes(mesh, io, &nodes, &num_nodes); if(err) return err;
    err = io_read_until(io, GMSH_MAGIC_ELEMS); if(err) return err;

    /* Allocate Memory to collect faces */
    face_t **facelist = malloc((num_nodes+1)*sizeof(face_t *));
    if(!facelist) {
        perror("mesh_parse_1D malloc");
        free(nodes);
        return -1;
    }
    memset(facelist, 0, (num_nodes+1)*sizeof(face_t *));

    /* Read elements */
    unsigned int num_elems, i;
    err = io_read_uint(io, &num_elems);
    /* Allocate memory for elements */
    mesh->volumes = malloc(num_elems * sizeof(volume_t *));
    if(!mesh->volumes) {
        perror("mesh_parse_1D malloc");
        free(nodes);
        free(facelist);
        return -1;
    }
    for(i = 0; i < num_elems; i++) {
        err = mesh_parse_elem_1D(mesh, io, nodes, facelist);
        if(err) {
            free(nodes);
            free(facelist);
            free(mesh->volumes);
            return err;
        }
    }
    free(nodes);
    
    /* Post processing */
    // Step 0: Realloc to the actual num_volumes
    mesh->volumes = realloc(mesh->volumes,mesh->num_volumes*sizeof(volume_t *));
    // Step 1: Collect all faces in facelist
    unsigned int num_faces = 0;
    for(i = 0; i < num_nodes+1; i++) {
        if(!facelist[i]) continue;
        num_faces++;
    }
    mesh->faces = malloc(num_faces * sizeof(face_t *));
    if(!mesh->faces) {
        perror("mesh_parse_1D malloc");
        free(facelist);
        return -1;
    }
    for(i = 0; i < num_nodes+1; i++) {
        if(!facelist[i]) continue;
        mesh->faces[mesh->num_faces++] = facelist[i];
    }
    free(facelist);
    // Step 2: Create the boundaries
    err = mesh_create_boundaries(mesh); if(err) return err;
    // Step 3: Final meshing and calculated distances
    return mesh_connect_volumes(mesh);
}

int mesh_parse_elem_1D(mesh_t *mesh, io_t *io, node_t *nodes, face_t **facelist)
{
    int err;
    size_t sze = IO_MAX_LINE_LEN;
    char *line = malloc(sze);
    if(!line) {
        perror("mesh_parse_elem_1D");
        return -1;
    }
    memset(line, 0, sze);
    err = io->getline(io->fh, &line, &sze); if(err < 0) goto elem_1D_err;

    /* Parse element in line */
    unsigned int index, type, num_tags, tag, i;
    /* Node indices */
    unsigned int a, b, tmp;
    char *t = strtok(line, " ");
    err = io_token_uint(&t, &index); if(err) goto elem_1D_err;
    err = io_token_uint(&t, &type); if(err) goto elem_1D_err;
    err = io_token_uint(&t, &num_tags); if(err) goto elem_1D_err;
    /* Consume tags */
    for(i = 0; i < num_tags; i++) {
        err = io_token_uint(&t, &tag); if(err) goto elem_1D_err;
    }
    switch(type) {
        case GMSH_LINE:
            /* Read two nodes making a line */
            err = io_token_uint(&t, &a); if(err) goto elem_1D_err;
            err = io_token_uint(&t, &b); if(err) goto elem_1D_err;
            /* Craft volume from line */
            if(a>b) {
                tmp = b; b = a; a = tmp;
            }
            err = mesh_add_vol_1D(mesh, nodes, facelist, a-1, b-1);
            if(err) goto elem_1D_err;
            break;
        case GMSH_NODE:
            /* Ignore nodes */
            free(line); return 0;
            break;
        default:
            fprintf(stderr, "Element type not supported for 1D mesh: %i\n",
                    type);
            goto elem_1D_err;
            break;
    }
    free(line);
    return 0;
elem_1D_err:
    free(line);
    return err;
}

int mesh_add_vol_1D(mesh_t *mesh, node_t *nodes, face_t **facelist,
        const unsigned int a, const unsigned int b)
{
    int err;
    double A = MESH_EXTENSION_LENGTH * MESH_EXTENSION_LENGTH;
    double V = A * node_distance(&nodes[a], &nodes[b]);

    /* Create the volume in memory */
    volume_t *vol;
    err = volume_new(&vol, 2, V, mesh->num_volumes); if(err) return err;
    mesh->volumes[mesh->num_volumes++] = vol;

    /* Calculate center of the volume */
    vol->center.x = (nodes[a].x+nodes[b].x)/2;
    vol->center.y = (nodes[a].y+nodes[b].y)/2;
    vol->center.z = (nodes[a].z+nodes[b].z)/2;

    /* Apply zone configuration to the volume */
    mesh_apply_zone_to_vol(mesh, vol);

    /* Identify or create faces for this volume */
    if(!facelist[a]) {
        err = face_new(&facelist[a], A); if(err) return err;
        facelist[a]->nvol[0] = vol;
        facelist[a]->center = nodes[a];
        facelist[a]->type = node_get_boundary(&nodes[a], &mesh->def->size);
    } else facelist[a]->nvol[1] = vol;
    if(!facelist[b]) {
        err = face_new(&facelist[b], A); if(err) return err;
        facelist[b]->nvol[0] = vol;
        facelist[b]->center = nodes[b];
        facelist[b]->type = node_get_boundary(&nodes[b], &mesh->def->size);
    } else facelist[b]->nvol[1] = vol;

    /* Connect the volume and face */
    vol->face[0] = facelist[a];
    vol->face[1] = facelist[b];
    return 0;
}
