#include "mesh.h"

/* UI Calls */
void mesh_set(mesh_t *mesh, meshdef_t *def)
{
    memset(mesh, 0, sizeof(mesh_t));
    mesh->def = def;
    mesh->num_volumes = 0;
    mesh->num_faces = 0;
    mesh->zones = NULL;
}

int mesh_parse(mesh_t *mesh, io_t *io)
{
    int err;
    err = io->fopen(mesh->def->filename, "r", &(io->fh)); if(err) return err;
    switch(mesh->def->dim) {
        case 1: return mesh_parse_1D(mesh, io);
        case 2:
        case 3:
        default:
            fprintf(stderr, "Dimension unsupported: %i\n", mesh->def->dim);
            return -1;
    }
}

void mesh_set_state(mesh_t *mesh, double *state, unsigned int len)
{
    unsigned int vi;
    for(vi = 0; vi < mesh->num_volumes; vi++) {
        memcpy(mesh->volumes[vi]->state, state, len * sizeof(double));
    }
}

void mesh_set_state_range(mesh_t *mesh, double a, double b, double *state, unsigned int len)
{
    unsigned int vi;
    for(vi = 0; vi < mesh->num_volumes; vi++) {
        double x = mesh->volumes[vi]->center.x;
        if(a <= x && x <= b) {
            memcpy(mesh->volumes[vi]->state, state, len * sizeof(double));
        }
    }
}

void boundary_set_state(mesh_t *mesh, int boundary, double *state,
        unsigned int len)
{
    volumelist_t *head = mesh->boundaries[boundary];
    if(!head) fprintf(stderr, "Boundary %i does not exist\n", boundary);
    while(head) {
        memcpy(head->volume->state, state, len * sizeof(double));
        head = head->next;
    }
}

void boundary_set_type(mesh_t *mesh, int boundary, int type)
{
    volumelist_t *head = mesh->boundaries[boundary];
    if(!head) fprintf(stderr, "Boundary %i does not exist\n", boundary);
    while(head) {
        head->volume->type = type;
        head = head->next;
    }
}

void boundary_modify_state(mesh_t *mesh, int boundary, double state, unsigned int index)
{
    volumelist_t *head = mesh->boundaries[boundary];
    if(!head) fprintf(stderr, "Boundary %i does not exist\n", boundary);
    while(head) {
        head->volume->state[index] = state;
        head = head->next;
    }
}

int mesh_add_zone(mesh_t *mesh, zone_t *zone)
{
    int err;
    if(mesh->zones == NULL) {
        err = zonelist_new(&mesh->zones); if(err) return err;
        mesh->zones->zone = zone;
        return 0;
    }
    return zonelist_add(mesh->zones, zone);
}

void mesh_apply_zones(mesh_t *mesh)
{
    unsigned int i;
    for(i = 0; i < mesh->num_volumes; i++) {
        mesh_apply_zone_to_vol(mesh, mesh->volumes[i]);
    }
}

int mesh_open_dumpfile(io_t *io, const char *filename)
{
    int err;
    err = io->fopen(filename, "w", &io->fh); if(err) return err;
    err = io->fprintf(io->fh, "{"); if(err) return err;
    return 0;
}

int mesh_dump(mesh_t *mesh, io_t *io, double indexval, unsigned int nv)
{
    int err;
    unsigned int i;
    err = io->fprintf(io->fh, "\"%.3f\":{\"stt\":[", indexval);
    /* Dump all volumes, except BC */
    for(i = 0; i < mesh->num_volumes; i++) {
        err = mesh_dump_volume(mesh->volumes[i], io, nv);
        if(err) return err;
        if(i < mesh->num_volumes-1) {
            err = io->fprintf(io->fh, ",");
            if(err) return err;
        }
    }
    err = io->fprintf(io->fh, "],\n\"flx\":["); if(err) return err;
    /* Dump all fluxes */
    for(i = 0; i < mesh->num_faces; i++) {
        if(mesh->faces[i]->type == FACETYPE_INTERNAL) {
            err = mesh_dump_flux(mesh->faces[i], io, nv);
            if(err) return err;
            if(i < mesh->num_faces-1) {
                err = io->fprintf(io->fh, ",");
                if(err) return err;
            }
        }
    }
    err = io->fprintf(io->fh, "]},\n\n"); if(err) return err;
    return 0;
}

int mesh_close_dumpfile(io_t *io)
{
    int err;
    err = io->fprintf(io->fh, "\"empty\":[]}"); if(err) return err;
    err = io->fclose(io->fh); if(err) return err;
    return 0;
}
