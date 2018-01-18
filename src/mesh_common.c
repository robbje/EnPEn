#include "mesh.h"

/* IO */
int mesh_read_nodes(mesh_t *mesh, io_t *io, node_t **nodes,
        unsigned int *num_nodes)
{
    int err;
    err = io_read_until(io, GMSH_MAGIC_NODES); if(err) return err;
    err = io_read_uint(io, num_nodes); if(err) return err;

    *nodes = malloc(*num_nodes*sizeof(node_t));
    if(!*nodes) {
        perror("mesh_read_nodes malloc");
        return -1;
    }
    unsigned int i;
    for(i = 0; i < *num_nodes; i++) {
        err = mesh_read_node(io, &(*nodes)[i]); if(err) return err;
    }
    return 0;
}

int mesh_read_node(io_t *io, node_t *node)
{
    int err;
    char *line;
    size_t sze = IO_MAX_LINE_LEN;
    line = malloc(sze);
    if(!line) {
        perror("mesh_read_node malloc");
        return -1;
    }
    memset(line, 0, sze);

    err = io->getline(io->fh, &line, &sze);
    if(err < 0) {
        free(line);
        return err;
    }
    int index;
    err = 0;
    if(sscanf(line, "%i %lf %lf %lf\n",
                &index,
                &(node->x),
                &(node->y),
                &(node->z)) != 4) err = -1;
    free(line);
    return err;
}

int mesh_dump_position(node_t *node, io_t *io)
{
    return io->fprintf(io->fh, "\"x\":%.8f,\"y\":%.8f,\"z\":%.8f,",
            node->x, node->y, node->z);
}

int mesh_dump_vector(double *vec, size_t len, io_t *io)
{
    int err;
    unsigned int i;
    err = io->fprintf(io->fh, "\"v\":["); if(err) return err;
    for(i = 0; i < len; i++) {
        err = io->fprintf(io->fh, "%.18f", vec[i]); if(err) return err;
        if(i < len-1) {
            err = io->fprintf(io->fh, ",");
            if(err) return err;
        }
    }
    err = io->fprintf(io->fh, "]"); if(err) return err;
    return 0;
}

int mesh_dump_volume(volume_t *v, io_t *io, unsigned int nv)
{
    int err;
    err = io->fprintf(io->fh, "{"); if(err) return err;
    err = mesh_dump_position(&(v->center), io); if(err) return err;
    err = io->fprintf(io->fh, "\"emp\":%.18f,", v->emp); if(err) return err;
    err = mesh_dump_vector(v->state, nv, io); if(err) return err;
    err = io->fprintf(io->fh, "}"); if(err) return err;
    return 0;
}

int mesh_dump_flux(face_t *f, io_t *io, unsigned int nv)
{
    int err;
    err = io->fprintf(io->fh, "{"); if(err) return err;
    err = mesh_dump_position(&(f->center), io); if(err) return err;
    err = mesh_dump_vector(f->flux, nv, io);
    if(err) return err;
    err = io->fprintf(io->fh, "}"); if(err) return err;
    return 0;
}

int mesh_find_volidx_for_node(mesh_t *mesh, node_t *n)
{
    double min_node_distance;
    min_node_distance = node_distance(&(mesh->def->size.a),
                                      &(mesh->def->size.b));
    int result = -1;
    unsigned int i;
    for(i = 0; i < mesh->num_volumes; i++) {
        double dist = node_distance(n, &(mesh->volumes[i]->center));
        if(dist < min_node_distance) {
            min_node_distance = dist;
            result = i;
        }
    }
    return result;
}

/* Geometry */
double node_distance(node_t *a, node_t *b)
{
    return sqrt(pow(a->x-b->x,2.0)+
                pow(a->y-b->y,2.0)+
                pow(a->z-b->z,2.0));
}

double triangle_area(node_t *a, node_t *b, node_t *c)
{
    double la = node_distance(a,b);
    double lb = node_distance(b,c);
    double lc = node_distance(c,a);
    double s = 0.5*(la+lb+lc);
    return sqrt(s*(s-la)*(s-lb)*(s-lc));
}

int node_in_box(node_t *a, box_t *s)
{
    if(a->x >= s->a.x && a->x <= s->b.x &&
       a->y >= s->a.y && a->y <= s->b.y &&
       a->z >= s->a.z && a->z <= s->b.z) return 1;
    return 0;
}

int node_get_boundary(node_t *a, box_t *s)
{
    if(a->x == s->a.x) return FACETYPE_BOUNDARY_LEFT;
    if(a->x == s->b.x) return FACETYPE_BOUNDARY_RIGHT;
    if(a->y == s->a.y) return FACETYPE_BOUNDARY_TOP;
    if(a->y == s->b.y) return FACETYPE_BOUNDARY_BOTTOM;
    return FACETYPE_INTERNAL;
}

/* Linked lists */
/* Zone list */
int zonelist_add(zonelist_t *zonelist, zone_t *zone)
{
    int err;
    zonelist_t *head = zonelist;
    while(head->next) head = head->next;
    err = zonelist_new(&head->next); if(err) return err;
    head->next->zone = zone;
    return 0;
}

int zonelist_new(zonelist_t **zonelist)
{
    *zonelist = malloc(sizeof(zonelist_t));
    if(!*zonelist) {
        perror("zonelist_new malloc");
        return -1;
    }
    (*zonelist)->next = NULL;
    (*zonelist)->zone = NULL;
    return 0;
}

void zonelist_destroy(zonelist_t *head)
{
    while(head) {
        zonelist_t *next = head->next;
        free(head);
        head = next;
    }
}


/* Volume list */
int volumelist_new(volumelist_t **list)
{
    *list = malloc(sizeof(volumelist_t));
    if(!*list) {
        perror("volumelist_new malloc");
        return -1;
    }
    (*list)->next = NULL;
    (*list)->volume = NULL;
    return 0;
}

int volumelist_add(volumelist_t **list, volume_t *v)
{
    int err;
    volumelist_t *head = *list;
    if(!*list) {
        err = volumelist_new(list); if(err) return err;
        (*list)->volume = v;
        return 0;
    }
    while(head->next) head = head->next;
    err = volumelist_new(&head->next); if(err) return err;
    head->next->volume = v;
    return 0; 
}

void volumelist_free(volumelist_t *head)
{
    while(head) {
        if(head->volume) volume_free(head->volume);
        volumelist_t *next = head->next;
        free(head);
        head = next;
    }
}

/* Face */
int face_new(face_t **face, const double area)
{
    if(area == 0.0) {
        fprintf(stderr, "new_face: can not create face with area 0\n");
        return -1;
    }
    *face = malloc(sizeof(face_t));
    if(!*face) {
        perror("new_face malloc");
        return -1;
    }
    (*face)->nvol[0] = NULL;
    (*face)->nvol[1] = NULL;
    (*face)->area = area;
    (*face)->type = FACETYPE_INTERNAL;
    (*face)->flux = NULL;
    return 0;
}

void face_free(face_t *face)
{
    free(face->flux);
    free(face);
}

/* Volume */
int volume_new(volume_t **vol,
               const unsigned int num_faces,
               const double volume,
               const unsigned int index)
{
    if(volume == 0.0) {
        fprintf(stderr, "new_volume: can not create volume with volume 0\n");
        return -1;
    }
    *vol = malloc(sizeof(volume_t));
    if(!*vol) {
        perror("new_volume malloc");
        return -1;
    }

    (*vol)->nvol = malloc(num_faces * sizeof(volume_t *));
    if(!(*vol)->nvol) {
        perror("new_volume malloc");
        free(*vol);
        return -1;
    }

    (*vol)->face = malloc(num_faces * sizeof(face_t *));
    if(!(*vol)->face) {
        perror("new_volume malloc");
        free((*vol)->nvol);
        free(*vol);
        return -1;
    }

    (*vol)->distance = malloc(num_faces * sizeof(double));
    if(!(*vol)->distance) {
        perror("new_volume malloc");
        free((*vol)->face);
        free((*vol)->nvol);
        free(*vol);
        return -1;
    }

    (*vol)->type = VOLTYPE_INTERNAL;
    (*vol)->num_faces = num_faces;
    (*vol)->volume = volume;
    (*vol)->index = index;
    (*vol)->state = NULL;
    (*vol)->charge = 0;
    (*vol)->alpha = 1.0;
    (*vol)->epsilon = 1.0;
    (*vol)->emp = 0.0;
    unsigned int i;
    for(i = 0; i < num_faces; i++) (*vol)->face[i] = NULL;
    return 0;
}

void volume_free(volume_t *vol)
{
    free(vol->distance);
    free(vol->state);
    free(vol);
}
/* Mesh */
void mesh_destroy(mesh_t *mesh)
{
    if(!mesh) return;
    unsigned int i;
    for(i = 0; i < mesh->num_volumes; i++) {
        if(mesh->volumes[i]) volume_free(mesh->volumes[i]);
    }
    for(i = 0; i < mesh->num_faces; i++) {
        if(mesh->faces[i]) face_free(mesh->faces[i]);
    }
    for(i = 0; i < 4; i++) {
        volumelist_free(mesh->boundaries[i]);
    }
    free(mesh->volumes);
    free(mesh->faces);
}

void mesh_apply_zone_to_vol(mesh_t *mesh, volume_t *volume)
{
    zonelist_t *head = mesh->zones;
    while(head) {
        if(node_in_box(&volume->center, &head->zone->box)) {
            volume->type = head->zone->type;
            volume->charge = head->zone->charge;
            volume->alpha = head->zone->alpha;
            volume->epsilon = head->zone->epsilon;
            return;
        }
        head = head->next;
    }
}

int mesh_connect_volumes(mesh_t *mesh)
{
    unsigned int i, j;
    for(i = 0; i < mesh->num_volumes; i++) {
        volume_t *v = mesh->volumes[i];
        for(j = 0; j < v->num_faces; j++) {
            face_t *f = v->face[j];
            if(!f->nvol[0]||!f->nvol[1]) {
                fprintf(stderr, "[!!] Face %i not properly connected\n", j);
                fprintf(stderr, "[!!] Check mesh size\n");
                return -1;
            }
            if(f->nvol[0] == v) {
                v->nvol[j] = f->nvol[1];
            } else {
                v->nvol[j] = f->nvol[0];
            }
            /* Calculate distance between volumes */
            if(v->nvol[j]->type == VOLTYPE_BULK ||
               v->nvol[j]->type == VOLTYPE_ELECTRODE) {
                v->distance[j] = node_distance(&v->center, &f->center);
            } else {
                v->distance[j] = node_distance(&v->center,
                    &(v->nvol[j]->center));
            }
        }
    }
    return 0;
}

int mesh_create_boundaries(mesh_t *mesh)
{
    int err;
    unsigned int i;
    for(i = 0; i < mesh->num_faces; i++) {
        face_t *f = mesh->faces[i];
        if(f->type == FACETYPE_INTERNAL) continue; // Ignore internal faces

        volume_t *v, *nv;
        nv = f->nvol[0];
        err = volume_new(&v, 1, nv->volume, 0); if(err) return err;
        v->type = VOLTYPE_BULK;
        f->nvol[1] = v;
        v->nvol[0] = nv;

        /* Dirty, dirty hack, which works for 1D only. */
        if(nv->index == 0) { 
            v->center.x = nv->center.x - 0.5;
        } else {
            v->center.x = nv->center.x + 0.5;
        }
        
        err = volumelist_add(&mesh->boundaries[f->type], v); if(err) return err;
    }
    return 0;
}

int mesh_init_states(mesh_t *mesh, unsigned int num_states)
{
    unsigned int i;
    for(i = 0; i < mesh->num_volumes; i++) {
        mesh->volumes[i]->state = malloc(num_states * sizeof(double));
        if(!mesh->volumes[i]->state) {
            perror("mesh_init_states malloc");
            return -1;
        }
        mesh->volumes[i]->bcond_type = malloc(num_states * sizeof(int));
        if(!mesh->volumes[i]->bcond_type) {
            perror("mesh_init_states malloc");
            return -1;
        }
    }
    for(i = 0; i < mesh->num_faces; i++) {
        mesh->faces[i]->flux = malloc(num_states * sizeof(double));
        if(!mesh->faces[i]->flux) {
            perror("mesh_init_states malloc");
            return -1;
        }
    }
    /* Allocate memory for boundary conditions */
    for(i = 0; i < 4; i++) {
        volumelist_t *head = mesh->boundaries[i];
        while(head) {
            volume_t *v = head->volume;
            v->state = malloc(num_states * sizeof(double));
            if(!v->state) {
                perror("mesh_init_states malloc");
                return -1;
            }
            head = head->next;
        }
    }
    return 0;
}
