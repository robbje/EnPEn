#include "mesh.h"

void mesh_print_volumes(mesh_t *mesh)
{
    unsigned int i;
    for(i = 0; i < mesh->num_volumes; i++) {
        volume_print(mesh->volumes[i]);
    }
}

void volume_print(volume_t *vol)
{
    if(!vol) return;
    printf("[%i] type: %i pos: (%g, %g, %g) V: %g charge: %g alpha: %g\n",
            vol->index, vol->type, vol->center.x, vol->center.y, vol->center.z,
            vol->volume, vol->charge, vol->alpha);
    unsigned int f;
    for(f = 0; f < vol->num_faces; f++)
        printf("\tDistance to neighbour %i: %g\n", f, vol->distance[f]);
}
