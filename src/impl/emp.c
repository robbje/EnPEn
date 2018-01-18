#include "emp.h"

/* Calculate a field (gradient) of a variable over a volume */
PetscScalar field(volume_t *v, unsigned int index)
{
    double ksi = 0;
    ksi = v->distance[0] + v->distance[1];
    /* This is to get the sign of the field right */
    volume_t *lo,*hi;
    if(v->nvol[0]->center.x < v->nvol[1]->center.x) {
        lo = v->nvol[0];
        hi = v->nvol[1];
    } else {
        lo = v->nvol[1];
        hi = v->nvol[0];
    }
    return (hi->state[index] - lo->state[index])/ksi;
}

volume_t *boundary = NULL;
PetscErrorCode update_EMP(SNES snes, PetscInt step)
{
    PetscErrorCode ierr;
    app_t *app;
    ierr = SNESGetApplicationContext(snes, &app); CHKERRQ(ierr);
    ierr = copy_solution_to_mesh(app); CHKERRQ(ierr);
    mesh_t *mesh = app->mesh;
    parameter_t *p = app->p;
    /* This function integrates the rho_el*E from one left boundary to v->center.x */
    // First find the boundary once
    unsigned int i;
    if(!boundary) {
        if(mesh->boundaries[FACETYPE_BOUNDARY_LEFT])
            boundary = mesh->boundaries[FACETYPE_BOUNDARY_LEFT]->volume;
    }
    if(!boundary)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,
                "update_EMP: Could not find a boundary volume\n");
    /* Integrate from boundary to boundary, traversing the volumes */
    volume_t *act = boundary->nvol[0];
    volume_t *prev = boundary;
    double emp = 0;
    act->emp = emp;
    while(act->type != VOLTYPE_BULK && act->type != VOLTYPE_ELECTRODE) {
        
        /* Calculate charge density of the moving parts */
        double rho = 0.0;
        for(i = 0; i < p->num_species; i++) {
            rho += act->state[i] * p->z[i];
        }
        /* Calculate volume width */
        double ksi = 0.5 * (act->distance[0] + act->distance[1]);
        /* Calculate the electric field over the volume */
        double E = field(act, p->num_species);
        /* Calculate electromechanical pressure acting on this volume,
         * by taking the boundary as reference
         */
        emp -= rho*E*ksi;
        act->emp = emp;
        /* Traverse to the next volume */
        volume_t *tmp = act;
        act = act->nvol[0]==prev?act->nvol[1]:act->nvol[0];
        prev = tmp;
    }
    act->emp = emp;
    return 0;
}
