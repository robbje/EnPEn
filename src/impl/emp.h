#ifndef EMP_H
#define EMP_H
#include <petscsnes.h>
#include "../app.h"
#include "../iface_petsc.h"

/* Calculate a field (gradient) of a variable over a volume */
PetscScalar field(volume_t *v, unsigned int index);
/* Update function to update the EMP in the system */
PetscErrorCode update_EMP(SNES snes, PetscInt step);
#endif
