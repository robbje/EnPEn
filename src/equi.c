#include "equi.h"

/* Selfmade solver */
int equi_create(equi_t **equi)
{
    *equi = malloc(sizeof(equi_t));
    if(!equi) {
        perror("malloc");
        return -1;
    }
    return 0;
}

int equi_parse(equi_t *equi, parameter_t *p)
{

    if(p->num_reactions < 2 || p->num_species < 4) {
        printf("Not enough reactions/species\n");
        return -1;
    }
    equi->num_components = 0;
    equi->tags = calloc(p->num_species, sizeof(unsigned int));
    if(!equi->tags) {
        perror("malloc");
        return -1;
    }

    /*** Identify tags */
    unsigned int i,j;
    for(j = 1; j < p->num_reactions; j++)
    {
        /* Find a tag for reaction j */
        unsigned int tag_j = 0;
        for(i = 3; i < p->num_species; i++) {
            int coeff = p->mu[i][j]+p->nu[i][j];
            if(coeff != 0 && equi->tags[i] != 0) {
                /* This reaction has tags[i] */
                tag_j = equi->tags[i];
            }
        }

        if(tag_j == 0) { /* Untagged reaction */
            tag_j = ++equi->num_components;
        }

        /* Tag all species to this component */
        for(i = 3; i < p->num_species; i++) {
            int coeff = p->mu[i][j]+p->nu[i][j];
            if(coeff != 0) {
                equi->tags[i] = tag_j;
            }
        }
    }

    equi->component = calloc(equi->num_components, sizeof(equi_component_t *));
    if(!equi->component) {
        perror("malloc");
        return -1;
    }

    /*** Add components with reactions */
    for(j = 1; j < p->num_reactions; j++) {
        /* Find a tag for reaction j */
        unsigned int tag_j = 0;
        for(i = 3; i < p->num_species; i++) {
            int coeff = p->mu[i][j]+p->nu[i][j];
            if(coeff != 0 && equi->tags[i] != 0) {
                /* This reaction has tags[i] */
                tag_j = equi->tags[i];
                if(p->z[i] == 0)
                    printf("Component %i: %s\n", tag_j-1, p->name[i]);
            }
        }

        /* j belongs to component tag_j */
        if(equi_add_reaction(equi, j, tag_j-1, p) < 0) return -1;
    }

    return 0;
}

int equi_add_reaction(equi_t *e, unsigned int reaction, unsigned int component,
        parameter_t *p)
{
    double K = p->k_f[reaction]/p->k_b[reaction];
    if(e->component[component] == NULL) {
        e->component[component] = malloc(sizeof(equi_component_t));
        if(!e->component[component]) {
            perror("malloc");
            return -1;
        }
        memset(e->component[component], 0, sizeof(equi_component_t));
    }

    equi_component_t *c = e->component[component];
    switch(equi_is_acid_reaction(p, reaction)) {
        case 1:
            c->num_acid_reacs++;
            c->Ka = realloc(c->Ka, c->num_acid_reacs*sizeof(double));
            if(c->Ka == NULL) {
                perror("realloc");
                return -1;
            }
            c->Ka[c->num_acid_reacs-1] = K; 
            break;
        case 0:
            c->num_base_reacs++;
            c->Kb = realloc(c->Kb, c->num_base_reacs*sizeof(double));
            if(c->Kb == NULL) {
                perror("realloc");
                return -1;
            }
            c->Kb[c->num_base_reacs-1] = 1e-14/K; 
            break;
        case 2:
            c->num_base_reacs++;
            c->Kb = realloc(c->Kb, c->num_base_reacs*sizeof(double));
            if(c->Kb == NULL) {
                perror("realloc");
                return -1;
            }
            c->Kb[c->num_base_reacs-1] = K; 
            break;
    }
    return 0;
}


int equi_is_acid_reaction(parameter_t *p, unsigned int j)
{
    /* Assuming Index of H+ is 1! */
    if(p->mu[1][j] != 0 || p->nu[1][j] != 0) {
        // We need to know if the educt is positively charged,
        // if so, this has the form of a base reaction
        int i;
        for(i = 0; i < p->num_species; i++) {
            if(p->mu[i][j] && p->z[i] > 0) {
                return 2;
            }
        }
        return 1;
    }
    return 0;
}

void equi_destroy(equi_t *equi)
{
    unsigned int i;
    free(equi->tags);
    for(i = 0; i < equi->num_components; i++) {
        free(equi->component[i]->Ka);
        free(equi->component[i]->Kb);
        free(equi->component[i]);
    }
    free(equi);
}

/* Actual solver */
double equi_charge_residual(equi_t *e, double x)
{
    double total = 0.0;
    int i, j, k;
    /* There is a problem here with acid reactions
     * which start positively charged somehow, i.e. AH+ -> H+ + A0
     * Maybe make charge of species availabe and use directly?
     */
    for(i = 0; i < e->num_components; i++) {
        double num = 0.0;
        double den = 0.0;
        double acid_num = 0.0;
        double acid_den = 0.0;
        for(j = 0; j < e->component[i]->num_acid_reacs; j++) {
            double acid_prod = 1.0;
            for(k = 0; k <= j; k++) {
                acid_prod *= e->component[i]->Ka[k]*x;
            }
            acid_num += acid_prod * (j+1); /* j for the charge of the species */
            acid_den += acid_prod;
        }
        double base_num = 0.0;
        double base_den = 0.0;
        for(j = 0; j < e->component[i]->num_base_reacs; j++) {
            double base_prod = 1.0;
            for(k = 0; k <= j; k++)
                base_prod /= (e->component[i]->Kb[k]*x);
            base_num -= base_prod * (j+1); /* j for the charge of the species */
            base_den += base_prod;
        }
        num = e->component[i]->c * x * (acid_num + base_num);
        den = 1 + acid_den + base_den;
        total += num/den;
    }
    return pow(10,-14)*x*x - 1 + total;
}

double equi_solve_ph(equi_t *e)
{
    double x = equi_bifurcate(&equi_charge_residual, 1, pow(10, 14), e, 0);
    return log(x)/log(10);
}

double equi_bifurcate(double (*func)(), double a, double b, equi_t *e, int n)
{
    double res_c = func((a+b)/2, e);
    double int_size = fabs(b-a);
    if(n > EQUI_MAX_ITER) {
        return (a+b)/2;
    }
    if(fabs(res_c) < EQUI_TOLF) {
        return (a+b)/2;
    }
    if(int_size < EQUI_TOLX) {
        return (a+b)/2;
    }
    if(res_c < 0)
        return equi_bifurcate(func, (a+b)/2, b, e, n+1);
    else if(res_c > 0)
        return equi_bifurcate(func, a, (a+b)/2, e, n+1);
    return (a+b)/2;

}

void equi_get_initial(equi_t *e, double *c)
{
    double pH = equi_solve_ph(e);
    double x = pow(10, pH);
    c[0] = 1000/18.5;
    c[1] = pow(10, -pH); 
    c[2] = pow(10,-14)/c[1];
    int i, j, k, index;
    index = 3;
    for(i = 0; i < e->num_components; i++) {
        double den1 = 1.0;
        // Acids
        for(j = e->component[i]->num_acid_reacs-1; j >= 0; j--) {
            double acid_prod = 1.0; 
            for(k = 0; k <= j; k++)
                acid_prod *= e->component[i]->Ka[k]*x;
            c[index++] = acid_prod;
            den1 += acid_prod;
        }
        int neutral_index = index++;
        // Bases
        for(j = 0; j < e->component[i]->num_base_reacs; j++) {
            double base_prod = 1.0;
            for(k = 0; k <= j; k++)
                base_prod /= (e->component[i]->Kb[k]*x);
            c[index++] = base_prod;
            den1 += base_prod;
        }
        c[neutral_index] = e->component[i]->c / den1;
        for(j = 1; j <= e->component[i]->num_acid_reacs; j++)
            c[neutral_index-j] *= c[neutral_index];
        for(j = 1; j <= e->component[i]->num_base_reacs; j++)
            c[neutral_index+j] *= c[neutral_index];
    }
}

/* Petsc solving */
PetscErrorCode equi_solve(PetscScalar *initial, parameter_t *p)
{
    PetscErrorCode ierr;
    SNES snes;
    Mat J;
    Vec x,r;
    PetscScalar *s;
    SNESConvergedReason reason;
    unsigned int n = p->num_species;
    unsigned int i;
    PetscPrintf(PETSC_COMM_WORLD, "{");
    PetscScalar rho_entry = 0.0;
    for(i = 0; i < n; i++) {
        if(i == n-1) {
            PetscPrintf(PETSC_COMM_WORLD, "%g}\n", initial[i]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "%g, ", initial[i]);
        }
        rho_entry += p->z[i]*initial[i];
    }
    PetscPrintf(PETSC_COMM_WORLD, "Charge: %g\n", rho_entry);
    ierr = SNESCreate(PETSC_COMM_SELF, &snes); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF, n, &x); CHKERRQ(ierr);
    ierr = VecDuplicate(x, &r); CHKERRQ(ierr);
    ierr = MatCreateSeqDense(PETSC_COMM_SELF, n, n, NULL, &J); CHKERRQ(ierr);
    ierr = MatSetUp(J);
    ierr = SNESSetFunction(snes,r,equi_function,(void*)p); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,J,J,equi_jacobian,(void*)p); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, 1e-8, 1e-30, 1e-15, 50, 1000); CHKERRQ(ierr);
    ierr = VecGetArray(x, &s); CHKERRQ(ierr);
    ierr = PetscMemcpy(s, initial, n*sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr = VecRestoreArray(x, &s); CHKERRQ(ierr);
    ierr = SNESSolve(snes, NULL, x); CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snes, &reason); CHKERRQ(ierr);
    if(reason < 0) goto equi_out;
    ierr = VecGetArray(x, &s); CHKERRQ(ierr);
    ierr = PetscMemcpy(initial, s, n*sizeof(PetscScalar)); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "{");
    PetscScalar rho = 0.0, Ix = 0.0, lambda, sigma = 0.0;
    for(i = 0; i < n; i++) {
        if(i == n-1) {
            PetscPrintf(PETSC_COMM_WORLD, "%g}\n", initial[i]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "%g, ", initial[i]);
        }
        rho += p->z[i]*initial[i];
        Ix += 0.5*pow(p->z[i],2)*initial[i]; 
        sigma += pow(p->F, 2.0)/(p->R*p->T) * initial[i]*p->c0 * p->D[i]*p->D0
            * pow(p->z[i], 2.0);
    }
    lambda = sqrt(p->epsilon*p->R*p->T/pow(p->F, 2.0)/Ix);
    PetscPrintf(PETSC_COMM_WORLD, "Charge: %g\n", rho);
    PetscPrintf(PETSC_COMM_WORLD, "Debye: %gm\n", lambda);
    PetscPrintf(PETSC_COMM_WORLD, "Conductivity: %g S/m\n", sigma);
    PetscPrintf(PETSC_COMM_WORLD, "pH: %g\n", -log(initial[1])/log(10));

    ierr = VecRestoreArray(x, &s); CHKERRQ(ierr);
equi_out:
    PetscPrintf(PETSC_COMM_WORLD, "Reason: %s\n", SNESConvergedReasons[reason]);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&r); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);
    return reason;
}

PetscErrorCode equi_function(SNES snes, Vec x, Vec r, void *context)
{
    PetscErrorCode ierr;
    parameter_t *p = (parameter_t *) context;

    PetscScalar *c;
    PetscScalar *f;
    ierr = VecGetArray(x, &c); CHKERRQ(ierr);
    ierr = VecGetArray(r, &f); CHKERRQ(ierr);

    unsigned int i, j, k;
    for(i = 0; i < p->num_species; i++) {
        f[i] = 0.0;
        for(j = 0; j < p->num_reactions; j++) {
            PetscInt r_coeff = p->mu[i][j] + p->nu[i][j];
            PetscScalar p_f = 1.0, p_b = 1.0;
            for(k = 0; k < p->num_species; k++) {
                p_f *= pow(c[k], p->mu[k][j]);
                p_b *= pow(c[k], -p->nu[k][j]);
            }
            f[i] += r_coeff * p->k_f[j] * (p_f - p->k_b[j]/p->k_f[j] * p_b);
        }
    }

    ierr = VecRestoreArray(x, &c); CHKERRQ(ierr);
    ierr = VecRestoreArray(r, &f); CHKERRQ(ierr);

    return 0;
}

PetscErrorCode equi_jacobian(SNES snes, Vec x, Mat jac, Mat B,
        void *context)
{
    PetscErrorCode ierr;
    parameter_t *p = (parameter_t *) context;

    PetscScalar *c;
    ierr = VecGetArray(x, &c); CHKERRQ(ierr);
    unsigned int i, j, k, m;
    for(i = 0; i < p->num_species; i++) {
        for(k = 0; k < p->num_species; k++) {
            PetscScalar ds_i = 0.0;
            for(j = 0; j < p->num_reactions; j++) {
                PetscInt r_coeff = p->mu[i][j] + p->nu[i][j];
                PetscScalar p_f = 1.0, p_b = 1.0;
                for(m = 0; m < p->num_species; m++) {
                    if(k == m) continue;
                    p_f *= pow(c[m], p->mu[m][j]);
                    p_b *= pow(c[m], -p->nu[m][j]);
                }
                ds_i += r_coeff * p->k_f[j] * (
                        p->mu[k][j]*pow(c[k],p->mu[k][j]-1) * p_f -
                        p->k_b[j]/p->k_f[j] *
                        (-p->nu[k][j]) * pow(c[k], -p->nu[k][j] - 1) * p_b);
            }
            // Set i, k
            ierr = MatSetValue(jac,i,k,ds_i,INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = VecRestoreArray(x, &c); CHKERRQ(ierr);
    return 0;
}
