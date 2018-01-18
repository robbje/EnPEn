#include "parameters.h"

void param_set(parameter_t *p, char *filename, param_norm_t n)
{
    memset(p, 0, sizeof(parameter_t));
    p->filename = filename;
    p->D0 = n.D0;
    p->c0 = n.c0;
}

int param_load(parameter_t *p, io_t *io)
{
    int error;
    char *line;
    unsigned int lineno = 1;
    unsigned int i, j;
    int ne;

    /* Allocate some memory */
    size_t sze = IO_MAX_LINE_LEN;
    line = malloc(sze);
    if(!line) goto param_err;
    memset(line, 0, sze);

    /* Start reading the file */
    error = io->fopen(p->filename, "r", &(io->fh)); if(error != 0) goto param_err;
    error = io->getline(io->fh, &line, &sze); if(error != 0) goto param_err;
    ne = sscanf(line, "%i %i %i", &(p->num_species),
                                  &(p->num_reactions),
                                  &(p->num_groups));
    if(ne != 3)
        goto param_fmterr;
    lineno++;

    /* Allocate memory */
    size_t N = p->num_species;
    size_t G = p->num_groups;
    size_t R = p->num_reactions;
    p->z = malloc(N * sizeof(int));
    p->mw = malloc(N * sizeof(double));
    p->D = malloc(N * sizeof(double));
    p->hydration = malloc(N * sizeof(double));
    p->group_valence = malloc(G * sizeof(int));
    p->nu = malloc(N * sizeof(int *));
    p->mu = malloc(N * sizeof(int *));
    p->group_nu = malloc(N * sizeof(int *));
    p->name = malloc(N * sizeof(char *));
    p->k_f = malloc(R * sizeof(double));
    p->k_b = malloc(R * sizeof(double));
    for(i = 0; i < N; i++) {
        p->nu[i] = malloc(R * sizeof(int));
        p->mu[i] = malloc(R * sizeof(int));
        p->group_nu[i] = malloc(G * sizeof(int));
        p->name[i] = malloc(256);
    }

    /* Parse species */
    for(i = 0; i < N; i++, lineno++) {
        error = io->getline(io->fh, &line, &sze); if(error != 0) goto param_err;
        ne = sscanf(line, "%s %i %lf %lf %lf", (p->name[i]),
                                               &(p->z[i]),
                                               &(p->mw[i]),
                                               &(p->D[i]),
                                               &(p->hydration[i]));
        if(ne != 5) goto param_fmterr;
    }

    /* Parse reactions */
    for(i = 0; i < R; i++, lineno++) {
        error = io->getline(io->fh, &line, &sze); if(error != 0) goto param_err; 
        char *token = strtok(line, " ");
        for(j = 0; j < N; j++) {
            int coeff;
            ne = sscanf(token, "%i", &coeff);
            if(ne != 1) goto param_fmterr;
            p->mu[j][i]=coeff>0?coeff:0; // forward
            p->nu[j][i]=coeff<0?coeff:0; // backward
            token = strtok(NULL, " ");
        }
        ne = sscanf(token, "%lf", &(p->k_f[i]));
        if(ne != 1) goto param_fmterr;
        token = strtok(NULL, " ");
        ne = sscanf(token, "%lf", &(p->k_b[i]));
        if(ne != 1) goto param_fmterr;
    }

    /* Parse group valences */
    for(i = 0; i < G; i++, lineno++) {
        error = io->getline(io->fh, &line, &sze); if(error != 0) goto param_err; 
        ne = sscanf(line, "%i", &(p->group_valence[i]));
        if(ne != 1) goto param_fmterr;
    }

    /* Parse group definitions */
    for(i = 0; i < N; i++, lineno++) {
        error = io->getline(io->fh, &line, &sze); if(error != 0) goto param_err; 
        char *token = strtok(line, " ");
        for(j = 0; j < G; j++) {
            ne = sscanf(token, "%i", &p->group_nu[i][j]);
            if(ne != 1) goto param_fmterr;
            token = strtok(NULL, " ");
        }
    }

    /* Parse parameters and constants */
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->T)); lineno++;
    if(ne != 1) goto param_fmterr;
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->eta)); lineno++;
    if(ne != 1) goto param_fmterr;
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->rho)); lineno++;
    if(ne != 1) goto param_fmterr;
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->closest_approach)); lineno++;
    if(ne != 1) goto param_fmterr;
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->F)); lineno++;
    if(ne != 1) goto param_fmterr;
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->R)); lineno++;
    if(ne != 1) goto param_fmterr;
    error = io->getline(io->fh, &line, &sze); if(error) goto param_err; 
    ne = sscanf(line, "%lf", &(p->epsilon)); lineno++;

    free(line);
    return 0;
param_err:
    perror("param_new");
    if(line) free(line);
    param_destroy(p);
    return -1;
param_fmterr:
    fprintf(stderr, "sscanf: %s in line %i\n", strerror(errno), lineno);
    free(line);
    param_destroy(p);
    return -1;
}

void param_destroy(parameter_t *p)
{
    unsigned int i;
    size_t N = p->num_species;
    if(p->z) free(p->z);
    if(p->mw) free(p->mw);
    if(p->D) free(p->D);
    if(p->hydration) free(p->hydration);
    if(p->group_valence) free(p->group_valence);
    if(p->k_f) free(p->k_f);
    if(p->k_b) free(p->k_b);
    for(i = 0; i < N; i++) {
        if(p->nu&&p->nu[i]) free(p->nu[i]);
        if(p->mu&&p->mu[i]) free(p->mu[i]);
        if(p->group_nu&&p->group_nu[i]) free(p->group_nu[i]);
        if(p->name&&p->name[i]) free(p->name[i]);
    }
    if(p->nu) free(p->nu);
    if(p->mu) free(p->mu);
    if(p->group_nu) free(p->group_nu);
    if(p->name) free(p->name);
}

double param_get_debye_length(parameter_t *p, double *c)
{
    double Ix = 0.0;
    unsigned int i;
    for(i = 0; i < p->num_species; i++) {
        Ix += 0.5*pow(p->z[i],2)*c[i]; 
    }
    return sqrt(p->epsilon*p->R*p->T/pow(p->F, 2.0)/Ix);
}
