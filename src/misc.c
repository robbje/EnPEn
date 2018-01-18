#include "misc.h"

int io_read_until(io_t *io, const char *string)
{
    int err = -1;
    size_t sze = IO_MAX_LINE_LEN;
    char *line;
    line = malloc(sze);
    if(!line) {
        perror("io_read_until malloc");
        return err;
    }
    memset(line, 0, sze);
    while((err = io->getline(io->fh, &line, &sze)) >= 0) {
        if(strcmp(line, string) == 0) {
            free(line);
            return 0;
        }
    }
    free(line);
    fprintf(stderr, "Did not find string: %s\n", string);
    return err;
}

int io_token_uint(char **token, unsigned int *dst)
{
    if(!*token) {
        fprintf(stderr, "token is zero\n"); 
        return -1;
    }
    if(sscanf(*token, "%u", dst) != 1) {
        perror("sscanf");
        return -1;
    }
    *token = strtok(NULL, " ");
    return 0;
}

int io_read_uint(io_t *io, unsigned int *dst)
{
    int err = 0;
    size_t sze = IO_MAX_LINE_LEN;
    char *line;
    line = malloc(sze);
    if(!line) {
        perror("io_read_int malloc");
        return -1;
    }
    memset(line, 0, sze);
    err = io->getline(io->fh, &line, &sze);
    if(err < 0) {
        free(line);
        return err;
    }
    err = 0;
    if(sscanf(line, "%u\n", dst) != 1) err = -1;
    free(line);
    return err;
}
