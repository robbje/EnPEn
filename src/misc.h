#ifndef MISC_H
#define MISC_H
#include "interface.h"
typedef struct io_s io_t;
int io_read_until(io_t *io, const char *string);
int io_token_uint(char **token, unsigned int *dst);
int io_read_uint(io_t *io, unsigned int *dst);

#endif
