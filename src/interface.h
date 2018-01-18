#ifndef INTERFACE_H
#define INTERFACE_H
#define IO_MAX_LINE_LEN 64
#include "app.h"
#include <stdio.h>
#include <stdarg.h>

typedef struct app_s app_t;
typedef struct interface_s {
    int (*init)(int *, char ***);
    int (*load)(app_t *);
    int (*finish)(app_t *);
    struct io_s *io;
} interface_t;

typedef struct io_s {
    int (*fopen)(const char *, const char *, FILE **);
    int (*fclose)(FILE *);
    int (*getline)(FILE *, char **, size_t *);
    int (*fprintf)(FILE *, const char *, ...);
    FILE *fh;
} io_t;

int libc_fopen(const char*, const char *, FILE **);
int libc_getline(FILE *, char **, size_t *);

#endif
