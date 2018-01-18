#ifndef INTERFACE_C
#define INTERFACE_C
#include "interface.h"

io_t IO_DEFAULT = {
    .fopen = libc_fopen,
    .fclose = fclose,
    .fprintf = fprintf,
    .getline = libc_getline
};

int libc_fopen(const char* filename, const char *mode, FILE **fp) {
    *fp = fopen(filename, mode);
    if(*fp == NULL)  {
        perror("fopen");
        return -1;
    }
    return 0;
}

int libc_getline(FILE *fp, char **dst, size_t *len) {
    int err = getline(dst, len, fp);
    if(err >= 0) return 0;
    return err;
}

#endif
