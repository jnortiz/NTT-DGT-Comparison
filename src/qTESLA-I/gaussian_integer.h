#ifndef GAUSSIAN_INTEGER_H
#define GAUSSIAN_INTEGER_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

/* Gaussian integer data type */
struct gauss {
    int32_t re;
    int32_t img;
};

typedef struct gauss gauss_t;

/* Auxiliary functions for gauss_t type */
void set_gauss(gauss_t*, const int32_t, const int32_t);
void print_gauss_t(const gauss_t);

int32_t mul_int32_t(const int32_t, const int32_t);
int32_t add_int32_t(const int32_t, const int32_t);

/* Operations over Gaussian integers */
void add(gauss_t*, const gauss_t, const gauss_t);
void sub(gauss_t*, const gauss_t, const gauss_t);
void mul(gauss_t*, const gauss_t, const gauss_t);

#endif