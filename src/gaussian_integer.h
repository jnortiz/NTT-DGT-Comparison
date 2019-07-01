#ifndef GAUSSIAN_INTEGER_H
#define GAUSSIAN_INTEGER_H

#include <inttypes.h>

struct gauss {
    int32_t re;
    int32_t img;
};

typedef struct gauss gauss_t;

int32_t reduce(int64_t);
int32_t barr_reduce(int32_t);

void set_gauss(gauss_t*, const int32_t, const int32_t);
void add(gauss_t*, const gauss_t, const gauss_t);
void sub(gauss_t*, const gauss_t, const gauss_t);
void mul(gauss_t*, const gauss_t, const gauss_t);

int32_t add_int32_t(const int32_t, const int32_t);

#endif