#ifndef GAUSSIAN_INTEGER_H
#define GAUSSIAN_INTEGER_H

struct gauss {
    int32_t re;
    int32_t img;
};

typedef struct gauss gauss_t;

void add(gauss_t*, const gauss_t, const gauss_t);
void sub(gauss_t*, const gauss_t, const gauss_t);
void mul(gauss_t*, const gauss_t, const gauss_t);

#endif