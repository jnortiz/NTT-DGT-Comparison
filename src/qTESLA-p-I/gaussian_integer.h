#ifndef GAUSSIAN_INTEGER_H
#define GAUSSIAN_INTEGER_H

struct gauss {
    int64_t re;
    int64_t img;
};

typedef struct gauss gauss_t;

void set_gauss(gauss_t*, const int64_t, const int64_t);
void add(gauss_t*, const gauss_t, const gauss_t);
void sub(gauss_t*, const gauss_t, const gauss_t);
void mul(gauss_t*, const gauss_t, const gauss_t);

#endif