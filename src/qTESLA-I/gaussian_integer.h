#ifndef GAUSSIAN_INTEGER_H
#define GAUSSIAN_INTEGER_H

struct gauss {
    int32_t re;
    int32_t img;
};

typedef struct gauss gauss_t;

void set_gauss(gauss_t*, const int32_t, const int32_t);
void mul(gauss_t*, const gauss_t, const gauss_t);

#endif