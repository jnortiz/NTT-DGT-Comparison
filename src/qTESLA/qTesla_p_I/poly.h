#ifndef POLY_H
#define POLY_H

#include "params.h"
#include "config.h"
#include <stdint.h>

typedef	int32_t poly[PARAM_N];
typedef	int32_t poly_k[PARAM_N*PARAM_K];

int32_t reduce(int64_t a);
sdigit_t barr_reduce(sdigit_t a);
int64_t barr_reduce64(int64_t a);
void dgt(poly x_dgt, const poly x);
void idgt(poly result, const poly x, const poly y);
void poly_dgt(poly x_dgt, const poly x);
void poly_pointwise(poly result, const poly x, const poly y);
void poly_invdgt(poly result, const poly x);
void poly_mul(poly result, const poly x, const poly y);
void poly_add(poly result, const poly x, const poly y);
void poly_add_correct(poly result, const poly x, const poly y);
void poly_sub(poly result, const poly x, const poly y);
void poly_sub_reduce(poly result, const poly x, const poly y);
void sparse_mul8(poly prod, const unsigned char *s, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H]);
void sparse_mul32(poly prod, const int32_t *pk, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H]);
void poly_uniform(poly_k a, const unsigned char *seed);

#endif
