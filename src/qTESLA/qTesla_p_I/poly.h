#ifndef POLY_H
#define POLY_H

#include "params.h"
#include "config.h"
#include <stdint.h>

typedef	int32_t poly[PARAM_N] __attribute__((aligned(32)));
typedef int64_t poly2x[PARAM_N] __attribute__((aligned(32)));
typedef	int32_t poly_k[PARAM_N*PARAM_K] __attribute__((aligned(32)));

int32_t reduce(int64_t a);
sdigit_t barr_reduce(sdigit_t a);
int64_t barr_reduce64(int64_t a);
void poly_pmul_asm(poly2x t, poly x, int32_t *nth);
void poly_dgt_asm(poly2x x_dgt, poly2x x, int32_t *gj);
void poly_idgt_asm(poly2x t_dgt, poly2x x_dgt, int32_t *invgj);
void poly_dgt(poly2x x_dgt, const poly x);
void poly_mul(poly result, const poly x, const poly2x y);
void poly_add(poly result, const poly x, const poly y);
void poly_add_correct(poly result, const poly x, const poly y);
void poly_sub(poly result, const poly x, const poly y);
void poly_sub_reduce(poly result, const poly x, const poly y);
void sparse_mul8(poly prod, const unsigned char *s, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H]);
void sparse_mul32(poly prod, const int32_t *pk, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H]);
void poly_uniform(poly_k a, const unsigned char *seed);

#endif
