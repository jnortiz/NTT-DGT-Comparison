#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

typedef	int32_t poly[PARAM_N];

void poly_dgt(poly x_dgt, const poly x);
void poly_mul(poly output, const poly poly_a, const poly poly_b);

#endif
