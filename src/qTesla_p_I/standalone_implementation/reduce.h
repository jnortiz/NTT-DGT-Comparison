#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

#define PARAM_Q 343576577
#define PARAM_QINV 2205847551
#define PARAM_BARR_MULT 3
#define PARAM_BARR_DIV 30

int32_t reduce(int64_t a);
int64_t barr_reduce64(int64_t a);

#endif