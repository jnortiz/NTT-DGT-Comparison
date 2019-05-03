#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "gauss.h"

unsigned less_than(const uint64_t x, const uint64_t y);
uint64_t select_uint64_t(uint64_t x, uint64_t y, uint64_t bit);
uint64_t mod_uint64_t(const uint64_t a);
void dgt(gauss_t* _x, gauss_t* input_signal);