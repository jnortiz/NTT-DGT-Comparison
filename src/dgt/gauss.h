#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

struct gauss {
    uint64_t re;
    uint64_t img;
};

typedef struct gauss gauss_t;

void Gauss(gauss_t *pointer, uint64_t re, uint64_t img);
void print_gauss_t(const gauss_t x);
gauss_t conjugate(const gauss_t x);
unsigned less_than(const uint64_t x, const uint64_t y);
uint64_t select_uint64_t(uint64_t x, uint64_t y, uint64_t bit);
uint64_t mod_uint64_t(const uint64_t a);
uint64_t mod_uint128(const unsigned __int128 a);
void add(gauss_t *z, const gauss_t x, const gauss_t y);
void sub(gauss_t *z, const gauss_t x, const gauss_t y);
void mul(gauss_t *z, const gauss_t x, const gauss_t y);
unsigned __int128 get_uint128_word();