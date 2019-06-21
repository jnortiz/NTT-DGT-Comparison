#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

extern uint64_t p; 
extern uint64_t gap;
extern int n;

/* Gaussian integer data type */
struct gauss {
    uint64_t re;
    uint64_t img;
};

typedef struct gauss gauss_t;

/* Auxiliary functions for gauss_t type */
void set_gauss(gauss_t*, const uint64_t, const uint64_t);
void print_gauss_t(const gauss_t);

/* Reduction mod p = 0xFFFFFFFF00000001 */
uint64_t mod_uint128(const unsigned __int128);
uint64_t mul_uint64_t(const uint64_t, const uint64_t);
uint64_t add_uint64_t(const uint64_t, const uint64_t);

/* Operations over Gaussian integers */
void add(gauss_t*, const gauss_t, const gauss_t);
void sub(gauss_t*, const gauss_t, const gauss_t);
void mul(gauss_t*, const gauss_t, const gauss_t);