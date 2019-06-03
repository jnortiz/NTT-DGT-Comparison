#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

extern uint64_t p; 
extern int64_t n;
extern int64_t q;
extern unsigned __int128 max_128;

/* Gaussian integer data type */
struct gauss {
    uint64_t re;
    uint64_t img;
};

typedef struct gauss gauss_t;

/* Auxiliary functions for gauss_t type */
void set_gauss(gauss_t*, uint64_t, uint64_t);
void print_gauss_t(const gauss_t);
gauss_t conjugate(const gauss_t);

/* Constant-time operations */
unsigned less_than(const uint64_t, const uint64_t);
uint64_t select_uint64_t(uint64_t, uint64_t, uint64_t);

/* Reduction mod p = 0xFFFFFFFF00000001 */
uint64_t mod_uint64_t(const uint64_t);
uint64_t mod_uint128(const unsigned __int128);

/* Operations over Gaussian integers */
void add(gauss_t*, const gauss_t, const gauss_t);
void sub(gauss_t*, const gauss_t, const gauss_t);
void mul(gauss_t*, const gauss_t, const gauss_t);

/* Random 128-bit unsigned integer */
unsigned __int128 get_uint128_word();