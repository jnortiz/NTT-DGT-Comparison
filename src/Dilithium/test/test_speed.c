#include <stdint.h>
#include "../sign.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../params.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 10000

uint64_t t[NTESTS];

static void bench_polymul() {
  poly a, b, c;
  poly_dgt(&a);
  poly_dgt(&b);
  poly_pointwise_montgomery(&c, &a, &b);
  poly_invdgt_tomont(&c);
}

int main(void)
{
  unsigned int i;
  size_t smlen;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sm[CRYPTO_BYTES + CRHBYTES];
  uint8_t seed[CRHBYTES];
  polyvecl mat[K];
  poly *a = &mat[0].vec[0];
  poly *b = &mat[0].vec[1];
  poly *c = &mat[0].vec[2];

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    expand_mat(mat, seed);
  }
  print_results("expand_mat:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    bench_polymul();
  }
  print_results("bench_polymul:", t, NTESTS);

/*
  unsigned int j;
  polyvecl *vl = &mat[0];
  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    for(j = 0; j < L; ++j)
      poly_uniform_eta(&vl->vec[j], seed, j);
  }
  print_results("polyvecl_uniform_eta:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    for(j = 0; j < L; ++j)
      poly_uniform_gamma1m1(&vl->vec[j], seed, j);
  }
  print_results("polyvecl_uniform_gamma1m1:", t, NTESTS);
*/

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    poly_dgt(a);
  }
  print_results("poly_dgt:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    poly_invdgt_tomont(a);
  }
  print_results("poly_invdgt_tomont:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    poly_pointwise_montgomery(c, a, b);
  }
  print_results("poly_pointwise_montgomery:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    challenge(c, seed, (polyveck *)mat);
  }
  print_results("challenge:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    crypto_sign_keypair(pk, sk);
  }
  print_results("Keypair:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    crypto_sign(sm, &smlen, sm, CRHBYTES, sk);
  }
  print_results("Sign:", t, NTESTS);

  for(i = 0; i < NTESTS; ++i) {
    t[i] = cpucycles();
    crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, CRHBYTES, pk);
  }
  print_results("Verify:", t, NTESTS);

  return 0;
}
