#include <stdio.h>
#include <stdlib.h>
#include "../random/random.h"
#include "cpucycles.h"
#include "../api.h"
#include "../poly.h"
#include "../pack.h"
#include "../sample.h"
#include "../params.h"
#include "../sha3/fips202.h"
  
#if (OS_TARGET == OS_LINUX)
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <unistd.h>
#endif

#define NRUNS 1000

static void poly_naivemul(poly c, const poly a, const poly b) {
  unsigned int i,j;
  uint32_t r[2*PARAM_N];

  for(i = 0; i < 2*PARAM_N; i++)
    r[i] = 0;

  for(i = 0; i < PARAM_N; i++)
    for(j = 0; j < PARAM_N; j++) {
      r[i+j] += ((uint64_t)a[i] * b[j]) % PARAM_Q;
      r[i+j] %= PARAM_Q;
    }

  for(i = PARAM_N; i < 2*PARAM_N-1; i++) {
    r[i-PARAM_N] = r[i-PARAM_N] + PARAM_Q - r[i];
    r[i-PARAM_N] %= PARAM_Q;
  }

  for(i = 0; i < PARAM_N; i++)
    c[i] = r[i];
}

static int cmp_llu(const void *a, const void*b)
{
  if (*(unsigned long long *)a < *(unsigned long long *)b) return -1;
  if (*(unsigned long long *)a > *(unsigned long long *)b) return 1;
  return 0;
}


static unsigned long long median(unsigned long long *l, size_t llen)
{
  qsort(l,llen,sizeof(unsigned long long),cmp_llu);

  if (llen%2) return l[llen/2];
  else return (l[llen/2-1]+l[llen/2])/2;
}


static unsigned long long average(unsigned long long *t, size_t tlen)
{
  unsigned long long acc=0;
  size_t i;
  for (i=0; i<tlen; i++)
    acc += t[i];
  return acc/(tlen);
}

static void print_results(const char *s, unsigned long long *t, size_t tlen)
{
  printf("%s", s);
  printf("\n");
  printf("median:  %llu ", median(t, tlen));  print_unit; printf("\n");
  printf("average: %llu ", average(t, tlen-1));  print_unit; printf("\n");
  printf("\n");
}

int main(void) {
  unsigned int i, j;
  unsigned long long cycles0[NRUNS], cycles1[NRUNS];
  uint16_t nonce = 0;
  unsigned char m1[PARAM_N], m2[PARAM_N];
  poly a, b;
  poly c1, c2;
  poly a_dgt, b_dgt;

  for(i = 0; i < NRUNS; ++i) {

    randombytes(m1, PARAM_N);
    randombytes(m2, PARAM_N);

    for(j = 0; j < PARAM_N; ++j) {
      a[j] = m1[j] - '0';
      b[j] = m2[j] - '0';
    }

    cycles0[i] = cpucycles();
    poly_naivemul(c1, a, b);
    cycles0[i] = cpucycles() - cycles0[i];

    cycles1[i] = cpucycles();
    poly_dgt(a_dgt, a);
    poly_dgt(b_dgt, b);
    poly_mul(c2, a_dgt, b_dgt);
    cycles1[i] = cpucycles() - cycles1[i];
  }

  print_results("naive: ", cycles0, NRUNS);
  print_results("dgt: ", cycles1, NRUNS);

  return 0;
}
