/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"

extern int32_t gj[256], invgj[256];
extern int32_t nth[1024];
extern int32_t nthroots[1024];
extern int32_t invnthroots[1024];
extern int32_t invnthavx[1024];

void poly_uniform(poly_k a, const unsigned char *seed)
{ // Generation of polynomials "a_i"
  unsigned int pos=0, i=0, nbytes = (PARAM_Q_LOG+7)/8;
  unsigned int nblocks=PARAM_GEN_A;
  uint32_t val1, val2, val3, val4, mask = (uint32_t)(1<<PARAM_Q_LOG)-1;
  unsigned char buf[SHAKE128_RATE*PARAM_GEN_A];
  uint16_t dmsp=0;

  cshake128_simple(buf, SHAKE128_RATE*PARAM_GEN_A, dmsp++, seed, CRYPTO_RANDOMBYTES);

  while (i < PARAM_K*PARAM_N) {
    if (pos > SHAKE128_RATE*nblocks - 4*nbytes) {
      nblocks = 1;
      cshake128_simple(buf, SHAKE128_RATE*nblocks, dmsp++, seed, CRYPTO_RANDOMBYTES);
      pos = 0;
    }
    val1  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    val2  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    val3  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    val4  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    if (val1 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val1*PARAM_R2_INVN);
    if (val2 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val2*PARAM_R2_INVN);
    if (val3 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val3*PARAM_R2_INVN);
    if (val4 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val4*PARAM_R2_INVN);
  }
}


int32_t reduce(int64_t a)
{ // Montgomery reduction
  int64_t u;

  u = (a*PARAM_QINV) & 0xFFFFFFFF;
  u *= PARAM_Q;
  a += u;
  return (int32_t)(a>>32);
}


int64_t barr_reduce64(int64_t a)
{ // Barrett reduction
  int64_t u = (a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - u*PARAM_Q;
}


sdigit_t barr_reduce(sdigit_t a)
{ // Barrett reduction
  digit_t u = ((int64_t)a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - (digit_t)u*PARAM_Q;
}


void dgt(poly2x x_dgt, poly x, int32_t *gj)
{
  int m, window;
  int i, index, j, k;
  int32_t a, mul_re, mul_img;

  window = 512;
  for(m = 1; m < 512; m <<= 1)
  {
    index = 0;
    for(j = 0; j < m; j++)
    {
      a = gj[j];
      for(i = index; i < index + window; i = i+2)
      {
        //printf("%d %d %d %d\n", i, j, m, i+window);
        mul_re = reduce((int64_t)a * x[i+window]);
        mul_img = reduce((int64_t)a * x[i+window+1]);

        x[i+window] = x[i] - mul_re;
        x[i+window] += (x[i+window] >> (RADIX32-1)) & PARAM_Q;

        x[i+window+1] = x[i + 1] - mul_img;
        x[i+window+1] += (x[i+window+1] >> (RADIX32-1)) & PARAM_Q;

        x[i] = x[i] + mul_re - PARAM_Q;
        x[i] += (x[i] >> (RADIX32-1)) & PARAM_Q;

        x[i+1] = x[i + 1] + mul_img - PARAM_Q;
        x[i+1] += (x[i+1] >> (RADIX32-1)) & PARAM_Q;
      }
      index += (window << 1);
    }
    window >>= 1;
  }

  for (int i = 0; i < PARAM_N; i++) {
      if (x_dgt[i] != x[i]) {
          printf("DIFF %d %ld %d %lx %x\n", i, (int32_t)x_dgt[i], x[i], x_dgt[i], x[i]);
          exit(0);
      } else x[i] = (int32_t)x_dgt[i];
  }

  for (i = 0; i < PARAM_N; i++) x_dgt[i] = x[i];
}


static void poly_pointwise(poly result, const poly x, const poly y)
{ // Pointwise polynomial multiplication result = x.y

  for (int i=0; i<PARAM_N; i++)
    result[i] = reduce((int64_t)x[i]*y[i]);
}


void poly_dgt(poly2x x_dgt, const poly x)
{
  poly_pmul_asm(x_dgt, x, nth); //OK
  poly_dgt_asm(x_dgt, x_dgt, gj);
}


void poly_mul(poly result, const poly x, const poly2x y)
{ /* It is assumed that both signals are already in the DGT domain.
     The DGT counterpart of poly_b was computed in sign.c. */

  poly2x t_dgt;
  /* Calculating the point-wise multiplication of input signals */
  poly_pmul_asm2(t_dgt, x, y);
  /* Recovering the multiplication result in Z[x]/<x^n+1> */
  poly_idgt_asm(t_dgt, t_dgt, invgj);
  poly_pmul_asm3(result, t_dgt, invnthavx);
}


void poly_add(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y

    for(int i=0; i<PARAM_N; i++)
      result[i] = x[i] + y[i];
}


void poly_add_correct(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y with correction

    for (int i=0; i<PARAM_N; i++) {
      result[i] = x[i] + y[i];
      result[i] += (result[i] >> (RADIX32-1)) & PARAM_Q;    // If result[i] < 0 then add q
      result[i] -= PARAM_Q;
      result[i] += (result[i] >> (RADIX32-1)) & PARAM_Q;    // If result[i] >= q then subtract q
    }
}


void poly_sub(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y

    for (int i=0; i<PARAM_N; i++)
      result[i] = x[i] - y[i];
}


void poly_sub_reduce(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y

    for (int i=0; i<PARAM_N; i++)
      result[i] = (int32_t)barr_reduce(x[i] - y[i]);
}


/********************************************************************************************
* Name:        sparse_mul8
* Description: performs sparse polynomial multiplication
* Parameters:  inputs:
*              - const unsigned char* s: part of the secret key
*              - const uint32_t pos_list[PARAM_H]: list of indices of nonzero elements in c
*              - const int16_t sign_list[PARAM_H]: list of signs of nonzero elements in c
*              outputs:
*              - poly prod: product of 2 polynomials
*
* Note: pos_list[] and sign_list[] contain public information since c is public
*********************************************************************************************/
void sparse_mul8(poly prod, const unsigned char *s, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H])
{
  int i, j, pos;
  int8_t *t = (int8_t*)s;

  for (i=0; i<PARAM_N; i++)
    prod[i] = 0;

  for (i=0; i<PARAM_H; i++) {
    pos = pos_list[i];
    for (j=0; j<pos; j++) {
        prod[j] = prod[j] - sign_list[i]*t[j+PARAM_N-pos];
    }
    for (j=pos; j<PARAM_N; j++) {
        prod[j] = prod[j] + sign_list[i]*t[j-pos];
    }
  }
}


/********************************************************************************************
* Name:        sparse_mul32
* Description: performs sparse polynomial multiplication
* Parameters:  inputs:
*              - const int32_t* pk: part of the public key
*              - const uint32_t pos_list[PARAM_H]: list of indices of nonzero elements in c
*              - const int16_t sign_list[PARAM_H]: list of signs of nonzero elements in c
*              outputs:
*              - poly prod: product of 2 polynomials
*********************************************************************************************/
void sparse_mul32(poly prod, const int32_t *pk, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H])
{
  int i, j, pos;
  int64_t temp[PARAM_N] = {0};

  for (i=0; i<PARAM_H; i++) {
    pos = pos_list[i];
    for (j=0; j<pos; j++) {
        temp[j] = temp[j] - sign_list[i]*pk[j+PARAM_N-pos];
    }
    for (j=pos; j<PARAM_N; j++) {
        temp[j] = temp[j] + sign_list[i]*pk[j-pos];
    }
  }
  for (i=0; i<PARAM_N; i++)
    prod[i] = (int32_t)barr_reduce64(temp[i]);
}
