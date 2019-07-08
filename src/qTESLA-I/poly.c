/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"
#include "params.h"
#include "params_dgt.h"

extern poly zeta;
extern poly zetainv;

int32_t reduce(int64_t a)
{ // Montgomery reduction
  int64_t u;

  u = (a*PARAM_QINV) & 0xFFFFFFFF;
  u *= PARAM_Q;
  a += u;
  return (int32_t)(a>>32);
}


void dgt(poly _x)
{    
    int i, index, j, l, m, stride, bound;
    int32_t sub_re, sub_img;

    m = PARAM_N;
    index = 0;
    bound = (PARAM_K2 >> 1);
    for(stride = 0; stride < PARAM_K2_LOG; stride++) {        
        m >>= 1;
        for(l = 0; l < bound; l++) {            
            j = _j_values[index];
            i = _i_values[index++];

            sub_re = _x[i]-_x[i+m];
            sub_img = _x[i+1]-_x[i+m+1];
            
            _x[i] += _x[i+m];
            _x[i+1] += _x[i+m+1];
            
            _x[i+m] = reduce((int64_t)__gj[j][stride]*sub_re);
            _x[i+m+1] = reduce((int64_t)__gj[j][stride]*sub_img);
        }
    }
}


void idgt(poly _output_signal)
{
    int i, index, j, l, m, stride, bound;
    int32_t mul_re, mul_img;

    index = 0;
    m = 2;
    bound = (PARAM_K2 >> 1);
    for(stride = 0; stride < PARAM_K2_LOG; stride++) {
        for(l = 0; l < bound; l++) {
            j = _idgt_j_values[index];
            i = _idgt_i_values[index++];
            
            mul_re = reduce((int64_t)_output_signal[i+m]*__invgj[j][stride]);
            mul_img = reduce((int64_t)_output_signal[i+m+1]*__invgj[j][stride]);
            
            _output_signal[i+m] =  _output_signal[i]-mul_re;
            _output_signal[i+m+1] = _output_signal[i+1]-mul_img;
            
            _output_signal[i] += mul_re;
            _output_signal[i+1] += mul_img;
        }
        m <<= 1;
    }
}

void poly_mul(poly _output, const poly _poly_a, const poly _poly_b)
{
    poly _folded_a, _folded_b;
    poly _mul, _output_gaussian;
    int64_t t1, t2, t3;
    int i, j;

    for(i = 0, j = 0; i < PARAM_N && j < PARAM_K2; i+=2, j++) {             
        /* Computing the folded signal. The same is done for both input signals */
        t1 = ((int64_t)_poly_a[j]*__nthroots[i]);
        t2 = ((int64_t)_poly_a[PARAM_K2+j]*__nthroots[i+1]);  
        t3 = ((int64_t)(_poly_a[j]+_poly_a[PARAM_K2+j])*(__nthroots[i]+__nthroots[i+1]));
        _folded_a[i] = reduce(t1-t2);
        _folded_a[i+1] = reduce(t3-t1-t2);           
        
        t1 = ((int64_t)_poly_b[j]*__nthroots[i]);
        t2 = ((int64_t)_poly_b[PARAM_K2+j]*__nthroots[i+1]);  
        t3 = ((int64_t)(_poly_b[j]+_poly_b[PARAM_K2+j])*(__nthroots[i]+__nthroots[i+1]));
        _folded_b[i] = reduce((int64_t)t1-t2);
        _folded_b[i+1] = reduce((int64_t)t3-t1-t2);
    } 

    /* Computing the DGT of signal b. The signal a is assumed to be in the DGT domain */
    dgt(_folded_b);

    /* Calculating the point-wise multiplication of input signals */
    for(i = 0; i < PARAM_N; i+=2) {
        t1 = ((int64_t)_folded_a[i]* _folded_b[i]);
        t2 = ((int64_t)_folded_a[i+1]*_folded_b[i+1]);  
        t3 = ((int64_t)(_folded_a[i]+_folded_a[i+1])*(_folded_b[i]+ _folded_b[i+1]));
        _mul[i] = reduce(t1-t2);
        _mul[i+1] = reduce(t3-t1-t2);
    }

    /* Recovering the multiplication result in Z[x]/<x^n+1> */
    idgt(_mul);

    /* Removing the twisting factors and writing the result from the Gaussian integer to the polynomial form */
    for(i = 0, j = 0; i < PARAM_N && j < PARAM_K2; i+=2, j++) {
        t1 = ((int64_t)_mul[i]*__invnthroots[i]);
        t2 = ((int64_t)_mul[i+1]*__invnthroots[i+1]);  
        t3 = ((int64_t)(_mul[i]+_mul[i+1])*(__invnthroots[i]+__invnthroots[i+1]));
        _output[j] = reduce(t1-t2);
        _output[j+PARAM_K2] = reduce(t3-t1-t2);   
    }
}

void poly_uniform(poly a, const unsigned char *seed)         
{ // Generation of polynomial "a"
  unsigned int pos=0, i=0, nbytes = (PARAM_Q_LOG+7)/8;
  unsigned int nblocks=PARAM_GEN_A;
  uint32_t val1, val2, val3, val4, mask = (1<<PARAM_Q_LOG)-1;
  unsigned char buf[SHAKE128_RATE*PARAM_GEN_A];
  uint16_t dmsp=0;

  cshake128_simple(buf, SHAKE128_RATE*PARAM_GEN_A, dmsp++, seed, CRYPTO_RANDOMBYTES);    
  
  while (i < PARAM_N) {  
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
    if (val1 < PARAM_Q && i < PARAM_N)
      a[i++] = reduce((int64_t)val1*PARAM_R2_INVN);
    if (val2 < PARAM_Q && i < PARAM_N)
      a[i++] = reduce((int64_t)val2*PARAM_R2_INVN);
    if (val3 < PARAM_Q && i < PARAM_N)
      a[i++] = reduce((int64_t)val3*PARAM_R2_INVN);
    if (val4 < PARAM_Q && i < PARAM_N)
      a[i++] = reduce((int64_t)val4*PARAM_R2_INVN);
  }
}


void ntt(poly a, const poly w)
{ // Forward NTT transform
  int NumoProblems = PARAM_N>>1, jTwiddle=0;

  for (; NumoProblems>0; NumoProblems>>=1) {
    int jFirst, j=0;
    for (jFirst=0; jFirst<PARAM_N; jFirst=j+NumoProblems) {
      sdigit_t W = (sdigit_t)w[jTwiddle++];
      for (j=jFirst; j<jFirst+NumoProblems; j++) {
        int32_t temp = reduce((int64_t)W * a[j+NumoProblems]);
        a[j + NumoProblems] = a[j] - temp;
        a[j] = temp + a[j];
      }
    }
  }
}

//#if !defined(_qTESLA_I_)

int32_t barr_reduce(int32_t a)
{ // Barrett reduction
  int32_t u = ((int64_t)a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - (int32_t)u*PARAM_Q;
}

//#endif

void nttinv(poly a, const poly w)
{ // Inverse NTT transform
  int NumoProblems = 1, jTwiddle=0;
  for (NumoProblems=1; NumoProblems<PARAM_N; NumoProblems*=2) {
    int jFirst, j=0;
    for (jFirst = 0; jFirst<PARAM_N; jFirst=j+NumoProblems) {
      sdigit_t W = (sdigit_t)w[jTwiddle++];
      for (j=jFirst; j<jFirst+NumoProblems; j++) {
        int32_t temp = a[j];
#if defined(_qTESLA_I_)
        a[j] = temp + a[j + NumoProblems];
#else
        if (NumoProblems == 16) 
          a[j] = barr_reduce(temp + a[j + NumoProblems]);
        else
          a[j] = temp + a[j + NumoProblems];
#endif
        a[j + NumoProblems] = reduce((int64_t)W * (temp - a[j + NumoProblems]));
      }
    }
  }

  for (int i = 0; i < PARAM_N/2; i++)
    a[i] = reduce((int64_t)PARAM_R*a[i]);
}


static void poly_pointwise(poly result, const poly x, const poly y)
{ // Pointwise polynomial multiplication result = x.y

  for (int i=0; i<PARAM_N; i++)
    result[i] = reduce((int64_t)x[i]*y[i]);
}

void poly_add(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y

    for (int i=0; i<PARAM_N; i++)
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


void poly_sub_correct(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y with correction

    for (int i=0; i<PARAM_N; i++) {
      result[i] = x[i] - y[i];
      result[i] += (result[i] >> (RADIX32-1)) & PARAM_Q;    // If result[i] < 0 then add q
    }
}


void poly_sub_reduce(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y with Montgomery reduction

    for (int i=0; i<PARAM_N; i++)
      result[i] = reduce((int64_t)PARAM_R*(x[i] - y[i]));
}


/********************************************************************************************
* Name:        sparse_mul16
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
void sparse_mul16(poly prod, const int16_t *s, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H])
{
  int i, j, pos;
  int16_t *t = (int16_t*)s;

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

  for (i=0; i<PARAM_N; i++)
    prod[i] = 0;
  
  for (i=0; i<PARAM_H; i++) {
    pos = pos_list[i];
    for (j=0; j<pos; j++) {
        prod[j] = prod[j] - sign_list[i]*pk[j+PARAM_N-pos];
    }
    for (j=pos; j<PARAM_N; j++) {
        prod[j] = prod[j] + sign_list[i]*pk[j-pos];
    }
  }
}
