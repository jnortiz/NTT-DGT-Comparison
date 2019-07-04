/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include <stdio.h>
#include <inttypes.h>
#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"
#include "gaussian_integer.h"
#include "params_dgt.h"

extern poly zeta;
extern poly zetainv;

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


int64_t reduce(int64_t a)
{ // Montgomery reduction
  int64_t u;

  u = (a*PARAM_QINV) & 0xFFFFFFFF;
  u *= PARAM_Q;
  a += u;
  return a>>32;
}


void dgt(gauss_t *_x, const gauss_t *_input_signal)
{    
    int i, j, k, l, m, stride;
    gauss_t xi, xim, aux_sub, aux_power;

    k = (int)(PARAM_N/2);

    for(i = 0; i < k; i++) {
        set_gauss(&_x[i], _input_signal[i].re, _input_signal[i].img);
    }
   
    for(stride = 0; stride < UPPERBOUND; stride++) {        
        m = k/(2 << stride);
        
        for(l = 0; l < k/2; l++) {            
            j = (2*m*l)/k;
            i = j + (l%(k/(2*m)))*2*m;

            set_gauss(&xi, _x[i].re, _x[i].img);
            set_gauss(&xim, _x[i+m].re, _x[i+m].img);
            set_gauss(&aux_power, __gj[j][stride], (int32_t) 0);
            add(&_x[i], xi, xim);
            sub(&aux_sub, xi, xim);
            mul(&_x[i+m], aux_power, aux_sub);
        }
    }
}

void idgt(gauss_t *_output_signal, const gauss_t *_x)
{
    int i, j, k, l, m, stride;
    gauss_t xi, xim, aux_inv, aux_mul, aux_power;

    k = PARAM_N >> 1;

    for(i = 0; i < k; i++) {
        set_gauss(&_output_signal[i], _x[i].re, _x[i].img);
    }

    m = 1;
    for(stride = 0; stride < UPPERBOUND; stride++) {
        for(l = 0; l < k >> 1; l++) {
            j = (m*l << 1)/k;
            i = j + (l % (k/(m << 1)))*(m << 1);

            set_gauss(&xi, _output_signal[i].re, _output_signal[i].img);
            set_gauss(&xim, _output_signal[i+m].re, _output_signal[i+m].img);
            set_gauss(&aux_power, __invgj[j][stride], (int32_t) 0);
            mul(&aux_mul, aux_power, xim);

            add(&_output_signal[i], xi, aux_mul);
            sub(&_output_signal[i+m], xi, aux_mul);
        }
        m = m << 1;
    }

    set_gauss(&aux_inv, invofkmodp, (int32_t) 0);

    for(i = 0; i < k; i++) {
        mul(&_output_signal[i], _output_signal[i], aux_inv);
    }
}

void print_signal(const int64_t *_signal_a, const int length) {    
    int i;

    printf("{"); 
    for(i = 0; i < length; i++) {
        printf("%" PRId64 ", ", _signal_a[i]); 
    }
    printf("}\n"); 
}

void poly_mul(int64_t *output, const int64_t * _poly_a, const int64_t *_poly_b)
{
    int k;
    k = PARAM_N >> 1;

    gauss_t _folded_a[k], _folded_b[k];
    gauss_t _dgt_a[k], _dgt_b[k];
    gauss_t _mul[k], _output_gaussian[k];
    gauss_t root;
    int i;

    for(i = 0; i < k; i++) {
        set_gauss(&root, __nthroots[i][0], __nthroots[i][1]);
        set_gauss(&_folded_a[i], barr_reduce(_poly_a[i]), barr_reduce(_poly_a[k+i]));
        mul(&_folded_a[i], _folded_a[i], root);
        set_gauss(&_folded_b[i], barr_reduce(_poly_b[i]), barr_reduce(_poly_b[k+i]));
        mul(&_folded_b[i], _folded_b[i], root);
    }

    dgt(_dgt_a, _folded_a);
    dgt(_dgt_b, _folded_b);

    for(i = 0; i < k; i++) {
        mul(&_mul[i], _dgt_a[i], _dgt_b[i]);
    }

    idgt(_output_gaussian, _mul);

    for(i = 0; i < k; i++) {
        set_gauss(&root, __invnthroots[i][0], __invnthroots[i][1]);
        mul(&_output_gaussian[i], _output_gaussian[i], root);
        output[i] = _output_gaussian[i].re;
        output[i+k] = _output_gaussian[i].img;
    }
}


int64_t barr_reduce(int64_t a)
{ // Barrett reduction
  int64_t u = (a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - u*PARAM_Q;
}


void ntt(poly a, const poly w)
{ // Forward NTT transform
  int NumoProblems = PARAM_N>>1, jTwiddle=0;

  for (; NumoProblems>0; NumoProblems>>=1) {
    int jFirst, j=0;
    for (jFirst=0; jFirst<PARAM_N; jFirst=j+NumoProblems) {
      sdigit_t W = (sdigit_t)w[jTwiddle++];
      for (j=jFirst; j<jFirst+NumoProblems; j++) {
#if defined(_qTESLA_p_I_)
        int64_t temp = reduce((int64_t)W * a[j+NumoProblems]);
        a[j + NumoProblems] = a[j] + (PARAM_Q - temp);
        a[j] = temp + a[j];
#else
        int64_t temp = barr_reduce(reduce((int64_t)W* a[j + NumoProblems]));
        a[j + NumoProblems] = barr_reduce(a[j] +(2LL*PARAM_Q - temp));
        a[j] = barr_reduce(temp + a[j]);
#endif
      }
    }
  }
}


void nttinv(poly a, const poly w)
{ // Inverse NTT transform
  int NumoProblems = 1, jTwiddle=0;
  for (NumoProblems=1; NumoProblems<PARAM_N; NumoProblems*=2) {
    int jFirst, j=0;
    for (jFirst = 0; jFirst<PARAM_N; jFirst=j+NumoProblems) {
      sdigit_t W = (sdigit_t)w[jTwiddle++];
      for (j=jFirst; j<jFirst+NumoProblems; j++) {
        int64_t temp = a[j];
#if defined(_qTESLA_p_I_)
        a[j] = (temp + a[j + NumoProblems]);
        a[j + NumoProblems] = reduce((int64_t)W * (temp + (2*PARAM_Q - a[j + NumoProblems])));
      }
    }
    NumoProblems*=2;
    for (jFirst = 0; jFirst<PARAM_N; jFirst=j+NumoProblems) {
      sdigit_t W = (sdigit_t)w[jTwiddle++];
      for (j=jFirst; j<jFirst+NumoProblems; j++) {
        int64_t temp = a[j];
        a[j] = barr_reduce(temp + a[j + NumoProblems]);
        a[j + NumoProblems] = reduce((int64_t)W * (temp + (2*PARAM_Q - a[j + NumoProblems])));
#else
        a[j] = barr_reduce((temp + a[j + NumoProblems]));
        a[j + NumoProblems] = barr_reduce(reduce((int64_t)W * (temp + (2LL*PARAM_Q - a[j + NumoProblems]))));
#endif
      }
    }
  }
}


void poly_pointwise(poly result, const poly x, const poly y)
{ // Pointwise polynomial multiplication result = x.y

  for (int i=0; i<PARAM_N; i++)
    result[i] = reduce(x[i]*y[i]);
}


void poly_ntt(poly x_ntt, const poly x)
{ // Call to NTT function. Avoids input destruction 

  for (int i=0; i<PARAM_N; i++)
    x_ntt[i] = x[i];
  ntt(x_ntt, zeta);
}

// void poly_mul(poly result, const poly x, const poly y)
// { // Polynomial multiplication result = x*y, with in place reduction for (X^N+1)
//   // The inputs x and y are assumed to be in NTT form
    
//   poly_pointwise(result, x, y);
//   nttinv(result, zetainv);
// }


void poly_add(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y

    for(int i=0; i<PARAM_N; i++)
      result[i] = x[i] + y[i];
}


void poly_add_correct(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y with correction

    for (int i=0; i<PARAM_N; i++) {
      result[i] = x[i] + y[i];
      result[i] -= PARAM_Q;
      result[i] += (result[i] >> (RADIX32-1)) & PARAM_Q;    // If result[i] >= q then subtract q
    }
}


void poly_sub(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y

    for (int i=0; i<PARAM_N; i++)
      result[i] = barr_reduce(x[i] - y[i]);
}


/********************************************************************************************
* Name:        sparse_mul8
* Description: performs sparse polynomial multiplication
* Parameters:  inputs:
*              - const unsigned char* s: part of the secret key
*              - const uint32_t pos_list[PARAM_H]: list of indexes of nonzero elements in c
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
  for (i=0; i<PARAM_N; i++)
    prod[i] = barr_reduce(prod[i]);
}