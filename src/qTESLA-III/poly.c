/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"
#include "gaussian_integer.h"
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


// void poly_mul(int32_t *output, const int32_t * _poly_a, const int32_t *_poly_b)
// {
//     int k;
//     k = PARAM_N >> 1;

//     gauss_t _folded_a[k], _folded_b[k];
//     gauss_t _dgt_a[k], _dgt_b[k];
//     gauss_t _mul[k], _output_gaussian[k];
//     gauss_t root;
//     int i;

//     for(i = 0; i < k; i++) {
//         set_gauss(&root, __nthroots[i][0], __nthroots[i][1]);
//         set_gauss(&_folded_a[i], _poly_a[i], _poly_a[k+i]);
//         mul(&_folded_a[i], _folded_a[i], root);
//         set_gauss(&_folded_b[i], _poly_b[i], _poly_b[k+i]);
//         mul(&_folded_b[i], _folded_b[i], root);
//     }

//     dgt(_dgt_a, _folded_a);
//     dgt(_dgt_b, _folded_b);

//     for(i = 0; i < k; i++) {
//         mul(&_mul[i], _dgt_a[i], _dgt_b[i]);
//     }

//     idgt(_output_gaussian, _mul);

//     for(i = 0; i < k; i++) {
//         set_gauss(&root, __invnthroots[i][0], __invnthroots[i][1]);
//         mul(&_output_gaussian[i], _output_gaussian[i], root);
//         output[i] = _output_gaussian[i].re;
//         output[i+k] = _output_gaussian[i].img;
//     }
// }


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

#if !defined(_qTESLA_I_)

int32_t barr_reduce(int32_t a)
{ // Barrett reduction
  int32_t u = ((int64_t)a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - (int32_t)u*PARAM_Q;
}

#endif

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


void poly_mul(poly result, const poly x, const poly y)
{ // Polynomial multiplication result = x*y, with in place reduction for (X^N+1)
  // The input x is assumed to be in NTT form
  poly y_ntt;
    
  for (int i=0; i<PARAM_N; i++)
    y_ntt[i] = y[i];
  
  ntt(y_ntt, zeta);
  poly_pointwise(result, x, y_ntt);
  nttinv(result, zetainv);
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
