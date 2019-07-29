/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"
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


int64_t barr_reduce(int64_t a)
{ // Barrett reduction
  int64_t u = (a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - u*PARAM_Q;
}


void dgt(poly _x)
{

  int i, index, j, m, window;
  int64_t a, sub_re, sub_img;

  window = 1;
  m = PARAM_K2;

  for(m = PARAM_K2; m >= 4; m >>= 1) {
    index = 0;

    for(j = 0; j < m; j += 2) {
      a = _gj[index];

      for(i = j; i < PARAM_N; i += (m << 1)) {
        sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
        sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
        
        _x[i] = (_x[i] + _x[i+m]);
        _x[i+1] = (_x[i+1] + _x[i+m+1]);
        
        _x[i+m] = reduce((int64_t)a * sub_re);
        _x[i+m+1] = reduce((int64_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;
  }

  // m >>= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }
  // window <<= 1;

  // m >>= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }
  // window <<= 1;

  // m >>= 1;        

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }
  // window <<= 1;

  // m >>= 1;      

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }
  // window <<= 1;  

  // m >>= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }
  // window <<= 1;  

  // m >>= 1;      

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }  
  // window <<= 1;  

  // m >>= 1;      

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _gj[index];

  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     sub_re = (_x[i] + (2LL*PARAM_Q - _x[i+m]));
  //     sub_img = (_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
  //     _x[i] = (_x[i] + _x[i+m]);
  //     _x[i+1] = (_x[i+1] + _x[i+m+1]);
      
  //     _x[i+m] = (reduce((int64_t)a * sub_re));
  //     _x[i+m+1] = (reduce((int64_t)a * sub_img));        
  //   }
  //   index += window;
  // }  
  // window <<= 1;  

  // m >>= 1;      

  index = 0;
  for(j = 0; j < m; j += 2) {
    a = _gj[index];

    for(i = j; i < PARAM_N; i += (m << 1)) {
      sub_re = barr_reduce(_x[i] + (2LL*PARAM_Q - _x[i+m]));
      sub_img = barr_reduce(_x[i+1] + (2LL*PARAM_Q - _x[i+m+1]));
      
      _x[i] = barr_reduce(_x[i] + _x[i+m]);
      _x[i+1] = barr_reduce(_x[i+1] + _x[i+m+1]);
      
      _x[i+m] = (reduce((int64_t)a * sub_re));
      _x[i+m+1] = (reduce((int64_t)a * sub_img));        
    }
    index += window;
  }  

}

void idgt(poly _output_signal)
{

  int i, index, j, m, window;
  int64_t a, mul_re, mul_img;

  window = (PARAM_K2 >> 1);

  for(m = 2; m <= PARAM_K2; m <<= 1) {
    index = 0;
    for(j = 0; j < m; j += 2) {
      a = _invgj[index];
      
      for(i = j; i < PARAM_N; i += (m << 1)) {
        mul_re = reduce((int64_t)_output_signal[i+m] * a);
        mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
        
        _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
        _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
        
        _output_signal[i] = (_output_signal[i] + mul_re);
        _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
      }
      index += window;
    }
    window >>= 1; 
  }

  // m = 2;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }
  // window >>= 1; 

  // m <<= 1;

  // index = 0;
  // for(j = 0; j < m; j += 2) {
  //   a = _invgj[index];
    
  //   for(i = j; i < PARAM_N; i += (m << 1)) {
  //     mul_re = reduce((int64_t)_output_signal[i+m] * a);
  //     mul_img = reduce((int64_t)_output_signal[i+m+1] * a);
      
  //     _output_signal[i+m] = (_output_signal[i] + (2LL*PARAM_Q - mul_re));
  //     _output_signal[i+m+1] = (_output_signal[i+1] + (2LL*PARAM_Q - mul_img));
      
  //     _output_signal[i] = (_output_signal[i] + mul_re);
  //     _output_signal[i+1] = (_output_signal[i+1] + mul_img);        
  //   }
  //   index += window;
  // }            

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


void poly_dgt(poly x_dgt, const poly x)
{

  int i, j;

  for(i = 0, j = 0; i < PARAM_N && j < PARAM_K2; i+=2, j++) {             
      x_dgt[i] = reduce(
        ((int64_t)x[j] * _nthroots[i]) + 
        (2LL*PARAM_Q - ((int64_t)x[PARAM_K2+j] * _nthroots[i+1])));
      
      x_dgt[i+1] = reduce(
        ((int64_t)x[j] * _nthroots[i+1]) + 
        ((int64_t)x[PARAM_K2+j] * _nthroots[i]));
  } 

  dgt(x_dgt);
}


void poly_mul_ntt(poly result, const poly x, const poly y)
{ // Polynomial multiplication result = x*y, with in place reduction for (X^N+1)
  // The inputs x and y are assumed to be in NTT form
    
  poly_pointwise(result, x, y);
  nttinv(result, zetainv);
}

void poly_mul(poly _output, const poly _poly_a, const poly _poly_b)
{ /* It is assumed that both signals are already in the DGT domain. 
     The DGT counterpart of poly_b was computed in sign.c. */
  
  poly _mul;
  int i, j;

  /* Calculating the point-wise multiplication of input signals */
  for(i = 0, j = 0; i < PARAM_N && j < PARAM_K2; i+=2, j++) {             
    _mul[i] = reduce(
      ((int64_t)_poly_a[j] * _poly_b[i]) +
      (2LL*PARAM_Q - ((int64_t)_poly_a[j+PARAM_K2] * _poly_b[i+1])));

    _mul[i+1] = reduce(
      ((int64_t)_poly_a[j] * _poly_b[i+1]) + 
      ((int64_t)_poly_a[j+PARAM_K2] * _poly_b[i]));
  }

  /* Recovering the multiplication result in Z[x]/<x^n+1> */
  idgt(_mul);

  /* Removing the twisting factors and writing the result from the Gaussian integer to the polynomial form */
  for(i = 0, j = 0; i < PARAM_N && j < PARAM_K2; i+=2, j++) {
      _output[j] = barr_reduce(
        reduce((int64_t)_mul[i] * _invnthroots[i]) +
        (2LL*PARAM_Q - (reduce((int64_t)_mul[i+1] * _invnthroots[i+1])))
      );
      _output[j+PARAM_K2] = barr_reduce(
        reduce((int64_t)_mul[i] * _invnthroots[i+1]) + 
        reduce((int64_t)_mul[i+1] * _invnthroots[i])
      );
  }
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