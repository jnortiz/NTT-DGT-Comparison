/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include <inttypes.h>
#include <stdio.h>
#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"
#include "params_dgt.h"

#define IFLESS32(val, gauge, expr) ((((val) - (gauge)) >> (31)) & (expr))
#define IFMORE32(val, gauge, expr) ((((gauge) - (val)) >> (31)) & (expr))
#define IFLESS64(val, gauge, expr) ((((val) - (gauge)) >> (63)) & (expr))
#define IFMORE64(val, gauge, expr) ((((gauge) - (val)) >> (63)) & (expr))
#define CENTER32(v) ((v) += IFLESS32((v), -(PARAM_Q >> 1), PARAM_Q) - IFMORE32((v), (PARAM_Q >> 1), PARAM_Q))

#define xstr(s) str(s)
#define str(s) #s


#if defined(USE_REFERENCE)

static inline int64_t mulh64(int64_t a, int64_t b)
{ // Computes high word result of 64x64 product
    uint64_t al, bl, mask_low = (uint64_t)(-1) >> 32;
    int64_t ah, bh, albl, albh, ahbl, ahbh, t0, t1;

    al = a & mask_low;         // low part
    ah = a >> 32;              // high part
    bl = b & mask_low;         
    bh = b >> 32;   

    albl = (int64_t)(al*bl);
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;

    t0 = (ahbl & mask_low) + (albh & mask_low) + ((uint64_t)albl >> 32);  
    t0 = (ahbh & mask_low) + (ahbl >> 32) + (albh >> 32) + (t0 >> 32);
    t1 = (ahbh >> 32) + (t0 >> 32);

    return t0 + (t1 << 32);
}


static inline int64_t Barrett(int64_t a)
{ // Barrett reduction
  int64_t u = mulh64(a, (int64_t)RING_QREC);
  a -= (int64_t)u*PARAM_Q;
  return a + IFLESS64(a, 0, PARAM_Q) - IFMORE64(a, PARAM_Q - 1, PARAM_Q);
}

#else

static inline int64_t Barrett(int64_t a)
{ // Barrett reduction
    __asm__ (
        "movq   %1, %%rax;"                     // RAX = a
        "movq   $" xstr(RING_QREC) ", %%rdx;"   // RDX = x = floor(2^64/q)
        "imulq  %%rdx;"                         // RDX = u = (a*x) >> 64
        "movq   $-" xstr(PARAM_Q) ", %%rax;"    // RAX = -q
        "imulq  %%rdx;"                         // RAX = -u*q
        "movq   %1, %%rdx;"                     // RDX = a
        "addq   %%rdx, %%rax;"                  // RAX = a - u*q
        "movq   %%rax, %0;"
        : "=r"(a)           // output
        : "r"(a)            // input
        : "%rax", "%rdx"    // clobbered registers
    );
    // Now -q <= a < 2*q
    // NB: the following reduction is only needed at the very end of a long list of operations:
    return a + IFLESS64(a, 0, PARAM_Q) - IFMORE64(a, PARAM_Q - 1, PARAM_Q);
    // Now 0 <= a < q
}

#endif


static inline int32_t MODQ(int64_t v) 
{
    return (int32_t)Barrett(v);
}


void poly_uniform(poly a, const unsigned char *seed)
{ // Generation of polynomial "a" (in NTT form, i.e. as the sequence of eigenvalues)
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
      a[i++] = MODQ((int64_t)val1);
    if (val2 < PARAM_Q && i < PARAM_N)
      a[i++] = MODQ((int64_t)val2);
    if (val3 < PARAM_Q && i < PARAM_N)
      a[i++] = MODQ((int64_t)val3);
    if (val4 < PARAM_Q && i < PARAM_N)
      a[i++] = MODQ((int64_t)val4);
  }
}

void dgt_6(poly _x) 
{

  size_t i, j, k;
  int32_t *temp;
  int32_t s0_re, s0_img;
  int32_t s1_re, s1_img;
  int32_t s2_re, s2_img;
  int32_t s3_re, s3_img;

  temp = _x;

  // For each of the PARAM_M/2 blocks containing six Gaussian integers
  for(j = 0; j < PARAM_M; j += 2) {
    
    // Extract a block
    int32_t block[12] = {
      temp[0 << PARAM_LGM], temp[(0 << PARAM_LGM)+1],
      temp[1 << PARAM_LGM], temp[(1 << PARAM_LGM)+1], 
      temp[2 << PARAM_LGM], temp[(2 << PARAM_LGM)+1],
      temp[3 << PARAM_LGM], temp[(3 << PARAM_LGM)+1], 
      temp[4 << PARAM_LGM], temp[(4 << PARAM_LGM)+1], 
      temp[5 << PARAM_LGM], temp[(5 << PARAM_LGM)+1],
    };

    // Apply the \psi transformation in each Gaussian integer
    k = 0;
    for(i = 0; i < 12; i += 2) {
      block[i] = MODQ((int64_t)block[i] * _theta[k]); // real part
      block[i+1] = MODQ((int64_t)block[i+1] * _theta[k]); // imaginary part
      k++;
    }

    // And then the \mu transformation
    s0_re = block[0] + block[6];
    s0_img = block[1] + block[7];

    s1_re = block[2] + block[8];
    s1_img = block[3] + block[9];

    s2_re = block[4] + block[10];
    s2_img = block[5] + block[11];

    s3_re = MODQ((int64_t)(s1_re - s2_re) * _theta[3]);
    s3_img = MODQ((int64_t)(s1_img - s2_img)*_theta[3]);

    temp[0 << PARAM_LGM] = (s0_re + s1_re + s2_re);
    temp[(0 << PARAM_LGM)+1] = (s0_img + s1_img + s2_img);

    temp[2 << PARAM_LGM] = (s0_re - s2_re + s3_re);
    temp[(2 << PARAM_LGM)+1] = (s0_img - s2_img + s3_img);

    temp[4 << PARAM_LGM] = (s0_re - s3_re - s1_re);
    temp[(4 << PARAM_LGM)+1] = (s0_img - s3_img - s1_img);

    s0_re = block[0] + MODQ((int64_t)block[6] * _theta[3]);
    s0_img = block[1] + MODQ((int64_t)block[7] * _theta[3]);

    s1_re = MODQ((int64_t)block[2] * _theta[1]) + MODQ((int64_t)block[8] * _theta[4]);
    s1_img = MODQ((int64_t)block[3] * _theta[1]) + MODQ((int64_t)block[9] * _theta[4]);

    s2_re = MODQ((int64_t)block[4] * _theta[2]) + MODQ((int64_t)block[10] * _theta[5]);
    s2_img = MODQ((int64_t)block[5] * _theta[2]) + MODQ((int64_t)block[11] * _theta[5]);

    s3_re = MODQ((int64_t)(s1_re - s2_re) * _theta[3]);
    s3_img = MODQ((int64_t)(s1_img - s2_img)*_theta[3]);

    temp[1 << PARAM_LGM] = (s0_re + s1_re + s2_re);
    temp[(1 << PARAM_LGM)+1] = (s0_img + s1_img + s2_img);

    temp[3 << PARAM_LGM] = (s0_re - s2_re + s3_re);
    temp[(3 << PARAM_LGM)+1] = (s0_img - s2_img + s3_img);

    temp[5 << PARAM_LGM] = (s0_re - s3_re - s1_re);
    temp[(5 << PARAM_LGM)+1] = (s0_img - s3_img - s1_img);

    temp += 2;

  }
}

void idgt_6(poly _output_signal) {
  
  size_t i, j, k;
  int32_t *temp;
  int32_t s0_re, s0_img;
  int32_t s1_re, s1_img;
  int32_t t0_re, t0_img;
  int32_t t1_re, t1_img;

  temp = _output_signal;

  for(j = 0; j < PARAM_M; j += 2) {

    int32_t block[12] = {
      temp[0 << PARAM_LGM], temp[(0 << PARAM_LGM)+1],
      temp[1 << PARAM_LGM], temp[(1 << PARAM_LGM)+1], 
      temp[2 << PARAM_LGM], temp[(2 << PARAM_LGM)+1],
      temp[3 << PARAM_LGM], temp[(3 << PARAM_LGM)+1], 
      temp[4 << PARAM_LGM], temp[(4 << PARAM_LGM)+1], 
      temp[5 << PARAM_LGM], temp[(5 << PARAM_LGM)+1],
    };

    t0_re = MODQ((int64_t)(block[4] - block[8])*_theta[3]);
    t0_img = MODQ((int64_t)(block[5] - block[9])*_theta[3]);

    t1_re = MODQ((int64_t)(block[6] - block[10])*_theta[3]);
    t1_img = MODQ((int64_t)(block[7] - block[11])*_theta[3]);

    s0_re = (block[0] + block[4] + block[8]);
    s0_img = (block[1] + block[5] + block[9]);

    s1_re = (block[2] + block[6] + block[10]);
    s1_img = (block[3] + block[7] + block[11]);

    temp[0 << PARAM_LGM] = s0_re - MODQ((int64_t)s1_re * _theta[6]);
    temp[(0 << PARAM_LGM)+1] = s0_img - MODQ((int64_t)s1_img * _theta[6]);

    temp[3 << PARAM_LGM] = (s0_re - s1_re);
    temp[(3 << PARAM_LGM)+1] = (s0_img - s1_img);

    s0_re = (block[0] - block[4] - t0_re);
    s0_img = (block[1] - block[5] - t0_img);

    s1_re = (block[2] - block[6] - t1_re);
    s1_img = (block[3] - block[7] - t1_img);    

    temp[1 << PARAM_LGM] = s0_re - MODQ((int64_t)s1_re * _theta[5]);
    temp[(1 << PARAM_LGM)+1] = s0_img - MODQ((int64_t)s1_img * _theta[5]);

    temp[4 << PARAM_LGM] = s0_re - MODQ((int64_t)s1_re * _theta[8]);
    temp[(4 << PARAM_LGM)+1] = s0_img - MODQ((int64_t)s1_img * _theta[8]);    

    s0_re = (block[0] - block[8] - t0_re);
    s0_img = (block[1] - block[9] - t0_img);

    s1_re = (block[2] - block[10] + t1_re);
    s1_img = (block[3] - block[11] + t1_img);

    temp[2 << PARAM_LGM] = s0_re - MODQ((int64_t)s1_re * _theta[4]);
    temp[(2 << PARAM_LGM)+1] = s0_img - MODQ((int64_t)s1_img * _theta[4]);

    temp[5 << PARAM_LGM] = s0_re - MODQ((int64_t)s1_re * _theta[7]);
    temp[(5 << PARAM_LGM)+1] = s0_img - MODQ((int64_t)s1_img * _theta[7]);

    temp[0 << PARAM_LGM] = MODQ((int64_t)temp[0 << PARAM_LGM] * one_minus_zeta_by_nine);
    temp[(0 << PARAM_LGM)+1] = MODQ((int64_t)temp[(0 << PARAM_LGM)+1] * one_minus_zeta_by_nine);

    temp[1 << PARAM_LGM] = MODQ((int64_t)MODQ((int64_t)temp[1 << PARAM_LGM] * one_minus_zeta_by_nine) * _theta[8]);
    temp[(1 << PARAM_LGM)+1] = MODQ((int64_t)MODQ((int64_t)temp[(1 << PARAM_LGM)+1] * one_minus_zeta_by_nine) * _theta[8]);

    temp[2 << PARAM_LGM] = MODQ((int64_t)MODQ((int64_t)temp[2 << PARAM_LGM] * one_minus_zeta_by_nine) * _theta[7]);
    temp[(2 << PARAM_LGM)+1] = MODQ((int64_t)MODQ((int64_t)temp[(2 << PARAM_LGM)+1] * one_minus_zeta_by_nine) * _theta[7]);

    temp[3 << PARAM_LGM] = MODQ((int64_t)MODQ((int64_t)temp[3 << PARAM_LGM] * zeta_minus_one_by_nine) * _theta[3]);
    temp[(3 << PARAM_LGM)+1] = MODQ((int64_t)MODQ((int64_t)temp[(3 << PARAM_LGM)+1] * zeta_minus_one_by_nine) * _theta[3]);

    temp[4 << PARAM_LGM] = MODQ((int64_t)MODQ((int64_t)temp[4 << PARAM_LGM] * zeta_minus_one_by_nine) * _theta[2]);
    temp[(4 << PARAM_LGM)+1] = MODQ((int64_t)MODQ((int64_t)temp[(4 << PARAM_LGM)+1] * zeta_minus_one_by_nine) * _theta[2]);

    temp[5 << PARAM_LGM] = MODQ((int64_t)MODQ((int64_t)temp[5 << PARAM_LGM] * zeta_minus_one_by_nine) * _theta[1]);
    temp[(5 << PARAM_LGM)+1] = MODQ((int64_t)MODQ((int64_t)temp[(5 << PARAM_LGM)+1] * zeta_minus_one_by_nine) * _theta[1]);
  
    temp += 2;

  }
}

void dgt_2(poly _input_signal)
{ // This function assumes that each PARAM_M-block is folded and twisted accordingly

  int i, index, j, k, m, window;
  int32_t a, sub_re, sub_img;

  int32_t *_x = _input_signal;

  // For each of the six blocks
  for(k = 0; k < PARAM_N; k += PARAM_M) {    
    // Apply the DGT in the ring Z_q[x]/<x^{PARAM_M}+1>
    window = 1;
    
    for(m = (PARAM_M >> 1); m > 1; m >>= 1) {
      index = 0;

      for(j = 0; j < m; j += 2) {
        a = _gj[index];

        for(i = j; i < PARAM_M; i += (m << 1)) {
          sub_re = _x[i] - _x[i+m];
          sub_img = _x[i+1] - _x[i+m+1];
          
          _x[i] += _x[i+m];
          _x[i+1] += _x[i+m+1];
          
          _x[i+m] = MODQ((int64_t)a * sub_re);
          _x[i+m+1] = MODQ((int64_t)a * sub_img);        
        }
        index += window;
      }
      window <<= 1;
    }
    _x += PARAM_M;
  }
}


void idgt_2(poly _a)
{
  int i, index, j, k, m, window;
  int32_t a, mul_re, mul_img;

  int32_t *_output_signal = _a;

  for(k = 0; k < PARAM_N; k += PARAM_M) {    

    window = (PARAM_M >> 2);
    for(m = 2; m <= (PARAM_M >> 1); m <<= 1) {
      index = 0;

      for(j = 0; j < m; j += 2) {
        a = _invgj[index];
        
        for(i = j; i < PARAM_M; i += (m << 1)) {
          mul_re = MODQ((int64_t)_output_signal[i+m] * a);
          mul_img = MODQ((int64_t)_output_signal[i+m+1] * a);
          
          _output_signal[i+m] =  _output_signal[i] - mul_re;
          _output_signal[i+m+1] = _output_signal[i+1] - mul_img;
          
          _output_signal[i] += mul_re;
          _output_signal[i+1] += mul_img;        
        }
        index += window;
      }
      window >>= 1;
    }
    _output_signal += PARAM_M;
  }
}


void poly_dgt(poly _a) {

  size_t i, j, l;
  int k;
  int32_t *_poly_a = _a;

  k = (PARAM_M >> 1);

  for(l = 0; l < PARAM_N; l += PARAM_M) {
    for(i = 0, j = 0; i < PARAM_M && j < k; i+=2, j++) {                       
        _poly_a[i] = MODQ(
          (int64_t)_poly_a[j] * _nthroots[i] -
          (int64_t)_poly_a[k+j] * _nthroots[i+1]
          );

        _poly_a[i+1] = MODQ(
          (int64_t)_poly_a[j] * _nthroots[i+1] + 
          (int64_t)_poly_a[k+j] * _nthroots[i]
          );
    }
    _poly_a += PARAM_M;
  } 

  /* Computing the DGT of signal b. The signal a is assumed to be in the DGT domain */
  dgt_2(_a);
  dgt_6(_a);

}

void poly_mul(poly _output, const poly a, const poly b)
{   // The input _poly_a is assumed to be in DGT form

    poly _poly_a, _poly_b, _folded_a, _folded_b, _mul; // Polynomials with 2^{l-1}*6 = PARAM_M*6 = PARAM_N coefficients
    size_t i, j, l;
    int k;

    k = (PARAM_M >> 1); // k = m/2 = 128/2 = 64

    for(l = 0; l < PARAM_N; ++l) {
      _poly_a[l] = 1;
      _poly_b[l] = 1;
    }


    // Twisting and folding the six blocks of PARAM_M coefficients
    for(l = 0; l < PARAM_N; l += PARAM_M) {
      for(i = 0, j = 0; i < PARAM_M && j < k; i+=2, j++) {                       
          _folded_b[l+i] = MODQ(
            (int64_t)_poly_b[l+j] * _nthroots[i] -
            (int64_t)_poly_b[l+k+j] * _nthroots[i+1]
            );

          _folded_b[l+i+1] = MODQ(
            (int64_t)_poly_b[l+j] * _nthroots[i+1] + 
            (int64_t)_poly_b[l+k+j] * _nthroots[i]
            );

          _folded_a[l+i] = MODQ(
            (int64_t)_poly_a[l+j] * _nthroots[i] -
            (int64_t)_poly_a[l+k+j] * _nthroots[i+1]
            );

          _folded_a[l+i+1] = MODQ(
            (int64_t)_poly_a[l+j] * _nthroots[i+1] + 
            (int64_t)_poly_a[l+k+j] * _nthroots[i]
            );          
      } 
    } 


    /* Computing the DGT of signal b. The signal a is assumed to be in the DGT domain */
    dgt_2(_folded_b);
    idgt_2(_folded_b);

    for(l = 0; l < PARAM_N; l += PARAM_M) {            
      for(i = 0, j = 0; i < PARAM_M && j < k; i+=2, j++) {
          _output[l+j] = MODQ(
            (int64_t)_folded_b[l+i] * _invnthroots[i] - 
            (int64_t)_folded_b[l+i+1] * _invnthroots[i+1]
            );
          _output[l+j+k] = MODQ(
            (int64_t)_folded_b[l+i] * _invnthroots[i+1] +
            (int64_t)_folded_b[l+i+1] * _invnthroots[i]
            );
      }
    }    

    printf("Depois: ");
    for(l = 0; l < PARAM_M; ++l) {
      printf("%" PRId32 ", ", _output[l]);
    }
    printf("\n\n\n");

    return;

    dgt_6(_folded_b);

    dgt_2(_folded_a);
    dgt_6(_folded_a);

    /* Calculating the point-wise multiplication of input signals */
    for(i = 0; i < PARAM_N; i += 2) {             
        _mul[i] = MODQ(
          (int64_t)_folded_a[i] *  _folded_b[i] -
          (int64_t)_folded_a[i+1] * _folded_b[i+1]
          );
        _mul[i+1] = MODQ(
          (int64_t)_folded_a[i] * _folded_b[i+1] +
          (int64_t)_folded_a[i+1] * _folded_b[i]
          );
    }

    /* Recovering the multiplication result in Z[x]/<x^n+1> */
    idgt_6(_mul);
    idgt_2(_mul);

    /* Removing the twisting factors and writing the result from the Gaussian integer to the polynomial form */
    for(l = 0; l < PARAM_N; l += PARAM_M) {            
      for(i = 0, j = 0; i < PARAM_M && j < k; i+=2, j++) {
          _output[l+j] = MODQ(
            (int64_t)_mul[l+i] * _invnthroots[i] - 
            (int64_t)_mul[l+i+1] * _invnthroots[i+1]
            );
          _output[l+j+k] = MODQ(
            (int64_t)_mul[l+i] * _invnthroots[i+1] +
            (int64_t)_mul[l+i+1] * _invnthroots[i]
            );
      }
    }

    // for(l = 0; l < PARAM_N; ++l) {
    //   printf("%" PRId32 ", ", _output[l]);
    // }
    // printf("\n\n\n");


}


//extern size_t rev_tab[PARAM_M];
//extern int32_t twiddle[PARAM_M], psi_phi_tab[PARAM_N], iph_ips_tab[PARAM_N], t_tab[9];
static size_t rev_tab[PARAM_M];
static int32_t twiddle[PARAM_M], psi_phi_tab[PARAM_N], iph_ips_tab[PARAM_N], t_tab[9];

void nttInit() 
{
    // Find roots of unity:
    int32_t gamma = 2;
    while (gamma < PARAM_Q) {
        bool ok = true;
        int32_t v = gamma; // gamma^(2^0)
        for (size_t i = 0; i <= PARAM_LGM; i++) { // v = gamma^(2^i)
            if (v == 1) {
                ok = false;
                break;
            }
            v = MODQ((int64_t)v*v);
        }
        if (ok && v == 1) {
            break;
        }
        gamma++;
    }
    int32_t theta = 2;
    while (theta < PARAM_Q) {
        int32_t v = MODQ((int64_t)MODQ((int64_t)theta*theta)*theta);
        if (v != 1 && MODQ((int64_t)MODQ((int64_t)v*v)*v) == 1) {
            break;
        }
        theta++;
    }

    // Prepare the bit-reversal table for lg(m)-bit integers
    for (size_t i = 0; i < PARAM_M; i++) {
        size_t v = i, r = 0;
        for (size_t s = 0; s < PARAM_LGM; s++) {
            r = (r << 1) + (v & 1);
            v >>= 1;
        }
        rev_tab[i] = r;
    }

    // Make a table of powers of the m-th root of unity omega in Z/qZ, 0 <= k < m.
    int32_t omega = MODQ((int64_t)gamma*gamma);
    twiddle[0] = 1;
    for (size_t k = 1; k < PARAM_M; k++) {
        twiddle[k] = MODQ((int64_t)twiddle[k - 1]*omega);
    }

    t_tab[0] = 1;
    for (size_t i = 1; i < 9; i++) {
        t_tab[i] = MODQ((int64_t)t_tab[i - 1]*theta);
    }

    // Make a table psi_phi_tab[] of size n = 6*m such that psi_phi_tab[j] = gamma^(j mod m)*theta^(j div m) % q
    psi_phi_tab[0] = 1;
    for (size_t j = 1; j < PARAM_M; j++) {
        psi_phi_tab[j] = MODQ((int64_t)psi_phi_tab[j - 1]*gamma);
    }
    for (size_t j = PARAM_M; j < PARAM_N; j++) {
        psi_phi_tab[j] = MODQ((int64_t)psi_phi_tab[j & (PARAM_M - 1)]*t_tab[j >> PARAM_LGM]);
    }

    // Make a table iph_ips_tab[] of size n = 6*m such that iph_ips_tab[j]*psi_phi_tab[j] == scale
    int32_t gamma_inv = MODQ((int64_t)MODQ((int64_t)psi_phi_tab[PARAM_M - 1]*psi_phi_tab[PARAM_M - 1])*psi_phi_tab[1]); // gamma^-1 (mod q)
    int32_t t_inv_tab[6] = { 1, t_tab[8], t_tab[7], -t_tab[3], -t_tab[2], -t_tab[1] };
    iph_ips_tab[0] = MODQ((int64_t)(1 - t_tab[3])*(PARAM_Q - (PARAM_Q - 1)/(9*PARAM_M))); // scale
    for (size_t j = 1; j < PARAM_M; j++) {
        iph_ips_tab[j] = MODQ((int64_t)iph_ips_tab[j - 1]*gamma_inv);
    }
    for (size_t j = PARAM_M; j < PARAM_N; j++) {
        iph_ips_tab[j] = MODQ((int64_t)iph_ips_tab[j & (PARAM_M - 1)]*t_inv_tab[j >> PARAM_LGM]);
    }
}


static void mu(int32_t Ai[], const int32_t ai[]) 
{ // The mu component of the sextic NTT
  // Input: ai, the vector from (Z/qZ)^6 to apply mu to.
    int32_t u[6] = {
        ai[0 << PARAM_LGM], ai[1 << PARAM_LGM], ai[2 << PARAM_LGM],
        ai[3 << PARAM_LGM], ai[4 << PARAM_LGM], ai[5 << PARAM_LGM],
    };
    int64_t s0, s1, s2, s3;
    s0 = u[0] + u[3];
    s1 = u[1] + u[4];
    s2 = u[2] + u[5];
    s3 = (s1 - s2)*t_tab[3]; //% q;
    Ai[0 << PARAM_LGM] = MODQ(s0 + s1 + s2);
    Ai[2 << PARAM_LGM] = MODQ(s0 - s2 + s3);
    Ai[4 << PARAM_LGM] = MODQ(s0 - s3 - s1);
    s0 = u[0] + (int64_t)u[3]*t_tab[3]; //% q;
    s1 = (int64_t)u[1]*t_tab[1] + (int64_t)u[4]*t_tab[4]; //% q;
    s2 = (int64_t)u[2]*t_tab[2] + (int64_t)u[5]*t_tab[5]; //% q;
    s3 = (int64_t)MODQ(s1 - s2)*t_tab[3]; //% q;
    Ai[1 << PARAM_LGM] = MODQ(s0 + s1 + s2);
    Ai[3 << PARAM_LGM] = MODQ(s0 - s2 + s3);
    Ai[5 << PARAM_LGM] = MODQ(s0 - s3 - s1);
}


static void mu_dag(int32_t ai[], const int32_t Ai[]) 
{ // The mu_dag component of the sextic NTT.
  // Input: Ai, the vector from (Z/qZ)^6 to apply mu_dag to.
    int32_t u[6] = {
        Ai[0 << PARAM_LGM], Ai[1 << PARAM_LGM], Ai[2 << PARAM_LGM],
        Ai[3 << PARAM_LGM], Ai[4 << PARAM_LGM], Ai[5 << PARAM_LGM],
    };
    int64_t s0, s1, t0, t1;
    t0 = (int64_t)(u[2] - u[4])*t_tab[3]; //% q;
    t1 = (int64_t)MODQ((int64_t)(u[3] - u[5])*t_tab[3]);
    s0 =  u[0] + u[2] + u[4];
    s1 =  u[1] + u[3] + u[5];
    ai[0 << PARAM_LGM] = MODQ(s0 - s1*t_tab[6]);
    ai[3 << PARAM_LGM] = MODQ(s0 - s1);
    s0 =  u[0] - u[2] - t0;
    s1 =  u[1] - u[3] - t1;
    ai[1 << PARAM_LGM] = MODQ(s0 - s1*t_tab[5]);
    ai[4 << PARAM_LGM] = MODQ(s0 - s1*t_tab[8]);
    s0 =  u[0] - u[4] + t0;
    s1 =  u[1] - u[5] + t1;
    ai[2 << PARAM_LGM] = MODQ(s0 - s1*t_tab[4]);
    ai[5 << PARAM_LGM] = MODQ(s0 - s1*t_tab[7]);
}


void NTT(int32_t A[], const int32_t a[]) 
{ // The number-theoretic transform on (Z/qZ)[x, y]/<x^m + 1, y^6 + y^3 + 1>.
  // Input:  a, vector from (Z/qZ)^n, n = 6*m 
  // Output: A, output vector NTT(a) in (Z/qZ)^n, n = 6*m
    int32_t Aa[PARAM_N];

    // Apply the fudge preprocessing:
    for (size_t j = 0; j < PARAM_N; j++) {
        Aa[j] = MODQ((int64_t)a[j]*psi_phi_tab[j]);
    }
    // Apply the binary NTT to each of the six m-entry blocks:
    int32_t *Ao = A;
    for (size_t o = 0; o < PARAM_N; o += PARAM_M) {
        // bit-reverse copy:
        int32_t *Aao = Aa + o;
        for (size_t v = 0; v < PARAM_M; v++) {
            Ao[rev_tab[v]] = Aao[v];
        }
        for (size_t s = 1; s <= PARAM_LGM; s++) {
            size_t t = 1 << s, NumoProblems = 1 << (s - 1);
            for (size_t jFirst = 0; jFirst < NumoProblems; jFirst++) {
                int32_t W = twiddle[jFirst << (PARAM_LGM - s)];
                for (size_t j = jFirst; j < jFirst + PARAM_M; j += t) {
                    int32_t temp = MODQ((int64_t)W*Ao[j + NumoProblems]);
                    Ao[j + NumoProblems] = Ao[j] - temp;
                    Ao[j] = Ao[j] + temp;
                }
            }
        }
        Ao += PARAM_M;
    }
    // Apply the mu transform to interlaced blocks:
    for (size_t i = 0; i < PARAM_M; i++) {
        mu(A + i, A + i);
    }
}


void invNTT(int32_t a[], const int32_t A[]) 
{ // The inverse number-theoretic transform on (Z/qZ)[x, y]/<x^m + 1, y^6 + y^3 + 1>.
  // Input:  A, vector from (Z/qZ)^n, n = 6*m 
  // Output: a, output vector NTT^{-1}(A) in (Z/qZ)^n, n = 6*m
    int32_t aA[PARAM_N];

    // Apply the mu transform to interlaced blocks:
    for (size_t i = 0; i < PARAM_M; i++) {
        mu_dag(aA + i, A + i);
    }
    // Apply the inverse binary NTT to each of the six m-entry blocks:
    int32_t *ao = a;
    for (size_t o = 0; o < PARAM_N; o += PARAM_M) {
        // bit-reverse copy:
        int32_t *aAo = aA + o;
        for (size_t v = 0; v < PARAM_M; v++) {
            ao[rev_tab[v]] = aAo[v];
        }
        for (size_t s = 1; s <= PARAM_LGM; s++) {
            size_t t = 1 << s, NumoProblems = 1 << (s - 1);
            for (size_t jFirst = 0; jFirst < NumoProblems; jFirst++) {
                int32_t W = twiddle[(-(jFirst << (PARAM_LGM - s))) & (PARAM_M - 1)];
                for (size_t k = jFirst; k < jFirst + PARAM_M; k += t) {
                    int32_t temp = MODQ((int64_t)W*ao[k + NumoProblems]);
                    ao[k + NumoProblems] = ao[k] - temp;
                    ao[k] = ao[k] + temp;
                }
            }
        }
        ao += PARAM_M;
    }
    // Apply the fudge postprocessing and centralize:
    for (size_t j = 0; j < PARAM_N; j++) {
        a[j] = MODQ((int64_t)a[j]*iph_ips_tab[j]);
        CENTER32(a[j]); aA[j] = 0;
    }
}


void ntt(poly a)
{ // Forward NTT transform
  if (twiddle[0] == 0) nttInit();
  NTT(a, a);
}


void nttinv(poly a)
{ // Inverse NTT transform
  if (twiddle[0] == 0) nttInit();
  invNTT(a, a);
}


static void poly_pointwise(poly result, const poly x, const poly y)
{ // Pointwise polynomial multiplication result = x.y

  for (int i=0; i<PARAM_N; i++)
    result[i] = MODQ((int64_t)x[i]*y[i]);
}


void _poly_mul(poly result, const poly x, const poly y)
{ // Polynomial multiplication result = x*y, with in place reduction for (X^N+1)
  // The input x is assumed to be in NTT form
  poly y_ntt;

  for (int i=0; i<PARAM_N; i++)
    y_ntt[i] = y[i];

  ntt(y_ntt);
  poly_pointwise(result, x, y_ntt);
  nttinv(result);
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
{ // Polynomial subtraction result = x-y with reduction

    for (int i=0; i<PARAM_N; i++)
      result[i] = MODQ((int64_t)(x[i] - y[i]));
}
