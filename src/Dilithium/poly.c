#include <stdint.h>
#include "test/cpucycles.h"
#include "params.h"
#include "symmetric.h"
#include "ntt.h"
#include "dgt.h"
#include "reduce.h"
#include "rounding.h"
#include "poly.h"

#ifdef DBENCH
extern const unsigned long long timing_overhead;
extern unsigned long long *tred, *tadd, *tmul, *tround, *tsample, *tpack;
#endif

static uint32_t nthroots[N] = {1, 0, 4444447, 5582209, 3402799, 2789044, 2735009, 5126828, 6903475, 4110000, 1707161, 4598158, 5878935, 7873703, 17974, 5012011, 7021488, 3835807, 2343650, 3710706, 5724507, 3745254, 4076560, 6071551, 6373184, 4138213, 3591409, 2345055, 19707, 1186345, 647099, 6414465, 3459653, 1343675, 5799730, 1533425, 6980770, 5488648, 69953, 1751674, 7215940, 5966514, 951145, 3690197, 4515765, 8222245, 7090491, 5875560, 3716310, 3120684, 7950546, 6851220, 7950045, 1441678, 6033113, 204350, 4019039, 4255387, 7898825, 2596014, 6830322, 5356789, 2710297, 6905117, 158852, 5206831, 6487343, 6804546, 4792887, 3214688, 4878768, 2321613, 1177025, 7037859, 1366297, 6225046, 1434700, 1168313, 5578680, 6678876, 6318369, 389684, 684138, 4165506, 5990791, 2072748, 5775919, 6250265, 3349367, 7479494, 1635156, 5492297, 3103407, 4031630, 923050, 6071654, 3735851, 1062415, 4190134, 3874566, 3835411, 3448631, 5250364, 5697886, 4607170, 6178426, 3210158, 3708788, 5900922, 5091275, 1775935, 512221, 1694455, 8098168, 400963, 30128, 2445100, 730463, 8068545, 5731169, 7958458, 6418974, 6033751, 1915146, 7685103, 2213633, 5641097, 1704551, 4616857, 3763560, 6841430, 5556766, 2628318, 1049673, 10087, 3561574, 3727680, 1655743, 3324495, 5569286, 1791610, 8287447, 7701372, 6704204, 2460517, 7720808, 2734790, 509867, 5817598, 3205226, 5126478, 7055267, 4741034, 6153026, 3552662, 7535332, 2672922, 6245862, 3874669, 5154562, 2273862, 1581592, 5931782, 1372982, 2165633, 6764288, 7125155, 1769988, 387962, 4216192, 4488940, 5178657, 2717487, 4585631, 3616914, 1253739, 4728585, 4825647, 4226731, 1009328, 6418696, 3301823, 7004614, 1281234, 4019177, 4342334, 4986702, 3867514, 3850876, 5668564, 3469305, 186667, 3042578, 8352509, 7553661, 3020221, 4595457, 938439, 1342346, 2432336, 1558117, 732057, 873235, 2213815, 2921299, 5791593, 2738888, 1436519, 5491743, 891096, 1902502, 467319, 5744078, 2039479, 1153989, 319248, 5980884, 8299063, 3208859, 1886691, 4070121, 4155072, 5095427, 2791322, 4726741, 7576956, 2448175, 3165674, 7176710, 8064734, 2753916, 3549506, 4901974, 7698852, 6326733, 2856010, 1492508, 1920452, 3548836, 825583, 5491376, 5468251, 3506823, 1886254, 8125258, 246878, 7217264, 6517413, 4910583, 5023418, 2020353, 5001009, 3240814, 4208723, 1916072, 6305306};
static uint32_t invnthroots[N] = {1, 0, 6305306, 6464345, 4208723, 5139603, 5001009, 6360064, 5023418, 3469834, 6517413, 1163153, 246878, 255159, 1886254, 4873594, 5468251, 2889041, 825583, 4831581, 1920452, 6887909, 2856010, 2053684, 7698852, 3478443, 3549506, 5626501, 8064734, 1203707, 3165674, 5932242, 7576956, 3653676, 2791322, 3284990, 4155072, 4310296, 1886691, 5171558, 8299063, 2399533, 319248, 7226428, 2039479, 2636339, 467319, 6477915, 891096, 2888674, 1436519, 5641529, 5791593, 5459118, 2213815, 7507182, 732057, 6822300, 2432336, 7038071, 938439, 3784960, 3020221, 826756, 8352509, 5337839, 186667, 4911112, 5668564, 4529541, 3867514, 3393715, 4342334, 4361240, 1281234, 1375803, 3301823, 1961721, 1009328, 4153686, 4825647, 3651832, 1253739, 4763503, 4585631, 5662930, 5178657, 3891477, 4216192, 7992455, 1769988, 1255262, 6764288, 6214784, 1372982, 2448635, 1581592, 6106555, 5154562, 4505748, 6245862, 5707495, 7535332, 4827755, 6153026, 3639383, 7055267, 3253939, 3205226, 2562819, 509867, 5645627, 7720808, 5919900, 6704204, 679045, 8287447, 6588807, 5569286, 5055922, 1655743, 4652737, 3561574, 8370330, 1049673, 5752099, 5556766, 1538987, 3763560, 3763560, 1704551, 2739320, 2213633, 695314, 1915146, 2346666, 6418974, 421959, 5731169, 311872, 730463, 5935317, 30128, 7979454, 8098168, 6685962, 512221, 6604482, 5091275, 2479495, 3708788, 5170259, 6178426, 3773247, 5697886, 3130053, 3448631, 4545006, 3874566, 4190283, 1062415, 4644566, 6071654, 7457367, 4031630, 5277010, 5492297, 6745261, 7479494, 5031050, 6250265, 2604498, 2072748, 2389626, 4165506, 7696279, 389684, 2062048, 6678876, 2801737, 1168313, 6945717, 6225046, 7014120, 7037859, 7203392, 2321613, 3501649, 3214688, 3587530, 6804546, 1893074, 5206831, 8221565, 6905117, 5670120, 5356789, 1550095, 2596014, 481592, 4255387, 4361378, 204350, 2347304, 1441678, 430372, 6851220, 429871, 3120684, 4664107, 5875560, 1289926, 8222245, 3864652, 3690197, 7429272, 5966514, 1164477, 1751674, 8310464, 5488648, 1399647, 1533425, 2580687, 1343675, 4920764, 6414465, 7733318, 1186345, 8360710, 2345055, 4789008, 4138213, 2007233, 6071551, 4303857, 3745254, 2655910, 3710706, 6036767, 3835807, 1358929, 5012011, 8362443, 7873703, 2501482, 4598158, 6673256, 4110000, 1476942, 5126828, 5645408, 2789044, 4977618, 5582209, 3935970};

/*************************************************
* Name:        poly_reduce
*
* Description: Reduce all coefficients of input polynomial to representative
*              in [0,2*Q[.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = reduce32(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_csubq
*
* Description: For all coefficients of input polynomial subtract Q if
*              coefficient is bigger than Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_csubq(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = csubq(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_freeze
*
* Description: Reduce all coefficients of the polynomial to standard
*              representatives.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_freeze(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = freeze(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_add
*
* Description: Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
**************************************************/
void poly_add(poly *c, const poly *a, const poly *b)  {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = a->coeffs[i] + b->coeffs[i];

  DBENCH_STOP(*tadd);
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract polynomials. Assumes coefficients of second input
*              polynomial to be less than 2*Q. No modular reduction is
*              performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first input polynomial
**************************************************/
void poly_sub(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = a->coeffs[i] + 2*Q - b->coeffs[i];

  DBENCH_STOP(*tadd);
}

/*************************************************
* Name:        poly_shiftl
*
* Description: Multiply polynomial by 2^D without modular reduction. Assumes
*              input coefficients to be less than 2^{32-D}.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_shiftl(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] <<= D;

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_ntt
*
* Description: Forward NTT. Output coefficients can be up to 16*Q larger than
*              input coefficients.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt(poly *a) {
  DBENCH_START();

  ntt(a->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_dgt
*
* Description: Forward DGT.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_dgt(poly *a)
{
  uint32_t r, s, t, u;
  unsigned int i, j;
  DBENCH_START();

  j = 0;
  for(i = 0; i < N; i += 2) 
  {
    r = montgomery_reduce((uint64_t)a->coeffs[j]    * nthroots[i]);
    s = montgomery_reduce((uint64_t)a->coeffs[j+K2]  * nthroots[i+1]);
    t = montgomery_reduce((uint64_t)a->coeffs[j]    * nthroots[i+1]);
    u = montgomery_reduce((uint64_t)a->coeffs[j+K2]  * nthroots[i]);

    a->coeffs[i] = r + 2*Q - s;
    a->coeffs[i+1] = t + u;

    j++;
  } 

  dgt(a->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_invntt_montgomery
*
* Description: Inverse NTT and multiplication with 2^{32}. Input coefficients
*              need to be less than 2*Q. Output coefficients are less than 2*Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_montgomery(poly *a) {
  DBENCH_START();

  invntt_frominvmont(a->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_invdgt_montgomery
*
* Description: Inverse DGT and multiplication with 2^{32}.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invdgt_montgomery(poly *a) {
  unsigned int i;
  uint16_t r, s, t, u;
  const uint32_t f = (((uint64_t)MONT*MONT % Q) * (Q-1) % Q) * ((Q-1) >> 8) % Q;
  DBENCH_START();

  idgt(a->coeffs);
 
  for(i = 0; i < N; i += 2) 
  {
    r = montgomery_reduce((uint32_t)a->coeffs[i]    * invnthroots[i]);
    s = montgomery_reduce((uint32_t)a->coeffs[i+1]  * invnthroots[i+1]);
    t = montgomery_reduce((uint32_t)a->coeffs[i]    * invnthroots[i+1]);
    u = montgomery_reduce((uint32_t)a->coeffs[i+1]  * invnthroots[i]);

    a->coeffs[i] = r + 3*Q - s;
    a->coeffs[i+1] = t + u;
  }

  for(i = 0; i < N; ++i) {
    a->coeffs[i] = montgomery_reduce((uint64_t)f * a->coeffs[i]);
  }  

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_pointwise_invmontgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              with 2^{-32}. Output coefficients are less than 2*Q if input
*              coefficient are less than 22*Q.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_invmontgomery(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = montgomery_reduce((uint64_t)a->coeffs[i] * b->coeffs[i]);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_mul_invmontgomery
*
* Description: Pointwise multiplication of polynomials in DGT domain
*              representation and multiplication of resulting polynomial
*              with 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_mul_invmontgomery(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; i += 2) 
  {
    c->coeffs[i]    = montgomery_reduce((uint32_t)a->coeffs[i] * b->coeffs[i]) + 
                    (3*Q - montgomery_reduce((uint32_t)a->coeffs[i+1] * b->coeffs[i+1]));

    c->coeffs[i+1]  = montgomery_reduce((uint32_t)a->coeffs[i] * b->coeffs[i+1]) + 
                    montgomery_reduce((uint32_t)a->coeffs[i+1] * b->coeffs[i]);
  }

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_power2round
*
* Description: For all coefficients c of the input polynomial,
*              compute c0, c1 such that c mod Q = c1*2^D + c0
*              with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients Q + a0
*              - const poly *v: pointer to input polynomial
**************************************************/
void poly_power2round(poly *a1, poly *a0, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a1->coeffs[i] = power2round(a->coeffs[i], &a0->coeffs[i]);

  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        poly_decompose
*
* Description: For all coefficients c of the input polynomial,
*              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
*              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
*              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients Q + a0
*              - const poly *c: pointer to input polynomial
**************************************************/
void poly_decompose(poly *a1, poly *a0, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a1->coeffs[i] = decompose(a->coeffs[i], &a0->coeffs[i]);

  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        poly_make_hint
*
* Description: Compute hint polynomial. The coefficients of which indicate
*              whether the low bits of the corresponding coefficient of
*              the input polynomial overflow into the high bits.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
**************************************************/
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1) {
  unsigned int i, s = 0;
  DBENCH_START();

  for(i = 0; i < N; ++i) {
    h->coeffs[i] = make_hint(a0->coeffs[i], a1->coeffs[i]);
    s += h->coeffs[i];
  }

  DBENCH_STOP(*tround);
  return s;
}

/*************************************************
* Name:        poly_use_hint
*
* Description: Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *a: pointer to output polynomial with corrected high bits
*              - const poly *b: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/
void poly_use_hint(poly *a, const poly *b, const poly *h) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = use_hint(b->coeffs[i], h->coeffs[i]);

  DBENCH_STOP(*tround);
}

/*************************************************
* Name:        poly_chknorm
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients to be standard representatives.
*
* Arguments:   - const poly *a: pointer to polynomial
*              - uint32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B and 1 otherwise.
**************************************************/
int poly_chknorm(const poly *a, uint32_t B) {
  unsigned int i;
  int32_t t;
  DBENCH_START();

  /* It is ok to leak which coefficient violates the bound since
     the probability for each coefficient is independent of secret
     data but we must not leak the sign of the centralized representative. */
  for(i = 0; i < N; ++i) {
    /* Absolute value of centralized representative */
    t = (Q-1)/2 - a->coeffs[i];
    t ^= (t >> 31);
    t = (Q-1)/2 - t;

    if((uint32_t)t >= B) {
      DBENCH_STOP(*tsample);
      return 1;
    }
  }

  DBENCH_STOP(*tsample);
  return 0;
}

/*************************************************
* Name:        rej_uniform
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling using array of random bytes.
*
* Arguments:   - uint32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const unsigned char *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_uniform(uint32_t *a,
                                unsigned int len,
                                const unsigned char *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    t  = buf[pos++];
    t |= (uint32_t)buf[pos++] << 8;
    t |= (uint32_t)buf[pos++] << 16;
    t &= 0x7FFFFF;

    if(t < Q)
      a[ctr++] = t;
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform
*
* Description: Sample polynomial with uniformly random coefficients
*              in [0,Q-1] by performing rejection sampling using the
*              output stream from SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const unsigned char seed[]: byte array with seed of length
*                                            SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
void poly_uniform(poly *a,
                  const unsigned char seed[SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int nblocks = (769 + STREAM128_BLOCKBYTES)/STREAM128_BLOCKBYTES;
  unsigned int buflen = nblocks*STREAM128_BLOCKBYTES;
  unsigned char buf[buflen + 2];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, nblocks, &state);

  ctr = rej_uniform(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    off = buflen % 3;
    for(i = 0; i < off; ++i)
      buf[i] = buf[buflen - off + i];

    buflen = STREAM128_BLOCKBYTES + off;
    stream128_squeezeblocks(buf + off, 1, &state);
    ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen);
  }
}

/*************************************************
* Name:        rej_eta
*
* Description: Sample uniformly random coefficients in [-ETA, ETA] by
*              performing rejection sampling using array of random bytes.
*
* Arguments:   - uint32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const unsigned char *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_eta(uint32_t *a,
                            unsigned int len,
                            const unsigned char *buf,
                            unsigned int buflen)
{
#if ETA > 7
#error "rej_eta() assumes ETA <= 7"
#endif
  unsigned int ctr, pos;
  uint32_t t0, t1;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos < buflen) {
#if ETA <= 3
    t0 = buf[pos] & 0x07;
    t1 = buf[pos++] >> 5;
#else
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;
#endif

    if(t0 <= 2*ETA)
      a[ctr++] = Q + ETA - t0;
    if(t1 <= 2*ETA && ctr < len)
      a[ctr++] = Q + ETA - t1;
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform_eta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-ETA,ETA] by performing rejection sampling using the
*              output stream from SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const unsigned char seed[]: byte array with seed of length
*                                            SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
void poly_uniform_eta(poly *a,
                      const unsigned char seed[SEEDBYTES],
                      uint16_t nonce)
{
  unsigned int ctr;
  unsigned int nblocks = ((N/2 * (1U << SETABITS)) / (2*ETA + 1)
                          + STREAM128_BLOCKBYTES) / STREAM128_BLOCKBYTES;
  unsigned int buflen = nblocks*STREAM128_BLOCKBYTES;
  unsigned char buf[buflen];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, nblocks, &state);

  ctr = rej_eta(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    stream128_squeezeblocks(buf, 1, &state);
    ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM128_BLOCKBYTES);
  }
}

/*************************************************
* Name:        rej_gamma1m1
*
* Description: Sample uniformly random coefficients
*              in [-(GAMMA1 - 1), GAMMA1 - 1] by performing rejection sampling
*              using array of random bytes.
*
* Arguments:   - uint32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const unsigned char *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_gamma1m1(uint32_t *a,
                                 unsigned int len,
                                 const unsigned char *buf,
                                 unsigned int buflen)
{
#if GAMMA1 > (1 << 19)
#error "rej_gamma1m1() assumes GAMMA1 - 1 fits in 19 bits"
#endif
  unsigned int ctr, pos;
  uint32_t t0, t1;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos + 5 <= buflen) {
    t0  = buf[pos];
    t0 |= (uint32_t)buf[pos + 1] << 8;
    t0 |= (uint32_t)buf[pos + 2] << 16;
    t0 &= 0xFFFFF;

    t1  = buf[pos + 2] >> 4;
    t1 |= (uint32_t)buf[pos + 3] << 4;
    t1 |= (uint32_t)buf[pos + 4] << 12;

    pos += 5;

    if(t0 <= 2*GAMMA1 - 2)
      a[ctr++] = Q + GAMMA1 - 1 - t0;
    if(t1 <= 2*GAMMA1 - 2 && ctr < len)
      a[ctr++] = Q + GAMMA1 - 1 - t1;
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform_gamma1m1
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-(GAMMA1 - 1), GAMMA1 - 1] by performing rejection
*              sampling on output stream of SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const unsigned char seed[]: byte array with seed of length
*                                            CRHBYTES
*              - uint16_t nonce: 16-bit nonce
**************************************************/
void poly_uniform_gamma1m1(poly *a,
                           const unsigned char seed[CRHBYTES],
                           uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int nblocks = (641 + STREAM256_BLOCKBYTES) / STREAM256_BLOCKBYTES;
  unsigned int buflen = nblocks * STREAM256_BLOCKBYTES;
  unsigned char buf[buflen + 4];
  stream256_state state;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, nblocks, &state);

  ctr = rej_gamma1m1(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    off = buflen % 5;
    for(i = 0; i < off; ++i)
      buf[i] = buf[buflen - off + i];

    buflen = STREAM256_BLOCKBYTES + off;
    stream256_squeezeblocks(buf + off, 1, &state);
    ctr += rej_gamma1m1(a->coeffs + ctr, N - ctr, buf, buflen);
  }
}

/*************************************************
* Name:        polyeta_pack
*
* Description: Bit-pack polynomial with coefficients in [-ETA,ETA].
*              Input coefficients are assumed to lie in [Q-ETA,Q+ETA].
*
* Arguments:   - unsigned char *r: pointer to output byte array with at least
*                                  POLETA_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyeta_pack(unsigned char *r, const poly *a) {
#if 2*ETA >= 16
#error "polyeta_pack() assumes 2*ETA < 16"
#endif
  unsigned int i;
  unsigned char t[8];
  DBENCH_START();

#if 2*ETA <= 7
  for(i = 0; i < N/8; ++i) {
    t[0] = Q + ETA - a->coeffs[8*i+0];
    t[1] = Q + ETA - a->coeffs[8*i+1];
    t[2] = Q + ETA - a->coeffs[8*i+2];
    t[3] = Q + ETA - a->coeffs[8*i+3];
    t[4] = Q + ETA - a->coeffs[8*i+4];
    t[5] = Q + ETA - a->coeffs[8*i+5];
    t[6] = Q + ETA - a->coeffs[8*i+6];
    t[7] = Q + ETA - a->coeffs[8*i+7];

    r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
    r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
  }
#else
  for(i = 0; i < N/2; ++i) {
    t[0] = Q + ETA - a->coeffs[2*i+0];
    t[1] = Q + ETA - a->coeffs[2*i+1];
    r[i] = t[0] | (t[1] << 4);
  }
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyeta_unpack
*
* Description: Unpack polynomial with coefficients in [-ETA,ETA].
*              Output coefficients lie in [Q-ETA,Q+ETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const unsigned char *a: byte array with bit-packed polynomial
**************************************************/
void polyeta_unpack(poly *r, const unsigned char *a) {
  unsigned int i;
  DBENCH_START();

#if 2*ETA <= 7
  for(i = 0; i < N/8; ++i) {
    r->coeffs[8*i+0] = a[3*i+0] & 0x07;
    r->coeffs[8*i+1] = (a[3*i+0] >> 3) & 0x07;
    r->coeffs[8*i+2] = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 0x07;
    r->coeffs[8*i+3] = (a[3*i+1] >> 1) & 0x07;
    r->coeffs[8*i+4] = (a[3*i+1] >> 4) & 0x07;
    r->coeffs[8*i+5] = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 0x07;
    r->coeffs[8*i+6] = (a[3*i+2] >> 2) & 0x07;
    r->coeffs[8*i+7] = (a[3*i+2] >> 5) & 0x07;

    r->coeffs[8*i+0] = Q + ETA - r->coeffs[8*i+0];
    r->coeffs[8*i+1] = Q + ETA - r->coeffs[8*i+1];
    r->coeffs[8*i+2] = Q + ETA - r->coeffs[8*i+2];
    r->coeffs[8*i+3] = Q + ETA - r->coeffs[8*i+3];
    r->coeffs[8*i+4] = Q + ETA - r->coeffs[8*i+4];
    r->coeffs[8*i+5] = Q + ETA - r->coeffs[8*i+5];
    r->coeffs[8*i+6] = Q + ETA - r->coeffs[8*i+6];
    r->coeffs[8*i+7] = Q + ETA - r->coeffs[8*i+7];
  }
#else
  for(i = 0; i < N/2; ++i) {
    r->coeffs[2*i+0] = a[i] & 0x0F;
    r->coeffs[2*i+1] = a[i] >> 4;
    r->coeffs[2*i+0] = Q + ETA - r->coeffs[2*i+0];
    r->coeffs[2*i+1] = Q + ETA - r->coeffs[2*i+1];
  }
#endif

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt1_pack
*
* Description: Bit-pack polynomial t1 with coefficients fitting in 9 bits.
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - unsigned char *r: pointer to output byte array with at least
*                                  POLT1_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt1_pack(unsigned char *r, const poly *a) {
#if D != 14
#error "polyt1_pack() assumes D == 14"
#endif
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    r[9*i+0]  = (a->coeffs[8*i+0] >> 0);
    r[9*i+1]  = (a->coeffs[8*i+0] >> 8) | (a->coeffs[8*i+1] << 1);
    r[9*i+2]  = (a->coeffs[8*i+1] >> 7) | (a->coeffs[8*i+2] << 2);
    r[9*i+3]  = (a->coeffs[8*i+2] >> 6) | (a->coeffs[8*i+3] << 3);
    r[9*i+4]  = (a->coeffs[8*i+3] >> 5) | (a->coeffs[8*i+4] << 4);
    r[9*i+5]  = (a->coeffs[8*i+4] >> 4) | (a->coeffs[8*i+5] << 5);
    r[9*i+6]  = (a->coeffs[8*i+5] >> 3) | (a->coeffs[8*i+6] << 6);
    r[9*i+7]  = (a->coeffs[8*i+6] >> 2) | (a->coeffs[8*i+7] << 7);
    r[9*i+8]  = (a->coeffs[8*i+7] >> 1);
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt1_unpack
*
* Description: Unpack polynomial t1 with 9-bit coefficients.
*              Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const unsigned char *a: byte array with bit-packed polynomial
**************************************************/
void polyt1_unpack(poly *r, const unsigned char *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    r->coeffs[8*i+0] = ((a[9*i+0] >> 0) | ((uint32_t)a[9*i+1] << 8)) & 0x1FF;
    r->coeffs[8*i+1] = ((a[9*i+1] >> 1) | ((uint32_t)a[9*i+2] << 7)) & 0x1FF;
    r->coeffs[8*i+2] = ((a[9*i+2] >> 2) | ((uint32_t)a[9*i+3] << 6)) & 0x1FF;
    r->coeffs[8*i+3] = ((a[9*i+3] >> 3) | ((uint32_t)a[9*i+4] << 5)) & 0x1FF;
    r->coeffs[8*i+4] = ((a[9*i+4] >> 4) | ((uint32_t)a[9*i+5] << 4)) & 0x1FF;
    r->coeffs[8*i+5] = ((a[9*i+5] >> 5) | ((uint32_t)a[9*i+6] << 3)) & 0x1FF;
    r->coeffs[8*i+6] = ((a[9*i+6] >> 6) | ((uint32_t)a[9*i+7] << 2)) & 0x1FF;
    r->coeffs[8*i+7] = ((a[9*i+7] >> 7) | ((uint32_t)a[9*i+8] << 1)) & 0x1FF;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt0_pack
*
* Description: Bit-pack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*              Input coefficients are assumed to lie in ]Q-2^{D-1}, Q+2^{D-1}].
*
* Arguments:   - unsigned char *r: pointer to output byte array with at least
*                                  POLT0_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt0_pack(unsigned char *r, const poly *a) {
  unsigned int i;
  uint32_t t[4];
  DBENCH_START();

  for(i = 0; i < N/4; ++i) {
    t[0] = Q + (1U << (D-1)) - a->coeffs[4*i+0];
    t[1] = Q + (1U << (D-1)) - a->coeffs[4*i+1];
    t[2] = Q + (1U << (D-1)) - a->coeffs[4*i+2];
    t[3] = Q + (1U << (D-1)) - a->coeffs[4*i+3];

    r[7*i+0]  =  t[0];
    r[7*i+1]  =  t[0] >> 8;
    r[7*i+1] |=  t[1] << 6;
    r[7*i+2]  =  t[1] >> 2;
    r[7*i+3]  =  t[1] >> 10;
    r[7*i+3] |=  t[2] << 4;
    r[7*i+4]  =  t[2] >> 4;
    r[7*i+5]  =  t[2] >> 12;
    r[7*i+5] |=  t[3] << 2;
    r[7*i+6]  =  t[3] >> 6;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyt0_unpack
*
* Description: Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*              Output coefficients lie in ]Q-2^{D-1},Q+2^{D-1}].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const unsigned char *a: byte array with bit-packed polynomial
**************************************************/
void polyt0_unpack(poly *r, const unsigned char *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/4; ++i) {
    r->coeffs[4*i+0]  = a[7*i+0];
    r->coeffs[4*i+0] |= (uint32_t)(a[7*i+1] & 0x3F) << 8;

    r->coeffs[4*i+1]  = a[7*i+1] >> 6;
    r->coeffs[4*i+1] |= (uint32_t)a[7*i+2] << 2;
    r->coeffs[4*i+1] |= (uint32_t)(a[7*i+3] & 0x0F) << 10;

    r->coeffs[4*i+2]  = a[7*i+3] >> 4;
    r->coeffs[4*i+2] |= (uint32_t)a[7*i+4] << 4;
    r->coeffs[4*i+2] |= (uint32_t)(a[7*i+5] & 0x03) << 12;

    r->coeffs[4*i+3]  = a[7*i+5] >> 2;
    r->coeffs[4*i+3] |= (uint32_t)a[7*i+6] << 6;

    r->coeffs[4*i+0] = Q + (1U << (D-1)) - r->coeffs[4*i+0];
    r->coeffs[4*i+1] = Q + (1U << (D-1)) - r->coeffs[4*i+1];
    r->coeffs[4*i+2] = Q + (1U << (D-1)) - r->coeffs[4*i+2];
    r->coeffs[4*i+3] = Q + (1U << (D-1)) - r->coeffs[4*i+3];
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyz_pack
*
* Description: Bit-pack polynomial z with coefficients
*              in [-(GAMMA1 - 1), GAMMA1 - 1].
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - unsigned char *r: pointer to output byte array with at least
*                                  POLZ_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyz_pack(unsigned char *r, const poly *a) {
#if GAMMA1 > (1 << 19)
#error "polyz_pack() assumes GAMMA1 <= 2^{19}"
#endif
  unsigned int i;
  uint32_t t[2];
  DBENCH_START();

  for(i = 0; i < N/2; ++i) {
    /* Map to {0,...,2*GAMMA1 - 2} */
    t[0] = GAMMA1 - 1 - a->coeffs[2*i+0];
    t[0] += ((int32_t)t[0] >> 31) & Q;
    t[1] = GAMMA1 - 1 - a->coeffs[2*i+1];
    t[1] += ((int32_t)t[1] >> 31) & Q;

    r[5*i+0]  = t[0];
    r[5*i+1]  = t[0] >> 8;
    r[5*i+2]  = t[0] >> 16;
    r[5*i+2] |= t[1] << 4;
    r[5*i+3]  = t[1] >> 4;
    r[5*i+4]  = t[1] >> 12;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyz_unpack
*
* Description: Unpack polynomial z with coefficients
*              in [-(GAMMA1 - 1), GAMMA1 - 1].
*              Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const unsigned char *a: byte array with bit-packed polynomial
**************************************************/
void polyz_unpack(poly *r, const unsigned char *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/2; ++i) {
    r->coeffs[2*i+0]  = a[5*i+0];
    r->coeffs[2*i+0] |= (uint32_t)a[5*i+1] << 8;
    r->coeffs[2*i+0] |= (uint32_t)(a[5*i+2] & 0x0F) << 16;

    r->coeffs[2*i+1]  = a[5*i+2] >> 4;
    r->coeffs[2*i+1] |= (uint32_t)a[5*i+3] << 4;
    r->coeffs[2*i+1] |= (uint32_t)a[5*i+4] << 12;

    r->coeffs[2*i+0] = GAMMA1 - 1 - r->coeffs[2*i+0];
    r->coeffs[2*i+0] += ((int32_t)r->coeffs[2*i+0] >> 31) & Q;
    r->coeffs[2*i+1] = GAMMA1 - 1 - r->coeffs[2*i+1];
    r->coeffs[2*i+1] += ((int32_t)r->coeffs[2*i+1] >> 31) & Q;
  }

  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyw1_pack
*
* Description: Bit-pack polynomial w1 with coefficients in [0, 15].
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - unsigned char *r: pointer to output byte array with at least
*                                  POLW1_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyw1_pack(unsigned char *r, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/2; ++i)
    r[i] = a->coeffs[2*i+0] | (a->coeffs[2*i+1] << 4);

  DBENCH_STOP(*tpack);
}
