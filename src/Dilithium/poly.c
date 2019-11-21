#include <stdint.h>
#include <inttypes.h>
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

// Powers of the k-th root of i mod p multiplied by the Montgomery constant R = 2^32 mod p
static uint64_t nthroots[N] = {4193792, 0, 2821147, 6334486, 994247, 8057540, 3169394, 1524137, 79814, 4840185, 8319695, 7164506, 2176280, 6562692, 3490149, 925184, 2948712, 354494, 6081549, 7588297, 3804884, 4028859, 3976594, 3529408, 5689750, 2245574, 966761, 684327, 3186524, 3410182, 6497410, 2433681, 396475, 1380241, 8365472, 4949503, 4006839, 4878683, 3800708, 1856811, 5808897, 7926557, 6901987, 283681, 5611478, 2577717, 7155222, 6568037, 2799534, 1353021, 858550, 6963767, 1770127, 333925, 451881, 14659, 374140, 2070761, 1290052, 3835938, 5470202, 4263250, 5738297, 5025759, 827104, 8028621, 5060504, 7314189, 6006066, 1501070, 2804249, 7835688, 4509232, 4066552, 639903, 3190697, 7517584, 5358416, 2803813, 3363468, 4674322, 7878072, 4525645, 4783630, 6731162, 7901980, 1067671, 5941422, 1409000, 5305903, 426198, 2579555, 1951082, 1162083, 243909, 4735201, 642575, 5863640, 3188531, 738588, 1573633, 7200435, 7426617, 256705, 7385172, 6057001, 6948088, 5234731, 8073504, 6503348, 4117964, 280822, 1497098, 2180552, 921056, 6273178, 6970319, 656410, 7972769, 4562757, 584190, 2609120, 134005, 4088090, 1623997, 5525712, 440562, 1717348, 5235201, 5235201, 1407147, 5412078, 2341446, 93998, 1086555, 4680342, 6434883, 2023281, 6986128, 6958064, 3558757, 3514611, 4038104, 4283292, 2820132, 1333837, 319700, 37601, 5934307, 4446217, 6456128, 4898684, 6082842, 6020335, 106375, 5911332, 5285408, 5775022, 3300818, 4126617, 2468922, 163326, 5159652, 5742676, 8016264, 5221773, 1540916, 794901, 2174431, 8205114, 3461907, 6949856, 8160647, 5932252, 4812353, 8027028, 5655507, 5843170, 2177452, 2693988, 1620873, 539894, 1711185, 360158, 3490549, 6777019, 7189685, 7760971, 1592271, 4255382, 8007689, 6474078, 351796, 7553313, 4063845, 3004020, 5112305, 2577386, 8203557, 752744, 7299493, 2343767, 475527, 2832112, 5043280, 7619882, 5877073, 4864212, 4590156, 7988421, 129024, 982333, 8245601, 5949101, 3872067, 8040950, 7183874, 6857856, 8202854, 6830104, 2339528, 4848582, 357275, 6729319, 2488380, 6640925, 4836462, 6387727, 3370736, 7096852, 5680720, 2282575, 4869488, 5161297, 4152944, 1890499, 3970796, 3848893, 6201472, 4633333, 7881116, 5845333, 2588715, 6176088, 6411708, 802608, 592797, 615132, 5105712, 19690, 6909590, 5961879, 273204, 5134176, 7271255, 7618406};
// Inverse of the k-th root of i mod p multiplied by the squared Montgomery constant R^2 = 2^64 mod p and by the inverse of 128 mod p (required by IDGT)
static uint64_t invnthroots[N] = {83956, 0, 7114256, 3095656, 4412440, 7409917, 4244120, 2778078, 8211468, 6336186, 7662380, 3405698, 7280383, 7092444, 8178767, 1542197, 7201128, 523980, 3648874, 6561974, 4995653, 6813781, 647189, 5773613, 4680682, 2474014, 7826409, 5951490, 6589263, 6721339, 3333687, 3464085, 2500129, 3794673, 7397280, 1684449, 8136413, 3378907, 7632722, 1664734, 3433797, 8324543, 6896988, 6729775, 4606578, 631665, 4357132, 4768249, 3822317, 3092298, 651879, 421637, 5191418, 7036486, 3340544, 7408992, 1421017, 8132111, 7717602, 3772893, 4393212, 8213976, 4094032, 247716, 3033122, 5209648, 8337322, 1792623, 6718636, 7309198, 1841230, 2302113, 3102101, 3203363, 589576, 8104807, 6407146, 419557, 3270988, 452593, 3375932, 2368939, 3299098, 5612163, 5673464, 1766077, 771677, 3023147, 5338570, 7248050, 6180745, 5440601, 8337934, 5775501, 4294397, 7314013, 4507018, 4404893, 3211927, 1280433, 8146199, 1829176, 7444778, 982972, 381011, 4841806, 7116609, 1527705, 7445494, 2420269, 37865, 870450, 6341230, 3672994, 7696423, 5584140, 5785224, 6007790, 1525245, 831729, 1680214, 2024274, 1855022, 123396, 4137433, 7580891, 80289, 5290026, 4130825, 4249592, 1070134, 4885123, 2279517, 7010242, 6356266, 798688, 4954280, 471268, 4291902, 6174791, 2467218, 7592368, 4677067, 402833, 550803, 8042246, 7534559, 3741204, 3591647, 7577549, 5492579, 6872573, 3306204, 4633, 5124369, 8108624, 6313790, 6194189, 4833353, 1248638, 3621652, 6660621, 5846060, 3483342, 2252981, 570228, 34575, 6203867, 7616061, 3241253, 4424332, 7068331, 4250339, 7642421, 294586, 4865618, 307408, 2634667, 6562419, 2002022, 1786091, 2713871, 2674850, 2021442, 4640262, 6274462, 2787254, 4477752, 4770524, 6205570, 4064281, 4197289, 5209648, 3033122, 5534660, 5092487, 4712861, 6280051, 7939300, 3559620, 6937789, 2227111, 2603707, 2767755, 4274515, 4425029, 4009163, 3527669, 6354531, 8112506, 2816542, 91550, 6857679, 3503371, 641831, 470260, 4971135, 4949179, 3108601, 6599708, 5476371, 7539726, 4447342, 3593794, 1485992, 7919867, 5836946, 6591811, 3483604, 8343067, 3674353, 3018856, 2305693, 3407165, 4529946, 1217283, 1588109, 3863916, 1131769, 5303573, 7763871, 6227625, 760287, 7928546, 3681719, 5330333, 2394414, 3336779, 1190449, 8044625, 6300182, 8102448, 5724643, 7552588, 2072299, 3939202};

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
* Name:        poly_dgt
* 
* Description: Forward DGT transform of a polynomial in place
*              Input is assumed to have coefficients in normal order
*              Output has coefficients in bitreversed order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_dgt(poly *r)
{
  uint64_t copy[N];
  unsigned int i, j;
  uint64_t a, b, c, d;
  DBENCH_START();

  j = 0;
  for(i = 0; i < 128; ++i)
  {
    copy[j] = r->coeffs[i];
    copy[j+1] = r->coeffs[i+128];
    
    j = j+2;
  }

  j = 0;
  for(i = 0; i < N; i += 2) 
  {
    a = copy[i] * nthroots[i];
    b = copy[i+1] * nthroots[i+1];
    c = copy[i] * nthroots[i+1];
    d = copy[i+1] * nthroots[i];

    // The nthroots are already multiplied by the Montgomery constant R
    r->coeffs[i] = 2*Q + reduce(a - b);
    r->coeffs[i+1] = montgomery_reduce(c + d);

    j++;
  } 

  dgt(r->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_invdgt_montgomery
*
* Description: Inverse DGT and multiplication with 2^{32}.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invdgt_montgomery(poly *r)
{
  uint64_t copy[N];
  uint64_t a, b, c, d;
  unsigned int i, j;
  DBENCH_START();

  idgt((uint32_t *)r->coeffs);

  for(i = 0; i < N; ++i)
  {
    copy[i] = r->coeffs[i];
  }

  j = 0;
  for(i = 0; i < N; i += 2) 
  {
    a = copy[i] * invnthroots[i];
    b = copy[i+1] * invnthroots[i+1];
    c = copy[i] * invnthroots[i+1];
    d = copy[i+1] * invnthroots[i];

    r->coeffs[j] = 2*Q + reduce(a - b);
    r->coeffs[j+128] = montgomery_reduce(c + d);

    j++;
  }

  for(i = 0; i < N; ++i) 
  {
    r->coeffs[i] = reduce32(r->coeffs[i]);
  }


  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_mul_invmontgomery
*
* Description: Pointwise multiplication of polynomials in DGT domain (Z[i])
*              representation. The output is multiplied by R^-1.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_invmontgomery(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  DBENCH_START();
  
  for(i = 0; i < N; i += 2) 
  {
    r->coeffs[i] = 2*Q + reduce(a->coeffs[i] * (uint64_t)b->coeffs[i] -
                                a->coeffs[i+1] * (uint64_t)b->coeffs[i+1]);

    r->coeffs[i+1] = montgomery_reduce(a->coeffs[i] * (uint64_t)b->coeffs[i+1] + 
                                       a->coeffs[i+1] * (uint64_t)b->coeffs[i]);
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
