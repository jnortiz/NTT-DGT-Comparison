#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "dgt.h"
#include "reduce.h"
#include "rounding.h"
#include "symmetric.h"

#ifdef DBENCH
#include "test/cpucycles.h"
extern const uint64_t timing_overhead;
extern uint64_t *tred, *tadd, *tmul, *tround, *tsample, *tpack;
#define DBENCH_START() uint64_t time = cpucycles()
#define DBENCH_STOP(t) t += cpucycles() - time - timing_overhead
#else
#define DBENCH_START()
#define DBENCH_STOP(t)
#endif

static const uint32_t nthroots[N] = {4193792, 0, 3823250, 6767456, 1126796, 5978524, 5964771, 4724345, 1104991, 2167712, 1445227, 194787, 4580618, 1667385, 2084137, 7475135, 875594, 250783, 1136551, 826930, 532253, 7292158, 771785, 2302759, 7588774, 17198, 7256068, 6388578, 7410027, 7302180, 684253, 8144978, 7866503, 1622053, 3811129, 768456, 3702788, 8224351, 645556, 4594305, 8240644, 6991373, 7166650, 6349893, 7548664, 1021874, 3652664, 8294587, 1783253, 5608777, 7858800, 7879041, 7211335, 152181, 8300364, 2619427, 7465698, 7520184, 1064429, 2707109, 2694524, 5568852, 6960384, 6758456, 8331273, 114672, 2397251, 2265945, 5727714, 399609, 5435961, 7162390, 176926, 3738675, 574460, 1912415, 8336461, 1981530, 7223002, 3724564, 4562262, 7558189, 7503303, 3373592, 5622715, 1949710, 1013159, 4600915, 3164108, 1615600, 7123029, 6646571, 1378334, 7231633, 6992518, 2059166, 5647943, 2253225, 5276943, 4145404, 3364814, 2226731, 1514203, 7667546, 6625043, 4861972, 4149951, 439470, 3462045, 1196038, 4291775, 3090293, 3765172, 2198628, 1432006, 7724871, 6595650, 1064431, 1102567, 1259115, 2180517, 7801920, 7521607, 1294529, 4311190, 6054641, 5267977, 7614849, 5235201, 5235201, 1269859, 963597, 693071, 3813872, 114329, 648335, 4861884, 5110559, 4756922, 3438559, 3468953, 7234618, 6744134, 443911, 3215684, 7301373, 3117505, 3354553, 1457068, 6787987, 8008237, 101599, 4755378, 1390102, 896374, 1339499, 5508450, 3330969, 2984301, 3123612, 6127192, 2732474, 6678617, 3276128, 3659093, 5681605, 4367215, 938293, 7149865, 3592858, 4181821, 3773413, 6425985, 4701912, 7003920, 6414426, 3538150, 3547692, 7911659, 7342838, 6131087, 3252931, 3140757, 6140379, 2587421, 4606183, 2018743, 6406188, 696682, 4931251, 5674234, 311887, 114672, 8331273, 7770379, 8113168, 4671999, 2854068, 5844606, 6099029, 8317830, 4823822, 6777943, 4660411, 6243844, 6103680, 3953949, 7947873, 6577742, 3851540, 2322336, 1763904, 6973231, 3017475, 2891240, 5922646, 6382167, 4264780, 5199263, 5351588, 6087160, 7794919, 8211358, 304839, 6758364, 513914, 5115176, 4175625, 100025, 1757007, 126370, 2235128, 2953390, 1530397, 5683612, 8283087, 4234689, 1818645, 4408944, 7179072, 2984143, 2368969, 793726, 2865610, 7918491, 7954817, 1997910, 1137555, 699400, 5788524, 8166008, 4787826, 7166972, 359787, 2068043, 2027477};

/*************************************************
* Name:        poly_reduce
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in [0,2*Q[.
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
* Description: For all coefficients of in/out polynomial subtract Q if
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
* Description: Inplace reduction of all coefficients of polynomial to
*              standard representatives.
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
* Description: Inplace forward dgt. Output coefficients can be up to
*              16*Q larger than input coefficients.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_dgt(poly *a) {
  DBENCH_START();
  unsigned int i, j;
  uint32_t copy[N];

  for(i = 0; i < N; ++i)
    copy[i] = a->coeffs[i];

  for(i = 0, j = 0; i < N; i += 2, ++j) {             
      a->coeffs[i] = montgomery_reduce((uint64_t)copy[j] * nthroots[i]) + (2*Q - 
                     montgomery_reduce((uint64_t)copy[j+N/2] * nthroots[i+1]));
      a->coeffs[i+1] = montgomery_reduce((uint64_t)copy[j] * nthroots[i+1] + (uint64_t)copy[j+N/2] * nthroots[i]);
  }

  dgt(a->coeffs);
  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_invdgt_tomont
*
* Description: Inplace inverse dgt and multiplication by 2^{32}.
*              Input coefficients need to be less than 2*Q.
*              Output coefficients are less than 2*Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invdgt_tomont(poly *a) {
  DBENCH_START(); 
  invdgt_tomont(a->coeffs);
  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in dgt domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}. Output coefficients are less than 2*Q if input
*              coefficient are less than 22*Q.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
  DBENCH_START();  
  unsigned int i;
  uint32_t aux_b;

  for(i = 0; i < N; i += 2) {
    aux_b = b->coeffs[i];
    c->coeffs[i] = montgomery_reduce((uint64_t)a->coeffs[i] * b->coeffs[i]) + (2*Q - 
                   montgomery_reduce((uint64_t)a->coeffs[i+1] * b->coeffs[i+1]));
    c->coeffs[i+1] = montgomery_reduce((uint64_t)a->coeffs[i] * b->coeffs[i+1] + (uint64_t)a->coeffs[i+1] * aux_b);
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
*              - poly *a0: pointer to output polynomial with coefficients Q + c0
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
*              - poly *a0: pointer to output polynomial with coefficients Q + c0
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
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/
void poly_use_hint(poly *b, const poly *a, const poly *h) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    b->coeffs[i] = use_hint(a->coeffs[i], h->coeffs[i]);

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
  uint32_t t;
  DBENCH_START();

  /* It is ok to leak which coefficient violates the bound since
     the probability for each coefficient is independent of secret
     data but we must not leak the sign of the centralized representative. */
  for(i = 0; i < N; ++i) {
    /* Absolute value of centralized representative */
    t = (Q-1)/2 - a->coeffs[i];
    t ^= (int32_t)t >> 31;
    t = (Q-1)/2 - t;

    if(t >= B) {
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
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_uniform(uint32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
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
*              output stream of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_NBLOCKS ((768+STREAM128_BLOCKBYTES-1)/STREAM128_BLOCKBYTES)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES + 2];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

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
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_eta(uint32_t *a,
                            unsigned int len,
                            const uint8_t *buf,
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
*              output stream from SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_ETA_NBLOCKS ((192 + STREAM128_BLOCKBYTES - 1) \
                                  /STREAM128_BLOCKBYTES)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce)
{
  unsigned int ctr;
  unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

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
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_gamma1m1(uint32_t *a,
                                 unsigned int len,
                                 const uint8_t *buf,
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
*              sampling on output stream of SHAKE256(seed|nonce)
*              or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 16-bit nonce
**************************************************/
#define POLY_UNIFORM_GAMMA1M1_NBLOCKS ((640 + STREAM256_BLOCKBYTES - 1) \
                                       /STREAM256_BLOCKBYTES)
void poly_uniform_gamma1m1(poly *a,
                           const uint8_t seed[CRHBYTES],
                           uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_GAMMA1M1_NBLOCKS*STREAM256_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_GAMMA1M1_NBLOCKS*STREAM256_BLOCKBYTES + 4];
  stream256_state state;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, POLY_UNIFORM_GAMMA1M1_NBLOCKS, &state);

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
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLETA_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyeta_pack(uint8_t *r, const poly *a) {
#if ETA > 7
#error "polyeta_pack() assumes ETA <= 7"
#endif
  unsigned int i;
  uint8_t t[8];
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
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyeta_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

#if ETA <= 3
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
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLT1_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt1_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/8; ++i) {
    r[9*i+0] = (a->coeffs[8*i+0] >> 0);
    r[9*i+1] = (a->coeffs[8*i+0] >> 8) | (a->coeffs[8*i+1] << 1);
    r[9*i+2] = (a->coeffs[8*i+1] >> 7) | (a->coeffs[8*i+2] << 2);
    r[9*i+3] = (a->coeffs[8*i+2] >> 6) | (a->coeffs[8*i+3] << 3);
    r[9*i+4] = (a->coeffs[8*i+3] >> 5) | (a->coeffs[8*i+4] << 4);
    r[9*i+5] = (a->coeffs[8*i+4] >> 4) | (a->coeffs[8*i+5] << 5);
    r[9*i+6] = (a->coeffs[8*i+5] >> 3) | (a->coeffs[8*i+6] << 6);
    r[9*i+7] = (a->coeffs[8*i+6] >> 2) | (a->coeffs[8*i+7] << 7);
    r[9*i+8] = (a->coeffs[8*i+7] >> 1);
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
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyt1_unpack(poly *r, const uint8_t *a) {
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
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLT0_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt0_pack(uint8_t *r, const poly *a) {
#if D != 14
#error "polyt0_pack() assumes D == 14"
#endif
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
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyt0_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/4; ++i) {
    r->coeffs[4*i+0]  = a[7*i+0];
    r->coeffs[4*i+0] |= (uint32_t)a[7*i+1] << 8;
    r->coeffs[4*i+0] &= 0x3FFF;

    r->coeffs[4*i+1]  = a[7*i+1] >> 6;
    r->coeffs[4*i+1] |= (uint32_t)a[7*i+2] << 2;
    r->coeffs[4*i+1] |= (uint32_t)a[7*i+3] << 10;
    r->coeffs[4*i+1] &= 0x3FFF;

    r->coeffs[4*i+2]  = a[7*i+3] >> 4;
    r->coeffs[4*i+2] |= (uint32_t)a[7*i+4] << 4;
    r->coeffs[4*i+2] |= (uint32_t)a[7*i+5] << 12;
    r->coeffs[4*i+2] &= 0x3FFF;

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
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLZ_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyz_pack(uint8_t *r, const poly *a) {
#if GAMMA1 > (1 << 19)
#error "polyz_pack() assumes GAMMA1 - 1 fits in 19 bits"
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
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyz_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/2; ++i) {
    r->coeffs[2*i+0]  = a[5*i+0];
    r->coeffs[2*i+0] |= (uint32_t)a[5*i+1] << 8;
    r->coeffs[2*i+0] |= (uint32_t)a[5*i+2] << 16;
    r->coeffs[2*i+0] &= 0xFFFFF;

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
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLW1_SIZE_PACKED bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyw1_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N/2; ++i)
    r[i] = a->coeffs[2*i+0] | (a->coeffs[2*i+1] << 4);

  DBENCH_STOP(*tpack);
}
