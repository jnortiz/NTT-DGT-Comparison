#include "poly.h"
#include "dgt.h"
#include "reduce.h"
#include "fips202.h"

// Powers of the n/2-th root of i times the Montgomery constant (2^18 mod p) times the inverse of n/2 mod p.
static uint16_t nthroots[512] = {4075, 0, 1772, 1747, 6053, 9324, 4238, 4249, 3423, 872, 12128, 4983, 3032, 5397, 3287, 52, 3126, 5889, 7351, 6404, 10611, 11847, 2183, 928, 4502, 5100, 9132, 10719, 10386, 10090, 9289, 3041, 8001, 7672, 2853, 2683, 4433, 6075, 427, 2781, 6207, 5738, 10004, 7237, 4969, 6703, 9637, 428, 10319, 7273, 1722, 12288, 10689, 1543, 10597, 8236, 4687, 5157, 8226, 7156, 11034, 4268, 3668, 6523, 4649, 11894, 3738, 6369, 5041, 3208, 2976, 10637, 2322, 7186, 3499, 1379, 674, 10411, 7971, 1333, 12171, 4398, 9351, 10336, 7491, 735, 10367, 8890, 5534, 6920, 5631, 6989, 10713, 1955, 99, 10452, 4365, 1427, 11015, 408, 4054, 7074, 10220, 10939, 243, 10248, 8725, 435, 10031, 9946, 3174, 5094, 282, 7831, 5140, 9636, 8080, 8360, 5780, 4928, 262, 12133, 6170, 11088, 6868, 1182, 12121, 11450, 3429, 8945, 1254, 8671, 872, 10038, 4951, 4196, 10713, 4683, 6903, 11204, 10783, 6215, 2887, 2262, 4960, 1009, 6857, 8485, 9881, 3490, 9230, 11761, 2705, 6529, 4867, 585, 2134, 10420, 918, 2379, 2781, 1033, 10857, 10773, 3248, 8163, 3326, 1465, 4805, 8794, 3437, 7612, 10671, 5833, 8174, 10905, 351, 4350, 4099, 9135, 6093, 6764, 6333, 92, 8390, 3132, 1223, 191, 4111, 8985, 9326, 9792, 1118, 7837, 8773, 1354, 1615, 5719, 7952, 4666, 8460, 4265, 5816, 581, 7467, 6196, 2376, 11176, 3555, 3101, 12228, 4115, 1219, 5557, 9529, 3699, 4699, 5561, 2678, 8715, 6780, 372, 1103, 12203, 6180, 2649, 10635, 6615, 4925, 8039, 4416, 7501, 6862, 8424, 54, 3731, 11681, 8642, 11260, 3250, 9444, 6904, 6741, 11864, 5321, 1143, 2822, 3315, 11911, 3930, 11824, 11360, 3447, 2147, 3368, 8482, 11993, 11993, 1883, 10423, 5307, 5979, 7528, 5620, 9281, 8901, 2361, 2388, 8661, 1743, 338, 1781, 6211, 657, 5459, 11124, 1284, 7452, 3588, 4762, 7064, 1121, 10239, 9683, 4107, 7602, 668, 2105, 6101, 10086, 5431, 3078, 11012, 4475, 1133, 7644, 1899, 7516, 12278, 7264, 3591, 2192, 6421, 6009, 2324, 10598, 3937, 837, 7083, 9962, 9717, 4775, 8840, 7551, 11774, 1440, 4783, 330, 11079, 4685, 2396, 6796, 1898, 9664, 2391, 5390, 8913, 12226, 6246, 2889, 4837, 3165, 10831, 1728, 3007, 6203, 4930, 10229, 9375, 1263, 1572, 6607, 3596, 10658, 2275, 9970, 8808, 9201, 7844, 9867, 2284, 9680, 6970, 6630, 10457, 9490, 1109, 6933, 3849, 9398, 968, 2703, 7495, 4205, 8860, 8359, 10227, 10262, 10467, 3503, 11352, 8393, 8404, 6903, 5677, 7250, 8872, 4543, 8569, 9036, 2087, 4165, 7774, 12142, 3344, 8860, 6154, 10211, 11073, 836, 2795, 8633, 6860, 10851, 6377, 10063, 1936, 11573, 10208, 5389, 4128, 10622, 8134, 885, 1116, 5796, 9430, 7067, 8387, 8304, 1040, 8057, 10994, 5213, 8915, 10415, 5820, 691, 7310, 2108, 930, 3179, 8924, 3494, 10706, 2815, 5053, 7373, 2610, 7791, 6655, 5116, 7495, 11462, 10004, 2359, 4515, 7881, 2107, 4660, 2293, 3183, 2030, 4258, 3035, 9661, 2751, 4121, 6622, 5969, 9235, 3797, 4502, 375, 5654, 6101, 1369, 12194, 3510, 10416, 3810, 3703, 790, 3699, 3227, 6329, 3503, 4021, 10472, 7994, 4085, 5348, 3401, 1034, 6868, 12146, 2888, 922, 11889, 3358, 4445, 7293, 477, 8925, 2877, 1685, 10960, 983, 12122, 108, 2746, 4897, 3763, 3943, 3610, 1630, 4520, 11771, 10333, 7756, 11886, 3665, 10057, 664, 7816, 128, 6438, 2164, 11893, 6955, 11994, 7734};

static uint16_t invnthroots[512] = {1024, 0, 9727, 10418, 10252, 5570, 6729, 1799, 6145, 6498, 4995, 3465, 8415, 5234, 8671, 4424, 286, 8047, 7783, 1234, 7360, 8578, 10724, 8918, 7105, 4273, 1972, 9942, 5143, 2917, 1715, 10607, 6317, 4447, 10862, 5378, 4900, 3445, 6864, 10150, 11813, 3491, 1365, 11745, 9536, 11096, 3616, 8387, 3433, 7428, 6783, 1053, 6591, 10834, 3881, 8723, 4560, 4267, 2088, 1034, 6578, 7183, 2079, 876, 8424, 10631, 11105, 9158, 3254, 10501, 4529, 11417, 6973, 11752, 9811, 2824, 2671, 7807, 9658, 921, 2829, 3379, 212, 12215, 6991, 2390, 2477, 9053, 59, 10039, 4334, 10526, 7165, 7773, 9417, 6788, 3699, 9002, 3929, 10094, 7845, 11574, 6512, 764, 6945, 9328, 4876, 10236, 4439, 4412, 6676, 9473, 6282, 1520, 11686, 10713, 9790, 6905, 8536, 11160, 7579, 9766, 3442, 11270, 9028, 3077, 1432, 456, 4835, 755, 7056, 4482, 8993, 1864, 8676, 5775, 3138, 8030, 8381, 2138, 459, 10144, 2673, 4180, 3902, 10856, 11273, 11625, 4305, 7454, 7073, 3379, 5435, 9597, 3589, 417, 11308, 4076, 11462, 10376, 1274, 2757, 2342, 11320, 5655, 7842, 756, 4958, 711, 10888, 4554, 562, 2378, 1722, 821, 7596, 568, 3149, 9481, 9157, 3079, 3750, 7837, 10974, 8796, 4872, 3024, 9998, 11638, 4167, 3110, 5081, 5595, 4407, 8611, 3365, 8738, 8382, 4614, 12147, 6222, 6494, 4291, 11723, 1095, 8181, 8980, 4641, 7434, 951, 6504, 983, 5385, 322, 7709, 11761, 7902, 5129, 1758, 5228, 6402, 149, 12013, 2619, 7432, 10201, 9561, 7486, 3774, 512, 2198, 12201, 7637, 7269, 4915, 178, 10974, 187, 6764, 3963, 5331, 3192, 535, 3935, 2359, 10191, 8266, 2727, 2867, 3084, 598, 4963, 7944, 8956, 3545, 4361, 1919, 10370, 10690, 1907, 7545, 5699, 7725, 2258, 7984, 6434, 637, 277, 6581, 9628, 8111, 4054, 411, 10908, 3757, 12053, 3010, 7683, 5247, 2592, 1185, 9862, 8622, 3055, 7376, 2909, 1994, 6631, 8027, 1704, 4128, 3788, 6722, 5926, 11795, 5654, 3430, 4350, 6783, 2699, 3622, 9356, 11393, 9361, 10909, 10883, 4268, 3447, 9817, 2035, 8979, 8810, 4193, 543, 9523, 737, 8135, 3786, 8742, 3278, 4783, 4508, 9255, 5244, 11124, 704, 3121, 9548, 9421, 9472, 7873, 9048, 7131, 9817, 3924, 128, 113, 4559, 4987, 11393, 2663, 8359, 3294, 5219, 8003, 9438, 3414, 12180, 1424, 8436, 11323, 4998, 11861, 10598, 8698, 7197, 3689, 4120, 8787, 125, 6122, 6950, 766, 636, 4526, 7306, 10546, 9622, 724, 4589, 2025, 3397, 8905, 1446, 2924, 11830, 8707, 10375, 7505, 4157, 9736, 4989, 1618, 11036, 755, 4835, 3405, 4225, 4709, 10150, 8492, 1224, 7488, 287, 9236, 7082, 4257, 6881, 4454, 940, 5071, 1247, 1268, 4884, 1863, 2217, 3698, 974, 11945, 11664, 3355, 11289, 4540, 10257, 4994, 293, 5238, 607, 2153, 4752, 4472, 10375, 8620, 12219, 11932, 7563, 3395, 6056, 1587, 3187, 7721, 6444, 10098, 6625, 9750, 1649, 4121, 7774, 7542, 8195, 11453, 855, 5562, 7669, 5773, 8477, 1513, 7378, 6671, 1950, 6410, 4018, 4049, 1205, 604, 1600, 10533, 3774, 10209, 4807, 11959, 9223, 48, 8922, 7277, 3752, 4034, 7883, 10059, 5021, 9005, 921, 7223, 3000, 1691, 8207, 3336, 3871, 6395, 1765, 414, 3089, 1500, 3468, 7240, 6968, 1626, 8221, 980, 7183, 4612, 6472, 8927, 5479, 12122, 8756, 12264, 2580, 9793, 10308, 11302, 10357, 6596, 4561, 7300, 4547, 4961, 6800, 7141, 7897, 2167, 11322};

/*************************************************
* Name:        coeff_freeze
* 
* Description: Fully reduces an integer modulo q in constant time
*
* Arguments:   uint16_t x: input integer to be reduced
*              
* Returns integer in {0,...,q-1} congruent to x modulo q
**************************************************/
static uint16_t coeff_freeze(uint16_t x)
{
  uint16_t m,r;
  int16_t c;
  r = x % NEWHOPE_Q;

  m = r - NEWHOPE_Q;
  c = m;
  c >>= 15;
  r = m ^ ((r^m)&c);

  return r;
}

/*************************************************
* Name:        flipabs
* 
* Description: Computes |(x mod q) - Q/2|
*
* Arguments:   uint16_t x: input coefficient
*              
* Returns |(x mod q) - Q/2|
**************************************************/
static uint16_t flipabs(uint16_t x)
{
  int16_t r,m;
  r = coeff_freeze(x);

  r = r - NEWHOPE_Q/2;
  m = r >> 15;
  return (r + m) ^ m;
}

/*************************************************
* Name:        poly_frombytes
* 
* Description: De-serialization of a polynomial
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array
**************************************************/
void poly_frombytes(poly *r, const unsigned char *a)
{
  int i;
  for(i=0;i<NEWHOPE_N/4;i++)
  {
    r->coeffs[4*i+0] =                               a[7*i+0]        | (((uint16_t)a[7*i+1] & 0x3f) << 8);
    r->coeffs[4*i+1] = (a[7*i+1] >> 6) | (((uint16_t)a[7*i+2]) << 2) | (((uint16_t)a[7*i+3] & 0x0f) << 10);
    r->coeffs[4*i+2] = (a[7*i+3] >> 4) | (((uint16_t)a[7*i+4]) << 4) | (((uint16_t)a[7*i+5] & 0x03) << 12);
    r->coeffs[4*i+3] = (a[7*i+5] >> 2) | (((uint16_t)a[7*i+6]) << 6);
  }
}

/*************************************************
* Name:        poly_tobytes
* 
* Description: Serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const poly *p:    pointer to input polynomial
**************************************************/
void poly_tobytes(unsigned char *r, const poly *p)
{
  int i;
  uint16_t t0,t1,t2,t3;
  for(i=0;i<NEWHOPE_N/4;i++)
  {
    t0 = coeff_freeze(p->coeffs[4*i+0]);
    t1 = coeff_freeze(p->coeffs[4*i+1]);
    t2 = coeff_freeze(p->coeffs[4*i+2]);
    t3 = coeff_freeze(p->coeffs[4*i+3]);

    r[7*i+0] =  t0 & 0xff;
    r[7*i+1] = (t0 >> 8) | (t1 << 6);
    r[7*i+2] = (t1 >> 2);
    r[7*i+3] = (t1 >> 10) | (t2 << 4);
    r[7*i+4] = (t2 >> 4);
    r[7*i+5] = (t2 >> 12) | (t3 << 2);
    r[7*i+6] = (t3 >> 6);
  }
}

/*************************************************
* Name:        poly_compress
* 
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const poly *p:    pointer to input polynomial
**************************************************/
void poly_compress(unsigned char *r, const poly *p)
{
  unsigned int i,j,k=0;

  uint32_t t[8];

  for(i=0;i<NEWHOPE_N;i+=8)
  {
    for(j=0;j<8;j++)
    {
      t[j] = coeff_freeze(p->coeffs[i+j]);
      t[j] = (((t[j] << 3) + NEWHOPE_Q/2)/NEWHOPE_Q) & 0x7;
    }

    r[k]   =  t[0]       | (t[1] << 3) | (t[2] << 6);
    r[k+1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[k+2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    k += 3;
  }
}

/*************************************************
* Name:        poly_decompress
* 
* Description: De-serialization and subsequent decompression of a polynomial; 
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array
**************************************************/
void poly_decompress(poly *r, const unsigned char *a)
{
  unsigned int i,j;
  for(i=0;i<NEWHOPE_N;i+=8)
  {
    r->coeffs[i+0] =  a[0] & 7;
    r->coeffs[i+1] = (a[0] >> 3) & 7;
    r->coeffs[i+2] = (a[0] >> 6) | ((a[1] << 2) & 4);
    r->coeffs[i+3] = (a[1] >> 1) & 7;
    r->coeffs[i+4] = (a[1] >> 4) & 7;
    r->coeffs[i+5] = (a[1] >> 7) | ((a[2] << 1) & 6);
    r->coeffs[i+6] = (a[2] >> 2) & 7;
    r->coeffs[i+7] = (a[2] >> 5);
    a += 3;
    for(j=0;j<8;j++)
      r->coeffs[i+j] = ((uint32_t)r->coeffs[i+j] * NEWHOPE_Q + 4) >> 3;
  }
}

/*************************************************
* Name:        poly_frommsg
* 
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r:                  pointer to output polynomial
*              - const unsigned char *msg: pointer to input message
**************************************************/
void poly_frommsg(poly *r, const unsigned char *msg)
{
  unsigned int i,j,mask;
  for(i=0;i<32;i++) // XXX: MACRO for 32
  {
    for(j=0;j<8;j++)
    {
      mask = -((msg[i] >> j)&1);
      r->coeffs[8*i+j+  0] = mask & (NEWHOPE_Q/2);
      r->coeffs[8*i+j+256] = mask & (NEWHOPE_Q/2);
#if (NEWHOPE_N == 1024)
      r->coeffs[8*i+j+512] = mask & (NEWHOPE_Q/2);
      r->coeffs[8*i+j+768] = mask & (NEWHOPE_Q/2);
#endif
    }
  }
}

/*************************************************
* Name:        poly_tomsg
* 
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - unsigned char *msg: pointer to output message
*              - const poly *x:      pointer to input polynomial
**************************************************/
void poly_tomsg(unsigned char *msg, const poly *x)
{
  unsigned int i;
  uint16_t t;

  for(i=0;i<32;i++)
    msg[i] = 0;

  for(i=0;i<256;i++)
  {
    t  = flipabs(x->coeffs[i+  0]);
    t += flipabs(x->coeffs[i+256]);
#if (NEWHOPE_N == 1024)
    t += flipabs(x->coeffs[i+512]);
    t += flipabs(x->coeffs[i+768]);
    t = ((t - NEWHOPE_Q));
#else
    t = ((t - NEWHOPE_Q/2));
#endif

    t >>= 15;
    msg[i>>3] |= t<<(i&7);
  }
}
 
/*************************************************
* Name:        poly_uniform
* 
* Description: Sample a polynomial deterministically from a seed,
*              with output polynomial looking uniformly random
*
* Arguments:   - poly *a:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed
**************************************************/
void poly_uniform(poly *a, const unsigned char *seed)
{
  unsigned int ctr=0;
  uint16_t val;
  uint64_t state[25];
  uint8_t buf[SHAKE128_RATE];
  uint8_t extseed[NEWHOPE_SYMBYTES+1];
  int i,j;

  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    extseed[i] = seed[i];

  for(i=0;i<NEWHOPE_N/64;i++) /* generate a in blocks of 64 coefficients */
  {
    ctr = 0;
    extseed[NEWHOPE_SYMBYTES] = i; /* domain-separate the 16 independent calls */
    shake128_absorb(state, extseed, NEWHOPE_SYMBYTES+1);
    while(ctr < 64) /* Very unlikely to run more than once */
    {
      shake128_squeezeblocks(buf,1,state);
      for(j=0;j<SHAKE128_RATE && ctr < 64;j+=2)
      {
        val = (buf[j] | ((uint16_t) buf[j+1] << 8));
        if(val < 5*NEWHOPE_Q)
        {
          a->coeffs[i*64+ctr] = val;
          ctr++;
        }
      }
    }
  }
}

/*************************************************
* Name:        hw
* 
* Description: Compute the Hamming weight of a byte
*
* Arguments:   - unsigned char a: input byte
**************************************************/
static unsigned char hw(unsigned char a)
{
  unsigned char i, r = 0;
  for(i=0;i<8;i++)
    r += (a >> i) & 1;
  return r;
}

/*************************************************
* Name:        poly_sample
* 
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter k=8
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed 
*              - unsigned char nonce:       one-byte input nonce
**************************************************/
void poly_sample(poly *r, const unsigned char *seed, unsigned char nonce)
{
#if NEWHOPE_K != 8
#error "poly_sample in poly.c only supports k=8"
#endif
  unsigned char buf[128], a, b;
  int i,j;

  unsigned char extseed[NEWHOPE_SYMBYTES+2];

  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[NEWHOPE_SYMBYTES] = nonce;

  for(i=0;i<NEWHOPE_N/64;i++) /* Generate noise in blocks of 64 coefficients */
  {
    extseed[NEWHOPE_SYMBYTES+1] = i;
    shake256(buf,128,extseed,NEWHOPE_SYMBYTES+2);
    for(j=0;j<64;j++)
    {
      a = buf[2*j];
      b = buf[2*j+1];
      r->coeffs[64*i+j] = hw(a) + NEWHOPE_Q - hw(b);
      /*
      t = buf[j] | ((uint32_t)buf[j+1] << 8) | ((uint32_t)buf[j+2] << 16) | ((uint32_t)buf[j+3] << 24);
      d = 0;
      for(k=0;k<8;k++)
        d += (t >> k) & 0x01010101;
      a = d & 0xff;
      b = ((d >>  8) & 0xff);
      c = ((d >> 16) & 0xff);
      d >>= 24;
      r->coeffs[64*i+j/2]   = a + NEWHOPE_Q - b;
      r->coeffs[64*i+j/2+1] = c + NEWHOPE_Q - d;
      */
    }
  }
}

/*************************************************
* Name:        poly_sample_dgt
* 
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter k=8
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed 
*              - unsigned char nonce:       one-byte input nonce
**************************************************/
void poly_sample_dgt(poly *r, const unsigned char *seed, unsigned char nonce)
{
#if NEWHOPE_K != 8
#error "poly_sample in poly.c only supports k=8"
#endif
  unsigned char buf[128], a, b, c, d;
  uint16_t x, y;
  int i, j, k;

  unsigned char extseed[NEWHOPE_SYMBYTES+2];

  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[NEWHOPE_SYMBYTES] = nonce;

  k = 0;
  for(i=0;i<NEWHOPE_N/64;i++) /* Generate noise in blocks of 64 coefficients */
  {
    extseed[NEWHOPE_SYMBYTES+1] = i;
    shake256(buf,128,extseed,NEWHOPE_SYMBYTES+2);
    for(j=0;j<64;j+=2)
    {
      a = buf[2*j];
      b = buf[2*j+1];
      c = buf[2*j+2];
      d = buf[2*j+3];
      x = hw(a) + NEWHOPE_Q - hw(b);
      y = hw(c) + NEWHOPE_Q - hw(d);

      r->coeffs[64*i+j] = montgomery_reduce((uint32_t)x * nthroots[k]) +
                          (NEWHOPE_3Q - montgomery_reduce((uint32_t)y * nthroots[k+1]));

      r->coeffs[64*i+j+1] = montgomery_reduce((uint32_t)x * nthroots[k+1] + (uint32_t)y * nthroots[k]);

      k += 2;
    }
  } 

  dgt((uint16_t *)r->coeffs);

}



/*************************************************
* Name:        poly_mul
* 
* Description: Multiply two polynomials point-wise (i.e., coefficient-wise).
               Inputs are assumed to be NEWHOPE_N/2 Gaussian integers.
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_mul(poly *r, const poly *a, const poly *b)
{
  /* Notice that the output of the function poly_mul_pointwise is only Montgomery reduced.
  But, after every addition, authors perform a reduction using the % operator. */
  int i;
  uint32_t t, s;

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
    t = montgomery_reduce(3186*a->coeffs[i]); // R^2 mod p = 3186
    s = montgomery_reduce(3186*a->coeffs[i+1]);

    r->coeffs[i]   = montgomery_reduce((uint32_t)t * b->coeffs[i]) +
                    (NEWHOPE_3Q - montgomery_reduce((uint32_t)s * b->coeffs[i+1]));

    r->coeffs[i+1] = montgomery_reduce((uint32_t)t * b->coeffs[i+1]) + 
                    montgomery_reduce((uint32_t)s * b->coeffs[i]);
  }  
}

/*************************************************
* Name:        poly_add
* 
* Description: Add two polynomials
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  int i;
  for(i=0;i<NEWHOPE_N;i++)
    r->coeffs[i] = (a->coeffs[i] + b->coeffs[i]) % NEWHOPE_Q;
}

/*************************************************
* Name:        poly_sub
* 
* Description: Subtract two polynomials
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  int i;
  for(i=0;i<NEWHOPE_N;i++)
    r->coeffs[i] = (a->coeffs[i] + NEWHOPE_3Q - b->coeffs[i]) % NEWHOPE_Q;
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

  int i, j;

  uint16_t copy[NEWHOPE_N];

  for(i = 0; i < NEWHOPE_N; i++) 
  {
    copy[i] = r->coeffs[i];
  }  

  j = 0;
  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
      // The nthroots are already multiplied by the Montgomery constant R
      r->coeffs[i]   = montgomery_reduce((uint32_t)copy[j] * nthroots[i]) + 
                       (NEWHOPE_3Q - montgomery_reduce((uint32_t)copy[NEWHOPE_K2+j] * nthroots[i+1]));
      
      r->coeffs[i+1] = montgomery_reduce(
                       (uint32_t)copy[j] * nthroots[i+1] + 
                       (uint32_t)copy[NEWHOPE_K2+j] * nthroots[i]);

      j++;
  } 

  dgt((uint16_t *)r->coeffs);

}

/*************************************************
* Name:        poly_invdgt
* 
* Description: Inverse DGT transform of a polynomial in place
*              Input is assumed to have coefficients in bitreversal order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_invdgt(poly *r)
{

  int i;
  
  idgt((uint16_t *)r->coeffs);

  uint16_t a, b;
  uint32_t c, d;

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
    a = montgomery_reduce((uint32_t)r->coeffs[i] * invnthroots[i]);
    b = montgomery_reduce((uint32_t)r->coeffs[i+1] * invnthroots[i+1]);
    c = (uint32_t)r->coeffs[i] * invnthroots[i+1];
    d = (uint32_t)r->coeffs[i+1] * invnthroots[i];

    r->coeffs[i] = a + (NEWHOPE_3Q - b);
    r->coeffs[i+1] = montgomery_reduce(c + d);
  }

}

