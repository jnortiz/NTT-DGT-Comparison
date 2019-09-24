#include "poly.h"
#include "ntt.h"
#include "dgt.h"
#include "reduce.h"
#include "fips202.h"

static uint16_t nthroots[512] = {4075, 0, 8303, 5115, 2563, 11993, 3525, 6646, 8596, 6606, 12142, 2772, 6969, 7276, 1665, 2110, 2086, 599, 8981, 1176, 4905, 5943, 10134, 7904, 2802, 3650, 1191, 5164, 3282, 141, 6272, 11337, 8127, 4980, 9316, 6310, 3501, 8326, 1381, 8280, 11685, 299, 4344, 7517, 3131, 11574, 3622, 6077, 10687, 4341, 5503, 88, 10496, 11743, 11812, 7409, 3130, 5004, 1167, 3621, 4998, 6771, 12131, 12238, 11823, 4544, 9321, 4895, 10987, 9837, 2940, 1886, 11102, 11538, 10026, 7375, 4983, 214, 8449, 7119, 339, 6178, 7719, 2169, 1627, 1735, 5498, 4929, 700, 2882, 8801, 5050, 1108, 4871, 7962, 10420, 7444, 1415, 6793, 8397, 6947, 684, 9213, 6531, 9018, 8837, 9345, 7374, 4809, 9119, 8672, 6815, 991, 9391, 1851, 6455, 2373, 9761, 466, 1148, 4403, 6374, 2204, 12172, 2216, 4621, 1815, 8346, 1455, 1892, 10044, 12138, 9201, 1027, 3565, 1009, 4972, 11000, 2599, 8313, 1969, 5526, 6943, 6819, 8012, 5139, 5022, 11939, 2246, 12014, 3287, 7657, 7026, 3717, 4050, 6082, 6408, 6728, 10082, 7922, 10843, 7192, 276, 2664, 8265, 8374, 1308, 4825, 3195, 2260, 7096, 3956, 123, 8608, 6587, 7980, 4614, 10565, 4801, 12104, 2243, 4624, 6619, 4161, 3354, 11699, 3253, 3807, 3475, 10275, 5869, 5394, 8888, 9452, 7509, 608, 9974, 10043, 10527, 8100, 4560, 1536, 9821, 7807, 11502, 1624, 4618, 3066, 2316, 2885, 6731, 7941, 2007, 7539, 1867, 546, 8927, 7491, 5261, 11384, 3119, 1744, 6660, 1805, 1820, 7333, 9703, 2476, 2953, 843, 4051, 9975, 854, 9396, 11439, 2026, 7818, 2075, 2963, 10935, 3165, 7227, 5942, 6774, 7392, 10350, 1289, 10177, 9807, 7746, 597, 2800, 2922, 6168, 11740, 8301, 4620, 7669, 2049, 11434, 2917, 2398, 11667, 6916, 2805, 8525, 3730, 6114, 9591, 3409, 11216, 5782, 12242, 9427, 3946, 2837, 10218, 7163, 6156, 11344, 6757, 690, 8562, 9812, 1423, 6632, 9198, 9889, 8723, 3323, 3150, 1206, 3276, 419, 760, 1353, 9983, 300, 9738, 11016, 5212, 4928, 9428, 8287, 5156, 2546, 5042, 5999, 6531, 5681, 1119, 11112, 416, 1235, 10736, 1591, 8156, 834, 1482, 3710, 130, 10279, 12031, 4572, 6510, 4803, 7115, 8591, 336, 6963, 5087, 6509, 8953, 5887, 8030, 3845, 11882, 6635, 9205, 11425, 10793, 10406, 5383, 11855, 9233, 11340, 1213, 2020, 10111, 4598, 11537, 3797, 861, 11470, 8476, 12087, 11803, 5309, 6782, 7846, 9688, 4635, 9923, 570, 667, 5188, 6370, 1996, 11288, 11936, 9251, 7572, 6685, 6823, 5313, 7453, 1784, 11719, 1959, 4169, 4027, 5290, 6152, 7314, 8665, 1370, 5818, 11964, 91, 2102, 3316, 3601, 9634, 6412, 8647, 5287, 10533, 1590, 7430, 5348, 8411, 6966, 8207, 11280, 11719, 4470, 1717, 1411, 7611, 3, 8592, 2349, 6283, 4434, 11983, 793, 10646, 8014, 11551, 2059, 2670, 4813, 9019, 4108, 3772, 6630, 12123, 1039, 4003, 3513, 2839, 177, 8554, 6958, 1901, 2957, 8283, 4343, 9583, 8492, 510, 37, 10266, 10595, 812, 2055, 11744, 2251, 663, 3821, 9474, 9833, 9998, 7331, 5331, 7155, 7047, 3568, 1243, 7578, 768, 8599, 10310, 4302, 11159, 434, 5446, 3709, 4927, 3226, 9702, 2664, 9823, 6831, 5979, 3504, 9087, 9750, 4772, 7614, 3245, 11181, 780, 2090, 4364, 3790, 5558, 10501, 12224, 4226, 9536, 833, 5818, 9442, 8700, 9385, 4514, 3595, 8144, 4155, 2669, 9189, 6953, 10710, 2892, 91, 11662, 4488, 11512, 8264, 10945, 2042};
static uint16_t invnthroots[512] = {1024, 0, 296, 9222, 8865, 11860, 5778, 6771, 7921, 3637, 2058, 1941, 1332, 5222, 9473, 9953, 11775, 7759, 4213, 12063, 1477, 8906, 9172, 3035, 6065, 9169, 12090, 8715, 2415, 559, 10281, 573, 4028, 8292, 3198, 7854, 11271, 6061, 3854, 4345, 3915, 4522, 7307, 11003, 4909, 3005, 6303, 3339, 3746, 7205, 2417, 3320, 5074, 12286, 4926, 10508, 782, 6453, 652, 10108, 4493, 633, 7287, 59, 927, 7246, 2553, 10707, 11961, 2109, 7578, 1208, 10513, 12191, 10210, 5291, 449, 4336, 5532, 5225, 10108, 5055, 3793, 1093, 3422, 7809, 11573, 4321, 1274, 9010, 11729, 2797, 2467, 5270, 11769, 1443, 8576, 7159, 11092, 9890, 8370, 6648, 10138, 6879, 12145, 8947, 6006, 8682, 6642, 9507, 11565, 688, 9724, 10480, 1365, 259, 9703, 1735, 4293, 9519, 11738, 7739, 11487, 11700, 9705, 4368, 3311, 8906, 7974, 10383, 5309, 360, 4149, 8961, 8801, 8009, 2782, 11898, 10926, 9244, 4299, 1366, 5214, 1644, 4655, 1108, 2504, 10824, 9045, 7438, 9507, 9322, 11011, 10331, 4351, 6022, 3237, 1250, 9696, 1311, 2445, 4461, 2079, 771, 498, 6057, 1352, 9068, 8685, 780, 8543, 315, 4361, 1926, 4605, 11725, 1034, 5042, 12064, 4481, 71, 11918, 7082, 10685, 9868, 3839, 5458, 9717, 2947, 5255, 1746, 12194, 10457, 6240, 6255, 9691, 9124, 10529, 9655, 11479, 2165, 7679, 7340, 4556, 9959, 6263, 6984, 8525, 682, 1708, 7761, 10140, 9236, 4396, 11948, 442, 10178, 12202, 8790, 11902, 4466, 9780, 3557, 3732, 253, 878, 4599, 11389, 1178, 6859, 8295, 5439, 3747, 4822, 8493, 552, 268, 11193, 11292, 5073, 2197, 10033, 5111, 9941, 8414, 5675, 1464, 6994, 8626, 11750, 12124, 7011, 7786, 4837, 4173, 40, 558, 558, 7089, 10515, 11161, 5077, 779, 4078, 9151, 3754, 3064, 427, 7049, 10724, 6651, 2569, 9485, 4452, 3547, 7045, 11001, 6594, 1064, 8356, 3685, 4125, 471, 10113, 8692, 6565, 4042, 11051, 4397, 1337, 11672, 166, 2311, 2244, 6573, 6748, 9102, 10670, 10659, 3593, 6798, 10313, 12080, 3574, 8988, 567, 300, 462, 8071, 11380, 6223, 4426, 6, 9967, 4448, 1447, 9496, 11770, 7683, 4051, 997, 8798, 11446, 11354, 10649, 7043, 1599, 8676, 3742, 1235, 9185, 10487, 11539, 9352, 8880, 9246, 9018, 270, 10208, 8951, 4642, 5904, 6736, 8805, 2121, 5892, 1891, 1339, 3585, 3472, 7307, 959, 11165, 4326, 703, 4665, 8859, 359, 3000, 10065, 5919, 5445, 1134, 10308, 911, 9496, 4511, 7565, 11397, 3617, 4491, 1461, 5110, 8489, 6513, 1862, 427, 5165, 724, 11363, 12149, 11533, 7248, 2841, 7496, 8395, 4929, 1097, 11683, 8056, 5616, 7480, 1273, 2431, 6341, 10079, 10743, 3303, 9674, 2825, 3925, 10701, 4683, 10719, 4692, 9630, 2429, 6156, 5939, 2749, 6026, 12109, 4035, 1653, 2481, 6550, 5814, 931, 3689, 1217, 11972, 4028, 3380, 4622, 9132, 9022, 9188, 5835, 2743, 4362, 6489, 1842, 10681, 3983, 2380, 15, 2017, 5693, 2381, 1977, 11470, 4469, 7784, 5941, 7095, 11238, 10820, 5004, 3090, 2210, 2448, 4705, 6795, 6413, 10527, 6860, 5588, 2772, 749, 1682, 1630, 12248, 8065, 6075, 545, 9127, 3240, 1810, 9742, 2820, 7854, 11888, 10226, 7875, 8097, 4843, 5889, 8291, 4345, 4764, 6740, 9137, 8829, 6120, 5521, 10068, 10197, 8012, 9135, 11606, 1567, 7161, 9672, 1949, 4997, 973, 8115, 1816, 9321, 6186, 7133, 2709, 2123, 5233, 2426, 7071, 506, 9443, 1919, 134, 260, 5296};

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
//  uint32_t t, d, a, b, c;
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
* Name:        poly_pointwise
* 
* Description: Multiply two polynomials pointwise (i.e., coefficient-wise).
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_mul_pointwise(poly *r, const poly *a, const poly *b)
{
  int i;
  uint16_t t;
  for(i=0;i<NEWHOPE_N;i++)
  {
    t            = montgomery_reduce(3186*b->coeffs[i]); /* t is now in Montgomery domain */
    r->coeffs[i] = montgomery_reduce(a->coeffs[i] * t);  /* r->coeffs[i] is back in normal domain */
  }
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
  uint16_t copy[NEWHOPE_N];

  for(i = 0; i < NEWHOPE_N; i++) 
  {
    copy[i] = montgomery_reduce(3186*a->coeffs[i]);
  }

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {             
    r->coeffs[i]   = montgomery_reduce((uint32_t)copy[i] * b->coeffs[i]) +
                     (NEWHOPE_3Q - montgomery_reduce((uint32_t)copy[i+1] * b->coeffs[i+1]));

    r->coeffs[i+1] = montgomery_reduce((uint32_t)copy[i] * b->coeffs[i+1]) + 
                      montgomery_reduce((uint32_t)copy[i+1] * b->coeffs[i]);
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
* Name:        poly_ntt
* 
* Description: Forward NTT transform of a polynomial in place
*              Input is assumed to have coefficients in bitreversed order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  mul_coefficients(r->coeffs, gammas_bitrev_montgomery);
  ntt((uint16_t *)r->coeffs, gammas_bitrev_montgomery);
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
      
      r->coeffs[i+1] = montgomery_reduce((uint32_t)copy[j] * nthroots[i+1]) + 
                       montgomery_reduce((uint32_t)copy[NEWHOPE_K2+j] * nthroots[i]);

      j++;
  } 

  dgt((uint16_t *)r->coeffs);

}

/*************************************************
* Name:        poly_invntt
* 
* Description: Inverse NTT transform of a polynomial in place
*              Input is assumed to have coefficients in normal order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_invntt(poly *r)
{
  bitrev_vector(r->coeffs);
  ntt((uint16_t *)r->coeffs, omegas_inv_bitrev_montgomery);
  mul_coefficients(r->coeffs, gammas_inv_montgomery);
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

  int i, j;
  
  idgt((uint16_t *)r->coeffs);

  uint16_t copy[NEWHOPE_N];

  for(i = 0; i < NEWHOPE_N; i++) 
  {
    copy[i] = r->coeffs[i];
  }

  j = 0;
  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
    r->coeffs[j] = 
      (montgomery_reduce((uint32_t)copy[i] * invnthroots[i]) + 
      (NEWHOPE_3Q - montgomery_reduce((uint32_t)copy[i+1] * invnthroots[i+1])));

    r->coeffs[j+NEWHOPE_K2] = 
      (montgomery_reduce((uint32_t)copy[i] * invnthroots[i+1]) + 
      montgomery_reduce((uint32_t)copy[i+1] * invnthroots[i]));

    j++;
  }
}
