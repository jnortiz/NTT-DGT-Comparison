#include "poly.h"
#include "ntt.h"
#include "dgt.h"
#include "reduce.h"
#include "fips202.h"

static uint16_t nthroots[512] = {4075, 0, 1673, 6583, 7664, 4766, 8276, 6738, 5040, 4260, 2998, 2427, 5929, 4196, 131, 2902, 7133, 2546, 9087, 6166, 8601, 2205, 3220, 7918, 1441, 7198, 1536, 7581, 9735, 3531, 1646, 8696, 10874, 4845, 8842, 6760, 5077, 6070, 8114, 8631, 2986, 3533, 10594, 2963, 1909, 9181, 11643, 4320, 6966, 8411, 11302, 219, 9802, 2057, 1425, 3959, 8508, 12106, 2712, 189, 3670, 3744, 4590, 2119, 3401, 9452, 12249, 596, 9627, 6311, 10805, 664, 2528, 2625, 8646, 8506, 5068, 6817, 10061, 1215, 4514, 8694, 7697, 7496, 8038, 1169, 3227, 10392, 10653, 4766, 7172, 6265, 3457, 3561, 923, 4239, 3323, 8723, 9046, 228, 5740, 3751, 7054, 10098, 6722, 7628, 6758, 9226, 9707, 7594, 7381, 4589, 11950, 6178, 7640, 5493, 8108, 4175, 12097, 6516, 9588, 2175, 11608, 64, 8061, 12147, 10095, 10482, 8665, 10919, 10837, 7962, 7359, 2235, 7725, 3252, 3479, 2765, 5110, 9566, 387, 6356, 10241, 10982, 10565, 4614, 964, 8473, 10330, 7165, 8559, 8049, 12194, 6905, 11801, 5486, 4909, 6467, 3815, 9177, 1108, 9044, 11062, 5702, 11042, 7310, 4187, 7965, 11156, 10806, 4509, 3161, 4673, 11296, 9156, 8031, 12242, 2862, 9133, 1853, 11869, 7306, 2039, 4049, 418, 11437, 10384, 1963, 11264, 1480, 8075, 1108, 4544, 11823, 8695, 3437, 50, 6338, 9417, 599, 5489, 7642, 6196, 8606, 1865, 2674, 1642, 4946, 353, 1001, 1245, 1515, 5799, 11997, 11396, 7487, 6658, 11770, 4555, 4080, 3346, 2414, 3237, 4783, 10843, 5097, 5892, 1364, 9531, 1191, 5261, 7837, 11475, 7793, 1760, 9671, 757, 6687, 6366, 461, 1130, 434, 6135, 2839, 535, 4086, 8236, 6999, 5573, 8232, 1559, 3622, 9026, 833, 7167, 7601, 4620, 7669, 6254, 4285, 10172, 11612, 9724, 12023, 4437, 12234, 11805, 10752, 10876, 8937, 7540, 10515, 10687, 7948, 209, 5072, 11633, 11838, 10396, 10257, 2421, 10568, 7969, 9643, 4244, 7125, 1517, 4751, 11428, 11470, 2203, 6943, 9215, 779, 11973, 8623, 696, 590, 11118, 486, 8502, 4294, 10815, 10591, 7150, 4277, 6563, 1503, 8635, 11349, 7627, 8844, 1296, 5768, 8476, 3367, 8873, 3944, 7169, 7467, 3821, 663, 11698, 6662, 3172, 2287, 2463, 4878, 9265, 10484, 9087, 8813, 10685, 12176, 480, 164, 9124, 7227, 3986, 5894, 8448, 2903, 8444, 584, 5735, 3760, 6383, 9783, 11681, 8200, 8084, 10002, 7309, 4162, 2387, 3317, 5317, 2932, 9087, 7462, 12198, 8099, 3963, 9312, 11382, 8280, 2262, 8626, 6635, 11882, 10773, 1818, 8980, 12248, 9899, 473, 4847, 10112, 9961, 11475, 10878, 4542, 11613, 10662, 10834, 1892, 1877, 5863, 6821, 584, 1483, 9374, 2646, 11766, 7652, 11973, 1342, 4596, 11501, 9269, 8554, 5331, 114, 5965, 7772, 5685, 4007, 1382, 7570, 1754, 3661, 1209, 6474, 11296, 601, 7783, 7333, 1820, 4018, 5609, 7281, 6354, 4323, 6726, 7378, 3599, 2589, 6207, 1835, 8896, 9956, 157, 10203, 599, 302, 3583, 328, 8895, 2635, 2505, 6429, 68, 7409, 7205, 11297, 11284, 7501, 10301, 130, 2010, 5582, 417, 8328, 2334, 6063, 3317, 9670, 1491, 3407, 1903, 2447, 9904, 8027, 3384, 9391, 991, 10593, 1437, 7872, 1855, 11336, 906, 39, 5087, 215, 2251, 6078, 2140, 11279, 4678, 4275, 1643, 7759, 9764, 2174, 8865, 9904, 7429, 8606, 9351, 6084, 6466, 840, 1521, 10801, 6122, 2316, 9404, 8983, 1028, 10803, 10570, 2571, 12123, 422, 10192, 129, 1012, 11645, 10057, 9234, 5860};
static uint16_t invnthroots[512] = {1024, 0, 1367, 828, 8824, 5955, 580, 6192, 2344, 7967, 7968, 518, 8778, 2406, 12101, 1069, 3301, 567, 1080, 2310, 726, 3453, 9146, 9385, 5845, 7551, 12078, 8410, 4595, 6040, 10599, 3762, 7159, 8576, 8947, 676, 7881, 9097, 2553, 10320, 1604, 1872, 5668, 3412, 9272, 9186, 4758, 4615, 1588, 8364, 9614, 4337, 3879, 6855, 6968, 3779, 2166, 9467, 541, 8377, 10858, 6496, 4562, 9867, 1832, 6240, 9401, 3667, 11373, 1540, 10541, 11540, 9025, 1367, 2650, 3590, 3155, 3455, 62, 2207, 8115, 10473, 4753, 10906, 3107, 2057, 9289, 1382, 11583, 10052, 8955, 10880, 2233, 5396, 1126, 8529, 10952, 7892, 7375, 4270, 10797, 3527, 3413, 3682, 1831, 6979, 7398, 8001, 9767, 4386, 8616, 5472, 2181, 5055, 9781, 11332, 594, 2971, 2879, 10915, 526, 4118, 4741, 9739, 8835, 7894, 1223, 4073, 7496, 3894, 4362, 4419, 3186, 6006, 2205, 11146, 6184, 11454, 1874, 8170, 1968, 925, 11048, 966, 7247, 11255, 3778, 10264, 8097, 5620, 7717, 5889, 4496, 7921, 10494, 6061, 6732, 9436, 541, 3975, 9137, 6740, 11464, 7073, 11937, 7683, 9687, 11448, 3855, 4922, 8835, 12064, 8124, 12256, 12024, 6993, 9485, 7837, 4417, 10751, 5424, 9031, 7091, 6061, 617, 2316, 11636, 7623, 825, 4788, 12027, 8499, 5043, 11362, 10254, 20, 7312, 8078, 10430, 1311, 5783, 763, 5603, 9715, 8253, 8943, 1590, 7799, 3617, 11397, 7770, 2982, 2801, 2559, 1250, 5237, 8547, 8830, 3922, 9410, 11764, 12205, 10828, 7432, 2445, 7828, 5443, 11371, 2092, 7088, 4118, 1553, 8874, 5607, 11513, 7448, 9359, 5379, 2324, 10032, 11744, 9127, 11418, 5539, 1139, 5910, 42, 1346, 2640, 4063, 479, 12059, 7918, 8985, 3233, 5256, 558, 558, 3822, 12213, 9172, 3133, 10479, 1098, 10401, 9435, 8140, 2080, 496, 1102, 11196, 11833, 3746, 5084, 2450, 10632, 10827, 11758, 2774, 10746, 6895, 10084, 4783, 6748, 4277, 2795, 8262, 169, 1124, 4326, 3907, 7908, 7018, 851, 784, 9727, 334, 70, 9294, 6292, 1727, 7994, 1014, 10604, 1108, 4655, 8372, 5082, 6827, 3497, 4738, 2472, 1854, 5403, 8115, 9612, 3001, 2400, 7070, 11823, 10079, 9199, 8261, 6641, 2694, 12245, 4088, 6872, 4029, 7775, 2272, 11849, 5693, 4418, 9368, 8269, 10092, 10033, 7760, 9373, 10797, 3102, 8029, 7519, 9739, 7061, 10928, 4352, 5501, 1589, 8951, 2549, 8292, 4028, 1908, 11074, 9098, 2141, 7030, 1154, 363, 7729, 6896, 5295, 172, 4280, 11122, 9405, 12019, 3271, 1291, 8, 2137, 6287, 7814, 11789, 2459, 7235, 3661, 2130, 3321, 9140, 11072, 4038, 4315, 10383, 713, 5289, 6816, 5969, 9217, 4179, 6201, 5531, 6746, 3073, 8513, 8225, 6694, 10339, 10681, 8306, 930, 10196, 4158, 11243, 11845, 4870, 2526, 3142, 6856, 6789, 4287, 5162, 1345, 4093, 11411, 12036, 5441, 7437, 1118, 6179, 6505, 164, 4723, 7495, 5033, 7428, 5333, 4865, 8862, 786, 514, 7759, 3125, 3657, 4587, 9773, 9538, 9471, 9179, 10743, 4995, 2502, 4297, 7403, 8259, 10369, 997, 3491, 8889, 11407, 4623, 4114, 3217, 7286, 8784, 2847, 6592, 6955, 11865, 3514, 1777, 1780, 1809, 2565, 1553, 5859, 1716, 5609, 5244, 4663, 2462, 8149, 3538, 8513, 3576, 10205, 7323, 6590, 931, 5814, 418, 5274, 2558, 298, 4782, 12283, 10877, 7723, 895, 7092, 4761, 7311, 11257, 6061, 682, 10581, 8172, 6288, 7505, 1945, 6394, 8725, 4433, 8429, 8379, 4000, 4723, 11491, 3530, 6570};

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
* Description: Multiply two polynomials pointwise (i.e., coefficient-wise).
               Inputs are assumed to be NEWHOPE_N/2 Gaussian integers.
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_mul(poly *r, const poly *a, const poly *b)
{

  int i;

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {             
    r->coeffs[i] = (montgomery_reduce((uint32_t)a->coeffs[i] * b->coeffs[i]) +
      (3*NEWHOPE_Q - montgomery_reduce((uint32_t)a->coeffs[i+1] * b->coeffs[i+1]))
    ) % NEWHOPE_Q;

    r->coeffs[i+1] = (montgomery_reduce((uint32_t)a->coeffs[i] * b->coeffs[i+1]) + 
      montgomery_reduce((uint32_t)a->coeffs[i+1] * b->coeffs[i])
    ) % NEWHOPE_Q;
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
    r->coeffs[i] = (a->coeffs[i] + 3*NEWHOPE_Q - b->coeffs[i]) % NEWHOPE_Q;
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
      /* TODO: maybe this procedure requires x/nthroots being in 
      the Montgomery domain by multiplying each coefficient by 3186. */
      r->coeffs[i] = (montgomery_reduce((uint32_t)copy[j] * nthroots[i]) + 
        (3*NEWHOPE_Q - montgomery_reduce((uint32_t)copy[NEWHOPE_K2+j] * nthroots[i+1]))
      ) % NEWHOPE_Q;
      
      r->coeffs[i+1] = (montgomery_reduce((uint32_t)copy[j] * nthroots[i+1]) + 
        montgomery_reduce((uint32_t)copy[NEWHOPE_K2+j] * nthroots[i])
      ) % NEWHOPE_Q;

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
      (3*NEWHOPE_Q - montgomery_reduce((uint32_t)copy[i+1] * invnthroots[i+1]))) % NEWHOPE_Q;

    r->coeffs[j+NEWHOPE_K2] = 
      (montgomery_reduce((uint32_t)copy[i] * invnthroots[i+1]) + 
      montgomery_reduce((uint32_t)copy[i+1] * invnthroots[i])) % NEWHOPE_Q;

    j++;
  }
}
