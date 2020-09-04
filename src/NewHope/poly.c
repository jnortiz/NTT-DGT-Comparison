#include "poly.h"
#include "dgt.h"
#include "reduce.h"
#include "fips202.h"

#if (NEWHOPE_N == 512)
static const uint16_t nthroots[NEWHOPE_N] = {4075, 0, 8633, 11790, 535, 4086, 1872, 9399, 8606, 9351, 8602, 5409, 9875, 3346, 3052, 3506, 11298, 9391, 3562, 6874, 12239, 5951, 4118, 9833, 5860, 12221, 9671, 7902, 11296, 7616, 382, 2565, 1820, 4956, 4731, 1721, 10330, 7165, 6100, 6752, 2646, 11766, 4278, 3369, 142, 8061, 2645, 4183, 407, 6635, 2688, 2372, 6549, 8538, 4407, 7102, 6554, 8529, 674, 10861, 6817, 7221, 7983, 11176, 663, 8468, 3877, 10657, 9802, 2057, 6675, 754, 696, 590, 636, 5254, 8758, 9735, 1091, 8719, 4341, 10687, 11428, 11907, 4625, 7523, 5840, 11659, 6716, 4057, 2469, 6723, 11449, 10768, 495, 8482, 5097, 1446, 9784, 7385, 1855, 4417, 2791, 1760, 5489, 7642, 9609, 1570, 11297, 11284, 6156, 5797, 9427, 12242, 840, 6032, 5935, 7281, 5199, 1212, 95, 5384, 7091, 1886, 10947, 7693, 11929, 877, 10919, 3624, 4571, 1066, 12248, 3309, 1021, 2555, 6722, 7628, 1960, 10922, 11681, 8200, 473, 8541, 3595, 4514, 8579, 9159, 10002, 3172, 9301, 9108, 3781, 183, 12149, 2096, 3787, 7995, 6336, 3437, 4845, 1415, 8990, 9235, 11838, 656, 10286, 5673, 5040, 4260, 2936, 9222, 9026, 833, 11383, 7867, 2316, 9404, 4125, 4464, 11098, 9531, 11328, 9291, 7202, 39, 10383, 4249, 10424, 9615, 11370, 5233, 12159, 10279, 7455, 3116, 7306, 420, 4893, 4386, 3599, 4911, 7083, 10467, 4909, 6467, 5477, 10144, 8554, 5331, 10002, 1145, 10054, 7359, 6875, 9733, 2177, 4847, 1448, 8865, 2582, 4695, 9893, 2082, 4980, 8127, 1090, 10760, 1169, 4251, 4704, 9110, 10484, 3024, 2430, 9996, 3670, 3744, 3728, 7413, 7150, 4277, 270, 5037, 6219, 5077, 1734, 5727, 1721, 2421, 8849, 11022, 6360, 8093, 7580, 8794, 7669, 4620, 11873, 10097, 1486, 1719, 8422, 8686, 7793, 814, 8859, 9757, 2140, 6211, 8518, 7688, 353, 1001, 8087, 4493, 8328, 2334, 5082, 306, 852, 418, 10656, 10917, 3393, 1835, 3022, 6542, 11181, 3245, 5532, 1000, 4517, 6604, 2875, 9092, 2765, 8810, 3951, 3609, 4542, 1411, 7081, 12023, 11950, 6178, 11302, 2228, 5317, 2932, 55, 6821, 7523, 10653, 10294, 4087, 113, 10685, 8496, 2532, 8888, 2837, 9216, 8975, 3654, 940, 1288, 2702, 3533, 9303, 8250, 9285, 7125, 8045, 9234, 8285, 7133, 2546, 3475, 45, 10172, 11612, 2311, 4395, 422, 10192, 5969, 1118, 5602, 757, 7796, 776, 10646, 4275, 12260, 1981, 6490, 292, 3224, 9447, 2619, 10798, 10004, 3324, 1480, 1025, 969, 7509, 599, 2086, 7393, 12238, 11042, 7310, 9710, 5244, 7570, 1754, 7796, 1076, 5933, 387, 3666, 12123, 10397, 10834, 6456, 9779, 4181, 8114, 6365, 1940, 91, 4190, 10908, 9924, 3561, 8832, 5824, 11405, 7227, 3165, 11826, 11933, 9627, 6311, 4295, 1983, 1296, 5768, 2698, 8035, 3108, 1909, 11206, 7373, 819, 11428, 7459, 6206, 3688, 10084, 3680, 6391, 7852, 55, 4657, 11058, 644, 2232, 1321, 9701, 434, 11159, 576, 3521, 8865, 10115, 2091, 11251, 6658, 11770, 4325, 3426, 2447, 9904, 4895, 11225, 466, 4544, 1224, 2748, 3394, 328, 1000, 3108, 1133, 1483, 4907, 1652, 5815, 993, 5752, 695, 4614, 1724, 3074, 12047, 584, 5468, 8282, 9094, 9588, 2175, 8332, 7906, 11382, 8280, 1125, 472, 3566, 3323, 2473, 421, 9386, 8448, 759, 3313, 9761, 9664, 4323, 5611, 3416, 8345, 6383, 1552, 8411, 5323, 4116, 6990, 779, 3074, 8218, 3736, 1441, 7198, 9076, 12057, 10876, 8937, 2154, 6374};
#elif (NEWHOPE_N == 1024)
static const uint16_t nthroots[NEWHOPE_N] = {4075, 0, 1650, 10001, 9474, 9833, 6082, 4793, 10187, 91, 7653, 10614, 5174, 3698, 2766, 6985, 8525, 9484, 550, 7196, 9821, 7807, 9509, 12068, 6763, 1969, 3473, 8900, 3840, 5170, 2603, 1988, 599, 10203, 10341, 6684, 3709, 6843, 8040, 4090, 11719, 4470, 1706, 11735, 434, 5383, 2788, 10012, 5532, 11599, 4054, 5342, 11384, 7028, 3275, 12220, 6408, 6728, 8631, 4060, 1869, 7962, 8834, 10914, 4162, 7309, 11736, 3575, 11509, 10199, 6442, 11172, 4813, 9619, 12030, 9172, 6782, 7846, 8035, 4553, 11989, 9983, 4772, 2552, 850, 10263, 238, 2849, 8608, 12166, 1806, 7590, 8672, 6815, 8529, 5502, 7948, 10687, 7171, 9337, 8134, 8144, 9757, 8993, 4006, 7946, 31, 10926, 7453, 6976, 2209, 10932, 416, 1235, 4345, 4781, 9489, 597, 1994, 306, 8814, 2014, 3228, 731, 8346, 10474, 6145, 10508, 11823, 4544, 1087, 10253, 7174, 8303, 2437, 112, 4958, 9998, 4331, 8083, 8973, 8688, 2115, 10817, 6963, 11953, 7869, 7779, 3730, 6114, 10483, 3201, 10665, 11502, 1099, 3057, 5346, 5470, 9865, 516, 6178, 11950, 164, 242, 8981, 1176, 5301, 12026, 4927, 3226, 12199, 10856, 10878, 1717, 1406, 10423, 3056, 949, 2877, 11265, 9812, 3727, 4229, 4079, 3119, 1744, 3756, 12218, 4367, 10082, 2245, 7327, 4845, 10874, 7063, 6553, 6310, 2973, 8021, 10161, 3790, 7925, 556, 11274, 9019, 4108, 8187, 2851, 7654, 9688, 4178, 7068, 2551, 1273, 12179, 1170, 2075, 4471, 11698, 1550, 6587, 7980, 2217, 7838, 2898, 991, 608, 2592, 6786, 12201, 9170, 337, 9620, 3100, 3374, 875, 8492, 2706, 6912, 3306, 1784, 11719, 6102, 2697, 10698, 10736, 5046, 9773, 9367, 6121, 9621, 6987, 5394, 6420, 2577, 545, 1455, 1892, 5982, 5106, 7394, 9321, 667, 5648, 9726, 296, 6345, 3305, 6958, 5134, 8782, 8911, 6412, 2655, 4283, 11133, 5087, 6509, 1031, 1827, 8880, 9591, 8085, 511, 7671, 9223, 10875, 6837, 5139, 4277, 1364, 4917, 7719, 2169, 1786, 6867, 6346, 4905, 11551, 1309, 9625, 9702, 10373, 9230, 4678, 12286, 11018, 3859, 2020, 11076, 12263, 3622, 1423, 6632, 3941, 8343, 10484, 6660, 9432, 11747, 1446, 5097, 8954, 1650, 8397, 5496, 8048, 8515, 3501, 8326, 3327, 11861, 5558, 10501, 2290, 7799, 5659, 3772, 9162, 6426, 2366, 11719, 11012, 9816, 4928, 7077, 5791, 2288, 2963, 10935, 2454, 2375, 1724, 4614, 83, 849, 10438, 5834, 11715, 3878, 11743, 1793, 4311, 9757, 10710, 5336, 7221, 7818, 510, 37, 366, 5368, 8120, 1959, 9807, 11584, 4133, 11455, 5964, 4688, 8301, 549, 8395, 5640, 8888, 9452, 10085, 2519, 151, 10044, 3671, 12025, 1302, 2452, 11904, 3757, 6646, 8764, 3049, 11901, 3568, 5242, 9774, 1042, 8647, 5287, 7287, 10401, 6402, 8953, 6259, 7597, 1073, 6507, 638, 11580, 2885, 9973, 9509, 4665, 5022, 11939, 3194, 2690, 10554, 1627, 787, 7981, 2155, 4385, 2445, 9869, 2466, 5458, 411, 10542, 2349, 3697, 8232, 10154, 10111, 4598, 4818, 1423, 2400, 9198, 1489, 10670, 10469, 4956, 7091, 6753, 2664, 12013, 4820, 1372, 6947, 684, 4168, 8587, 4009, 1381, 734, 9400, 8063, 12224, 2791, 8402, 166, 11250, 4122, 4076, 5188, 11622, 9760, 2384, 9428, 8287, 11851, 9600, 5062, 3165, 10220, 8931, 7488, 185, 6242, 9580, 9761, 9916, 4553, 1888, 11812, 7409, 7652, 3693, 2892, 91, 10806, 11649, 1694, 10266, 6836, 7497, 8262, 6999, 8985, 5671, 3710, 10807, 4874, 9440, 4620, 7669, 4524, 3346, 11681, 7509, 6675, 8023, 3088, 11262, 3415, 4699, 1886, 9349, 1239, 10302, 8596, 6606, 7636, 3147, 1243, 7578, 10108, 4198, 10699, 10533, 9459, 619, 4259, 8444, 6547, 11064, 9427, 47, 6636, 4705, 6731, 7941, 4479, 9017, 275, 2246, 8589, 7088, 6791, 7360, 5314, 2400, 3650, 9487, 8225, 9859, 3504, 6310, 10782, 11423, 6283, 4434, 501, 4483, 8492, 11537, 2968, 4243, 3566, 8966, 7618, 3127, 2476, 2586, 5811, 3327, 8265, 8374, 6064, 9094, 5758, 9213, 11286, 9050, 604, 11990, 8394, 4570, 2753, 11456, 454, 8298, 3513, 8286, 8225, 7025, 6370, 1996, 2171, 8058, 9743, 5156, 10302, 2516, 6347, 5515, 5576, 8769, 4624, 10046, 2742, 5416, 466, 1148, 2896, 11880, 7285, 3130, 6788, 6963, 7801, 11662, 2683, 3032, 11477, 10234, 2091, 5417, 7314, 6137, 10391, 3695, 130, 10279, 11198, 4270, 855, 2049, 10972, 1963, 2315, 2246, 1348, 4676, 1009, 8724, 4378, 6443, 11102, 11538, 2994, 11810, 9517, 12142, 1280, 2480, 3690, 768, 6280, 10694, 4859, 6941, 4592, 11551, 6635, 407, 2707, 6150, 3946, 2837, 1143, 3704, 4750, 2007, 1862, 2820, 9002, 4632, 4379, 1975, 2882, 11589, 8087, 5318, 1191, 5164, 10354, 3298, 9087, 9750, 10737, 9262, 11496, 11983, 6977, 831, 11428, 819, 645, 1382, 1206, 9139, 11838, 466, 2953, 843, 11075, 12257, 7464, 1308, 4402, 6556, 3271, 3452, 8513, 8554, 7517, 7945, 7007, 5447, 9442, 6471, 6606, 7162, 2839, 177, 9138, 9298, 353, 11288, 7473, 455, 7247, 6290, 1852, 1189, 10350, 4897, 2269, 4314, 6619, 4161, 6865, 2109, 5915, 4403, 4638, 8504, 11122, 8668, 6638, 7817, 777, 4025, 4007, 9740, 2251, 545, 1323, 9940, 8665, 1370, 7399, 8446, 7717, 12031, 11080, 430, 9372, 9891, 838, 7033, 8100, 1762, 6368, 3554, 4972, 11000, 10225, 5052, 4914, 10026, 3195, 8879, 5320, 5013, 9429, 1199, 1979, 7987, 10702, 9641, 6966, 3878, 4618, 11673, 9205, 11425, 10142, 6424, 5126, 10218, 1537, 4477, 10422, 11743, 2879, 2519, 3717, 5263, 1758, 2725, 8801, 5050, 4665, 8765, 12148, 3282, 3387, 6770, 4675, 4772, 11813, 4811, 1643, 4275, 10831, 2636, 12087, 3813, 9912, 5626, 3276, 419, 7373, 12103, 2314, 4051, 3311, 8975, 9094, 10029, 6771, 2643, 7374, 2944, 8716, 4012, 3131, 11574, 806, 3613, 8700, 9385, 2625, 11732, 5331, 8554, 10149, 4623, 3038, 4717, 3384, 9767, 5681, 5758, 4143, 4953, 1289, 10177, 5078, 3587, 590, 3354, 4028, 5867, 10085, 117, 5841, 3040, 6771, 7291, 436, 9048, 2042, 1344, 4377, 10598, 663, 3821, 7821, 5611, 325, 5818, 4888, 11145, 5779, 7486, 4019, 10848, 6916, 622, 11295, 3145, 4560, 1536, 4859, 1590, 3976, 2599, 9755, 9200, 7306, 12075, 10115, 567, 2110, 10624, 3855, 614, 434, 1130, 1157, 7536, 8207, 11280, 730, 4691, 1883, 10793, 5975, 10433, 6133, 945, 7643, 2099, 7491, 3362, 536, 1896, 4050, 6082, 1454, 3129, 7418, 1108, 7136, 1730, 6017, 952, 2039, 6398, 9044, 1108, 3497, 10217, 2059, 738, 10835, 3223, 11803, 5309, 11219, 10001, 10936, 760, 9695, 11167, 11435, 2893, 10006, 3783, 3956, 5193, 11134, 11498, 4809, 9119, 12200, 12084, 6212, 3622, 2064, 5907, 8694, 4514, 11670, 8023, 10388, 9332, 643, 9846, 6823, 5604, 6835, 11077, 1119, 11112, 7811, 10198, 4543, 9807, 2765, 4077, 9036, 8482, 8967, 2016, 4621, 10073, 8911, 3390, 12131, 12238, 4633, 11193};
#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif

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
* Name:        poly_mul_pointwise
* 
* Description: Multiply two polynomials point-wise (i.e., coefficient-wise).
               Inputs are assumed to be NEWHOPE_N/2 Gaussian integers.
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_mul_pointwise(poly *r, const poly *a, const poly *b)
{
  int i;
  uint32_t t, s;

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
    t = montgomery_reduce(3186*a->coeffs[i]); // R^2 mod p = 3186
    s = montgomery_reduce(3186*a->coeffs[i+1]);

    r->coeffs[i]   = (montgomery_reduce((uint32_t)t * b->coeffs[i]) +
                      3*NEWHOPE_Q - montgomery_reduce((uint32_t)s * b->coeffs[i+1])) % NEWHOPE_Q;

    r->coeffs[i+1] = (montgomery_reduce((uint32_t)t * b->coeffs[i+1]) + 
                      montgomery_reduce((uint32_t)s * b->coeffs[i])) % NEWHOPE_Q;
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

  for(i = 0; i < NEWHOPE_N; ++i) 
  {
    copy[i] = r->coeffs[i];
  }  

  j = 0;
  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
      r->coeffs[i]   = (montgomery_reduce((uint32_t)copy[j] * nthroots[i]) + (3*NEWHOPE_Q - 
                        montgomery_reduce((uint32_t)copy[j+NEWHOPE_N/2] * nthroots[i+1]))) % NEWHOPE_Q;      
      r->coeffs[i+1] = (montgomery_reduce((uint32_t)copy[j] * nthroots[i+1] + (uint32_t)copy[j+NEWHOPE_N/2] * nthroots[i])) % NEWHOPE_Q;
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
  idgt((uint16_t *)r->coeffs);
}