/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: DGT
**************************************************************************************/

#include "dgt.h"
#include "params.h"
#include "reduce.h"

static uint16_t gj[256] = {4075, 5297, 6415, 10314, 7083, 8049, 11286, 6142, 3789, 3728, 5241, 350, 10256, 6507, 3600, 156, 1973, 10695, 12138, 2738, 6427, 1958, 8851, 9928, 9606, 8527, 2049, 11026, 6950, 10542, 8076, 4774, 10120, 11089, 12237, 7535, 8724, 8243, 7280, 1954, 7540, 1146, 787, 9087, 1254, 11606, 421, 5876, 8775, 9597, 2505, 723, 400, 8210, 5681, 9381, 5445, 5766, 3445, 1583, 11907, 3834, 9260, 11871, 4324, 3956, 6234, 9364, 9090, 11454, 12048, 3963, 5456, 6299, 9162, 10474, 10367, 2948, 7665, 8320, 11011, 5106, 8332, 2655, 6874, 10211, 975, 9259, 8471, 8273, 10968, 6374, 6093, 9235, 605, 4737, 7210, 9734, 1323, 426, 10587, 1319, 11404, 1805, 4789, 11964, 1010, 5369, 5435, 8633, 6068, 10258, 1018, 7991, 10710, 1693, 4948, 11848, 12147, 8760, 7753, 295, 7591, 2500, 8301, 7856, 6403, 6381, 5315, 6170, 677, 3757, 5529, 8719, 3532, 2447, 147, 8240, 9369, 1512, 3998, 1566, 3263, 9522, 5574, 1962, 10162, 6421, 6136, 7967, 2844, 10446, 1190, 2919, 7377, 12240, 5446, 9166, 11785, 6860, 11767, 7105, 9115, 10431, 11635, 709, 1956, 2051, 5537, 11341, 8807, 7796, 11316, 9830, 8209, 2281, 1041, 168, 5906, 174, 1728, 1058, 8812, 218, 3860, 11637, 7509, 6347, 316, 5257, 5594, 8517, 4916, 1360, 3336, 11942, 12233, 6224, 12231, 11713, 7840, 1159, 8120, 6906, 8410, 9786, 6077, 3991, 2344, 6328, 9450, 6554, 3643, 11177, 4212, 4115, 6118, 8212, 192, 1483, 3710, 5486, 9987, 1293, 9027, 6167, 2766, 3315, 1987, 9139, 6008, 2882, 4467, 10885, 6821, 2057, 1359, 12225, 3602, 6956, 6364, 8960, 11858, 9280, 6137, 11367, 11184, 3434, 1050, 6190, 7232, 10800, 468, 5919, 7507, 11836};
static uint16_t invgj[256] = {4075, 453, 4782, 6370, 11821, 1489, 5057, 6099, 11239, 8855, 1105, 922, 6152, 3009, 431, 3329, 5925, 5333, 8687, 64, 10930, 10232, 5468, 1404, 7822, 9407, 6281, 3150, 10302, 8974, 9523, 6122, 3262, 10996, 2302, 6803, 8579, 10806, 12097, 4077, 6171, 8174, 8077, 1112, 8646, 5735, 2839, 5961, 9945, 8298, 6212, 2503, 3879, 5383, 4169, 11130, 4449, 576, 58, 6065, 56, 347, 8953, 10929, 7373, 3772, 6695, 7032, 11973, 5942, 4780, 652, 8429, 12071, 3477, 11231, 10561, 12115, 6383, 12121, 11248, 10008, 4080, 2459, 973, 4493, 3482, 948, 6752, 10238, 10333, 11580, 654, 1858, 3174, 5184, 522, 5429, 504, 3123, 6843, 49, 4912, 9370, 11099, 1843, 9445, 4322, 6153, 5868, 2127, 10327, 6715, 2767, 9026, 10723, 8291, 10777, 2920, 4049, 12142, 9842, 8757, 3570, 6760, 8532, 11612, 6119, 6974, 5908, 5886, 4433, 3988, 9789, 4698, 11994, 4536, 3529, 142, 441, 7341, 10596, 1579, 4298, 11271, 2031, 6221, 3656, 6854, 6920, 11279, 325, 7500, 10484, 885, 10970, 1702, 11863, 10966, 2555, 5079, 7552, 11684, 3054, 6196, 5915, 1321, 4016, 3818, 3030, 11314, 2078, 5415, 9634, 3957, 7183, 1278, 3969, 4624, 9341, 1922, 1815, 3127, 5990, 6833, 8326, 241, 835, 3199, 2925, 6055, 8333, 7965, 418, 3029, 8455, 382, 10706, 8844, 6523, 6844, 2908, 6608, 4079, 11889, 11566, 9784, 2692, 3514, 6413, 11868, 683, 11035, 3202, 11502, 11143, 4749, 10335, 5009, 4046, 3565, 4754, 52, 1200, 2169, 7515, 4213, 1747, 5339, 1263, 10240, 3762, 2683, 2361, 3438, 10331, 5862, 9551, 151, 1594, 10316, 12133, 8689, 5782, 2033, 11939, 7048, 8561, 8500, 6147, 1003, 4240, 5206, 1975, 5874, 6992};


void dgt(uint16_t *poly)
{
  int i, index, j, m, window;
  uint64_t a, sub_re, sub_img;
  uint64_t copy[NEWHOPE_N];

  for(i = 0; i < NEWHOPE_N; ++i)
  {
    copy[i] = (uint64_t) poly[i];
  }

  window = 1;
  for(m = NEWHOPE_K2; m >= 2; m >>= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = (uint64_t) gj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        sub_re = copy[i] + (NEWHOPE_3Q - copy[i+m]);
        sub_img = copy[i+1] + (NEWHOPE_3Q - copy[i+m+1]);
        
        copy[i] = (copy[i] + copy[i+m]) % NEWHOPE_Q;
        copy[i+1] = (copy[i+1] + copy[i+m+1]) % NEWHOPE_Q;
        
        copy[i+m] = montgomery_reduce((uint32_t)a * sub_re);
        copy[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;

    m >>= 1;

    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = (uint64_t)gj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        sub_re = copy[i] + (NEWHOPE_3Q - copy[i+m]);
        sub_img = copy[i+1] + (NEWHOPE_3Q - copy[i+m+1]);
        
        copy[i] = (copy[i] + copy[i+m]) % NEWHOPE_Q;
        copy[i+1] = (copy[i+1] + copy[i+m+1]) % NEWHOPE_Q;
        
        copy[i+m] = montgomery_reduce((uint32_t)a * sub_re);
        copy[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;

    m >>= 1;

    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = (uint64_t)gj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        sub_re = copy[i] + (NEWHOPE_3Q - copy[i+m]);
        sub_img = copy[i+1] + (NEWHOPE_3Q - copy[i+m+1]);
        
        copy[i] = copy[i] + copy[i+m];
        copy[i+1] = copy[i+1] + copy[i+m+1];
        
        copy[i+m] = montgomery_reduce((uint32_t)a * sub_re);
        copy[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;
  }

  for(i = 0; i < NEWHOPE_N; ++i)
  {
    poly[i] = (uint16_t) copy[i];
  }

}

void idgt(uint16_t *poly)
{
  int i, index, j, m, window;
  uint64_t a, mul_re, mul_img;
  uint64_t copy[NEWHOPE_N];

  for(i = 0; i < NEWHOPE_N; ++i)
  {
    copy[i] = (uint64_t) poly[i];
  }

  window = (NEWHOPE_K2 >> 1);
  for(m = 2; m <= 256; m <<= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = (uint64_t)invgj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        mul_re = montgomery_reduce((uint32_t)copy[i+m] * a);
        mul_img = montgomery_reduce((uint32_t)copy[i+m+1] * a);
        
        copy[i+m] = (copy[i] + (NEWHOPE_3Q - mul_re)) % NEWHOPE_Q;
        copy[i+m+1] = (copy[i+1] + (NEWHOPE_3Q - mul_img)) % NEWHOPE_Q;
        
        copy[i] = (copy[i] + mul_re) % NEWHOPE_Q;
        copy[i+1] = (copy[i+1] + mul_img) % NEWHOPE_Q;        
      }
      index += window;
    }
    window >>= 1;

    m <<= 1;

    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = (uint64_t)invgj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        mul_re = montgomery_reduce((uint32_t)copy[i+m] * a);
        mul_img = montgomery_reduce((uint32_t)copy[i+m+1] * a);
        
        copy[i+m] = copy[i] + (NEWHOPE_3Q - mul_re);
        copy[i+m+1] = copy[i+1] + (NEWHOPE_3Q - mul_img);
        
        copy[i] = copy[i] + mul_re;
        copy[i+1] = copy[i+1] + mul_img;        
      }
      index += window;
    }
    window >>= 1;
  }

  index = 0;
  for(j = 0; j < m; j += 2) 
  {
    a = (uint64_t)invgj[index];
    for(i = j; i < NEWHOPE_N; i += (m << 1)) 
    {
      mul_re = montgomery_reduce((uint32_t)copy[i+m] * a);
      mul_img = montgomery_reduce((uint32_t)copy[i+m+1] * a);
      
      copy[i+m] = (copy[i] + (NEWHOPE_3Q - mul_re)) % NEWHOPE_Q;
      copy[i+m+1] = (copy[i+1] + (NEWHOPE_3Q - mul_img)) % NEWHOPE_Q;
      
      copy[i] = (copy[i] + mul_re) % NEWHOPE_Q;
      copy[i+1] = (copy[i+1] + mul_img) % NEWHOPE_Q;        
    }
    index += window;
  }  

  for(i = 0; i < NEWHOPE_N; ++i)
  {
    poly[i] = (uint16_t) copy[i];
  }

}