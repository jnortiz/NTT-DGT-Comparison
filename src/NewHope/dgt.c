#include "inttypes.h"
#include "dgt.h"
#include "params.h"
#include "reduce.h"

#if (NEWHOPE_N == 512)
static uint16_t gj[128] = {4075, 5315, 4324, 4916, 10120, 11767, 7210, 9027, 1973, 5574, 11011, 2344, 8775, 1041, 1018, 6364, 3789, 147, 5456, 7840, 7540, 5537, 4789, 4467, 9606, 1190, 8471, 6118, 5445, 3860, 7753, 1050, 7083, 5529, 9090, 12233, 8724, 11635, 10587, 1987, 6427, 6136, 6874, 3643, 400, 1728, 4948, 6137, 10256, 3998, 10367, 8410, 1254, 11316, 5435, 1359, 6950, 5446, 6093, 3710, 11907, 316, 8301, 468, 6415, 677, 6234, 3336, 12237, 9115, 1323, 2766, 12138, 10162, 8332, 9450, 2505, 5906, 10710, 11858, 5241, 9369, 9162, 8120, 787, 8807, 1010, 6821, 2049, 7377, 10968, 192, 3445, 7509, 7591, 7232, 11286, 3532, 12048, 12231, 7280, 1956, 11404, 6008, 8851, 2844, 975, 4212, 5681, 8812, 12147, 11184, 3600, 3263, 7665, 6077, 421, 8209, 6068, 3602, 8076, 11785, 605, 9987, 9260, 5594, 6403, 7507};
static uint16_t invgj[128] = {4075, 6974, 7373, 7965, 3262, 5079, 522, 2169, 5925, 11271, 11248, 3514, 9945, 1278, 6715, 10316, 11239, 4536, 8429, 6844, 6171, 3818, 11099, 2683, 7822, 7500, 6752, 4749, 4449, 6833, 12142, 8500, 11821, 3988, 11973, 382, 8579, 6196, 6843, 5339, 10930, 6854, 973, 11035, 3879, 1922, 8291, 2033, 6152, 7341, 10561, 11889, 8646, 5415, 6153, 5862, 10302, 1702, 654, 3565, 56, 3199, 6760, 5206, 4782, 5886, 6695, 3029, 2302, 11684, 504, 4213, 8687, 6221, 4080, 11868, 6212, 4624, 9026, 8689, 1105, 142, 3477, 6608, 8077, 11314, 9445, 3438, 6281, 885, 10333, 5009, 58, 241, 8757, 1003, 5057, 4698, 4780, 8844, 12097, 1321, 4912, 10240, 5468, 11279, 3482, 11502, 4169, 3127, 2920, 7048, 431, 1579, 6383, 9784, 2839, 3957, 2127, 151, 9523, 10966, 3174, 52, 8953, 6055, 11612, 5874};
#elif (NEWHOPE_N == 1024)
static uint16_t gj[256] = {4075, 5315, 4324, 4916, 10120, 11767, 7210, 9027, 1973, 5574, 11011, 2344, 8775, 1041, 1018, 6364, 3789, 147, 5456, 7840, 7540, 5537, 4789, 4467, 9606, 1190, 8471, 6118, 5445, 3860, 7753, 1050, 7083, 5529, 9090, 12233, 8724, 11635, 10587, 1987, 6427, 6136, 6874, 3643, 400, 1728, 4948, 6137, 10256, 3998, 10367, 8410, 1254, 11316, 5435, 1359, 6950, 5446, 6093, 3710, 11907, 316, 8301, 468, 6415, 677, 6234, 3336, 12237, 9115, 1323, 2766, 12138, 10162, 8332, 9450, 2505, 5906, 10710, 11858, 5241, 9369, 9162, 8120, 787, 8807, 1010, 6821, 2049, 7377, 10968, 192, 3445, 7509, 7591, 7232, 11286, 3532, 12048, 12231, 7280, 1956, 11404, 6008, 8851, 2844, 975, 4212, 5681, 8812, 12147, 11184, 3600, 3263, 7665, 6077, 421, 8209, 6068, 3602, 8076, 11785, 605, 9987, 9260, 5594, 6403, 7507, 5297, 6170, 3956, 1360, 11089, 7105, 9734, 6167, 10695, 1962, 5106, 6328, 9597, 168, 7991, 8960, 3728, 8240, 6299, 1159, 1146, 11341, 11964, 10885, 8527, 2919, 8273, 8212, 5766, 11637, 295, 6190, 8049, 8719, 11454, 6224, 8243, 709, 1319, 9139, 1958, 7967, 10211, 11177, 8210, 1058, 11848, 11367, 6507, 1566, 2948, 9786, 11606, 9830, 8633, 12225, 10542, 9166, 9235, 5486, 3834, 5257, 7856, 5919, 10314, 3757, 9364, 11942, 7535, 10431, 426, 3315, 2738, 6421, 2655, 6554, 723, 174, 1693, 9280, 350, 1512, 10474, 6906, 9087, 7796, 5369, 2057, 11026, 12240, 6374, 1483, 1583, 6347, 2500, 10800, 6142, 2447, 3963, 11713, 1954, 2051, 1805, 2882, 9928, 10446, 9259, 4115, 9381, 218, 8760, 3434, 156, 9522, 8320, 3991, 5876, 2281, 10258, 6956, 4774, 6860, 4737, 1293, 11871, 8517, 6381, 11836};
static uint16_t invgj[256] = {4075, 6974, 7373, 7965, 3262, 5079, 522, 2169, 5925, 11271, 11248, 3514, 9945, 1278, 6715, 10316, 11239, 4536, 8429, 6844, 6171, 3818, 11099, 2683, 7822, 7500, 6752, 4749, 4449, 6833, 12142, 8500, 11821, 3988, 11973, 382, 8579, 6196, 6843, 5339, 10930, 6854, 973, 11035, 3879, 1922, 8291, 2033, 6152, 7341, 10561, 11889, 8646, 5415, 6153, 5862, 10302, 1702, 654, 3565, 56, 3199, 6760, 5206, 4782, 5886, 6695, 3029, 2302, 11684, 504, 4213, 8687, 6221, 4080, 11868, 6212, 4624, 9026, 8689, 1105, 142, 3477, 6608, 8077, 11314, 9445, 3438, 6281, 885, 10333, 5009, 58, 241, 8757, 1003, 5057, 4698, 4780, 8844, 12097, 1321, 4912, 10240, 5468, 11279, 3482, 11502, 4169, 3127, 2920, 7048, 431, 1579, 6383, 9784, 2839, 3957, 2127, 151, 9523, 10966, 3174, 52, 8953, 6055, 11612, 5874, 453, 5908, 3772, 418, 10996, 7552, 5429, 7515, 5333, 2031, 10008, 6413, 8298, 3969, 2767, 12133, 8855, 3529, 12071, 2908, 8174, 3030, 1843, 2361, 9407, 10484, 10238, 10335, 576, 8326, 9842, 6147, 1489, 9789, 5942, 10706, 10806, 5915, 49, 1263, 10232, 6920, 4493, 3202, 5383, 1815, 10777, 11939, 3009, 10596, 12115, 11566, 5735, 9634, 5868, 9551, 8974, 11863, 1858, 4754, 347, 2925, 8532, 1975, 6370, 4433, 7032, 8455, 6803, 3054, 3123, 1747, 64, 3656, 2459, 683, 2503, 9341, 10723, 5782, 922, 441, 11231, 4079, 1112, 2078, 4322, 10331, 3150, 10970, 11580, 4046, 6065, 835, 3570, 4240, 6099, 11994, 652, 6523, 4077, 4016, 9370, 3762, 1404, 325, 948, 11143, 11130, 5990, 4049, 8561, 3329, 4298, 12121, 2692, 5961, 7183, 10327, 1594, 6122, 2555, 5184, 1200, 10929, 8333, 6119, 6992};
#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif

#if (NEWHOPE_N == 512)

void dgt(uint16_t *x)
{
  int m, distance;
  int j1, j2, j, k;
  uint16_t a, temp_re, temp_img;
  
  distance = 256;
  for(m = 1; m < 256; m <<= 1)
  {
    // Even level
    for(k = 0; k < m; ++k)
    {
      j1 = k*distance << 1;
      j2 = j1+distance-1;

      a = gj[k];
      for(j = j1; j <= j2; j = j+2)
      {
        temp_re =  montgomery_reduce((uint32_t)a * x[j+distance]);
        temp_img = montgomery_reduce((uint32_t)a * x[j+distance+1]);

        x[j+distance] = x[j] + 3*NEWHOPE_Q - temp_re;
        x[j+distance+1] = x[j+1] + 3*NEWHOPE_Q - temp_img;
        
        x[j] = x[j] + temp_re;
        x[j+1] = x[j+1] + temp_img;
      }
    }
    distance >>= 1;
    m <<= 1;
    // Odd level
    for(k = 0; k < m; ++k)
    {
      j1 = k*distance << 1;
      j2 = j1+distance-1;

      a = gj[k];
      for(j = j1; j <= j2; j = j+2)
      {
        temp_re =  montgomery_reduce((uint32_t)a * x[j+distance]);
        temp_img = montgomery_reduce((uint32_t)a * x[j+distance+1]);
        
        // NewHope forward NTT adopts the Gentleman-Sande butterflies, which support laziness after subtractions in all levels, differently of the Cooley-Tukey butterfly
        x[j+distance] = (x[j] + 3*NEWHOPE_Q - temp_re) % NEWHOPE_Q;         
        x[j+distance+1] = (x[j+1] + 3*NEWHOPE_Q - temp_img) % NEWHOPE_Q;
        
        x[j] = (x[j] + temp_re) % NEWHOPE_Q;
        x[j+1] = (x[j+1] + temp_img) % NEWHOPE_Q;
      }
    }
    distance >>= 1;
  }  
}

void idgt(uint16_t *poly)
{
  int distance, j1, jtwiddle, j;
  uint16_t sub_re, sub_img;

  for(distance = 2; distance < NEWHOPE_N; distance <<= 1)
  {
    // Even level
    for(j1 = 0; j1 < distance; j1 += 2)
    {
      jtwiddle = 0;
      for(j = j1; j < NEWHOPE_N; j += distance << 1)
      {
        sub_re = poly[j] + 3*NEWHOPE_Q - poly[j+distance]; // Omit reduction (be lazy)
        sub_img = poly[j+1] + 3*NEWHOPE_Q - poly[j+distance+1]; // Omit reduction (be lazy)

        poly[j] = poly[j] + poly[j+distance]; // Omit reduction (be lazy)
        poly[j+1] = poly[j+1] + poly[j+distance+1]; // Omit reduction (be lazy)

        poly[j+distance] = montgomery_reduce((uint32_t)invgj[jtwiddle] * sub_re);
        poly[j+distance+1] = montgomery_reduce((uint32_t)invgj[jtwiddle++] * sub_img);
      }
    }
    distance <<= 1;
    // Odd level
    for(j1 = 0; j1 < distance; j1 += 2)
    {
      jtwiddle = 0;
      for(j = j1; j < NEWHOPE_N; j += distance << 1)
      {
        sub_re = poly[j] + 3*NEWHOPE_Q - poly[j+distance]; // Omit reduction (be lazy)
        sub_img = poly[j+1] + 3*NEWHOPE_Q - poly[j+distance+1]; // Omit reduction (be lazy)

        poly[j] = (poly[j] + poly[j+distance]) % NEWHOPE_Q;
        poly[j+1] = (poly[j+1] + poly[j+distance+1]) % NEWHOPE_Q;

        poly[j+distance] = montgomery_reduce((uint32_t)invgj[jtwiddle] * sub_re);
        poly[j+distance+1] = montgomery_reduce((uint32_t)invgj[jtwiddle++] * sub_img);
      }
    }    
  }
}

#elif (NEWHOPE_N == 1024)

void dgt(uint16_t *x)
{
  int m, distance;
  int j1, j2, j, k;
  uint16_t a, temp_re, temp_img;
  
  distance = 512;
  for(m = 1; m < 512; m <<= 1)
  {
    // Even level
    for(k = 0; k < m; ++k)
    {
      j1 = k*distance << 1;
      j2 = j1+distance-1;

      a = gj[k];
      for(j = j1; j <= j2; j = j+2)
      {
        temp_re = montgomery_reduce((uint32_t)a * x[j+distance]);
        temp_img = montgomery_reduce((uint32_t)a * x[j+distance+1]);

        x[j+distance] = x[j] + 3*NEWHOPE_Q - temp_re; // Omit reduction (be lazy)
        x[j+distance+1] = x[j+1] + 3*NEWHOPE_Q - temp_img; // Omit reduction (be lazy)
        
        x[j] = x[j] + temp_re; // Omit reduction (be lazy)
        x[j+1] = x[j+1] + temp_img; // Omit reduction (be lazy)
      }
    }
    distance >>= 1;
    m <<= 1;
    if(m < 512)
    {
      // Odd level
      for(k = 0; k < m; ++k)
      {
        j1 = k*distance << 1;
        j2 = j1+distance-1;

        a = gj[k];
        for(j = j1; j <= j2; j = j+2)
        {
          temp_re = montgomery_reduce((uint32_t)a * x[j+distance]);
          temp_img = montgomery_reduce((uint32_t)a * x[j+distance+1]);

          x[j+distance] = (x[j] + 3*NEWHOPE_Q - temp_re) % NEWHOPE_Q;
          x[j+distance+1] = (x[j+1] + 3*NEWHOPE_Q - temp_img) % NEWHOPE_Q;
          
          x[j] = (x[j] + temp_re) % NEWHOPE_Q;
          x[j+1] = (x[j+1] + temp_img) % NEWHOPE_Q;
        }
      }
      distance >>= 1;      
    }
  }  
}

void idgt(uint16_t *poly)
{
  int distance, j1, jtwiddle, j;
  uint16_t sub_re, sub_img;

  for(distance = 2; distance < NEWHOPE_N; distance <<= 1)
  {
    // Even level
    for(j1 = 0; j1 < distance; j1 += 2)
    {
      jtwiddle = 0;
      for(j = j1; j < NEWHOPE_N; j += distance << 1)
      {
        sub_re = poly[j] + 3*NEWHOPE_Q - poly[j+distance]; // Omit reduction (be lazy)
        sub_img = poly[j+1] + 3*NEWHOPE_Q - poly[j+distance+1]; // Omit reduction (be lazy)

        poly[j] = poly[j] + poly[j+distance]; // Omit reduction (be lazy)
        poly[j+1] = poly[j+1] + poly[j+distance+1]; // Omit reduction (be lazy)

        poly[j+distance] = montgomery_reduce((uint32_t)invgj[jtwiddle] * sub_re);
        poly[j+distance+1] = montgomery_reduce((uint32_t)invgj[jtwiddle++] * sub_img);
      }
    }
    distance <<= 1;
    if(distance < NEWHOPE_N)
    {
      // Odd level
      for(j1 = 0; j1 < distance; j1 += 2)
      {
        jtwiddle = 0;
        for(j = j1; j < NEWHOPE_N; j += distance << 1)
        {
          sub_re = poly[j] + 3*NEWHOPE_Q - poly[j+distance]; // Omit reduction (be lazy)
          sub_img = poly[j+1] + 3*NEWHOPE_Q - poly[j+distance+1]; // Omit reduction (be lazy)

          poly[j] = (poly[j] + poly[j+distance]) % NEWHOPE_Q;
          poly[j+1] = (poly[j+1] + poly[j+distance+1]) % NEWHOPE_Q;

          poly[j+distance] = montgomery_reduce((uint32_t)invgj[jtwiddle] * sub_re);
          poly[j+distance+1] = montgomery_reduce((uint32_t)invgj[jtwiddle++] * sub_img);
        }
      }
    }
  }
}

#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif
