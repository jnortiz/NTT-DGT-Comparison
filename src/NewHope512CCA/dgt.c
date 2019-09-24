/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include "dgt.h"
#include "params.h"
#include "reduce.h"

static uint16_t gj[128] = {4075, 6415, 7083, 11286, 3789, 5241, 10256, 3600, 1973, 12138, 6427, 8851, 9606, 2049, 6950, 8076, 10120, 12237, 8724, 7280, 7540, 787, 1254, 421, 8775, 2505, 400, 5681, 5445, 3445, 11907, 9260, 4324, 6234, 9090, 12048, 5456, 9162, 10367, 7665, 11011, 8332, 6874, 975, 8471, 10968, 6093, 605, 7210, 1323, 10587, 11404, 4789, 1010, 5435, 6068, 1018, 10710, 4948, 12147, 7753, 7591, 8301, 6403, 5315, 677, 5529, 3532, 147, 9369, 3998, 3263, 5574, 10162, 6136, 2844, 1190, 7377, 5446, 11785, 11767, 9115, 11635, 1956, 5537, 8807, 11316, 8209, 1041, 5906, 1728, 8812, 3860, 7509, 316, 5594, 4916, 3336, 12233, 12231, 7840, 8120, 8410, 6077, 2344, 9450, 3643, 4212, 6118, 192, 3710, 9987, 9027, 2766, 1987, 6008, 4467, 6821, 1359, 3602, 6364, 11858, 6137, 11184, 1050, 7232, 468, 7507, };
static uint16_t invgj[128] = {4075, 4782, 11821, 5057, 11239, 1105, 6152, 431, 5925, 8687, 10930, 5468, 7822, 6281, 10302, 9523, 3262, 2302, 8579, 12097, 6171, 8077, 8646, 2839, 9945, 6212, 3879, 4169, 4449, 58, 56, 8953, 7373, 6695, 11973, 4780, 8429, 3477, 10561, 6383, 11248, 4080, 973, 3482, 6752, 10333, 654, 3174, 522, 504, 6843, 4912, 11099, 9445, 6153, 2127, 6715, 9026, 8291, 2920, 12142, 8757, 6760, 11612, 6974, 5886, 3988, 4698, 4536, 142, 7341, 1579, 11271, 6221, 6854, 11279, 7500, 885, 1702, 10966, 5079, 11684, 6196, 1321, 3818, 11314, 5415, 3957, 1278, 4624, 1922, 3127, 6833, 241, 3199, 6055, 7965, 3029, 382, 8844, 6844, 6608, 11889, 9784, 3514, 11868, 11035, 11502, 4749, 5009, 3565, 52, 2169, 4213, 5339, 10240, 2683, 3438, 5862, 151, 10316, 8689, 2033, 7048, 8500, 1003, 5206, 5874};

int32_t reduce(int64_t a)
{ // Montgomery reduction
  int64_t u;

  u = (a*12287) & 0xFFFFFFFF;
  u *= NEWHOPE_Q;
  a += u;
  return (int32_t)(a>>32);
}

void dgt(uint16_t *poly)
{

  int i, index, j, m, window;
  uint32_t a, sub_re, sub_img;

  window = 1;
  for(m = NEWHOPE_K2; m >= 8; m >>= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        sub_re = poly[i] + (NEWHOPE_3Q - poly[i+m]);
        sub_img = poly[i+1] + (NEWHOPE_3Q - poly[i+m+1]);
        
        poly[i] = (poly[i] + poly[i+m]) % NEWHOPE_Q;
        poly[i+1] = (poly[i+1] + poly[i+m+1]) % NEWHOPE_Q;
        
        poly[i+m] = montgomery_reduce((uint32_t)a * sub_re);
        poly[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;
  
    m >>= 1;
   
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        sub_re = poly[i] + (NEWHOPE_3Q - poly[i+m]);
        sub_img = poly[i+1] + (NEWHOPE_3Q - poly[i+m+1]);
        
        poly[i] = (poly[i] + poly[i+m]);
        poly[i+1] = (poly[i+1] + poly[i+m+1]);
        
        poly[i+m] = montgomery_reduce((uint32_t)a * sub_re);
        poly[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;   

    m >>= 1;
   
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        sub_re = poly[i] + (NEWHOPE_3Q - poly[i+m]);
        sub_img = poly[i+1] + (NEWHOPE_3Q - poly[i+m+1]);
        
        poly[i] = poly[i] + poly[i+m];
        poly[i+1] = poly[i+1] + poly[i+m+1];
        
        poly[i+m] = montgomery_reduce((uint32_t)a * sub_re);
        poly[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;

  }

  index = 0;
  for(j = 0; j < m; j += 2) 
  {
    a = gj[index];
    for(i = j; i < NEWHOPE_N; i += (m << 1)) 
    {
      sub_re = poly[i] + (NEWHOPE_3Q - poly[i+m]);
      sub_img = poly[i+1] + (NEWHOPE_3Q - poly[i+m+1]);
      
      poly[i] = (poly[i] + poly[i+m]) % NEWHOPE_Q;
      poly[i+1] = (poly[i+1] + poly[i+m+1]) % NEWHOPE_Q;
      
      poly[i+m] = montgomery_reduce((uint32_t)a * sub_re);
      poly[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
    }
    index += window;
  }
  window <<= 1;

  m >>= 1;  

  index = 0;
  for(j = 0; j < m; j += 2) 
  {
    a = gj[index];
    for(i = j; i < NEWHOPE_N; i += (m << 1)) 
    {
      sub_re = poly[i] + (NEWHOPE_3Q - poly[i+m]);
      sub_img = poly[i+1] + (NEWHOPE_3Q - poly[i+m+1]);
      
      poly[i] = (poly[i] + poly[i+m]) % NEWHOPE_Q;
      poly[i+1] = (poly[i+1] + poly[i+m+1]) % NEWHOPE_Q;
      
      poly[i+m] = montgomery_reduce((uint32_t)a * sub_re);
      poly[i+m+1] = montgomery_reduce((uint32_t)a * sub_img);        
    }
    index += window;
  }

}


void idgt(uint16_t *poly)
{

  int i, index, j, m, window;
  uint32_t a, mul_re, mul_img;

  window = (NEWHOPE_K2 >> 1);
  for(m = 2; m <= NEWHOPE_K2; m <<= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = invgj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        mul_re = montgomery_reduce((uint32_t)poly[i+m] * a);
        mul_img = montgomery_reduce((uint32_t)poly[i+m+1] * a);
        
        poly[i+m] = (poly[i] + (NEWHOPE_3Q - mul_re)) % NEWHOPE_Q;
        poly[i+m+1] = (poly[i+1] + (NEWHOPE_3Q - mul_img)) % NEWHOPE_Q;
        
        poly[i] = (poly[i] + mul_re) % NEWHOPE_Q;
        poly[i+1] = (poly[i+1] + mul_img) % NEWHOPE_Q;        
      }
      index += window;
    }
    window >>= 1; 

    m <<= 1;

    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = invgj[index];
      for(i = j; i < NEWHOPE_N; i += (m << 1)) 
      {
        mul_re = (montgomery_reduce((uint32_t)poly[i+m] * a));
        mul_img = (montgomery_reduce((uint32_t)poly[i+m+1] * a));
        
        poly[i+m] = (poly[i] + (NEWHOPE_3Q - mul_re));
        poly[i+m+1] = (poly[i+1] + (NEWHOPE_3Q - mul_img));
        
        poly[i] = (poly[i] + mul_re);
        poly[i+1] = (poly[i+1] + mul_img);        
      }
      index += window;
    }
    window >>= 1;
  }
}