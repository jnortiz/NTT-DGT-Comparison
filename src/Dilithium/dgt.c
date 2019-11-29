#include "dgt.h"
#include "params.h"
#include "reduce.h"

static uint32_t gj[64] = {4193792, 2348700, 6271868, 7426187, 2091905, 3077325, 8360995, 4873154, 777960, 4519302, 1024112, 5582638, 5268920, 4464978, 6980856, 4558682, 2608894, 300467, 7830929, 531354, 6554070, 8284641, 8100412, 6779997, 466468, 4805995, 4794489, 3900724, 5495562, 1661693, 3859737, 6681150, 25847, 7841118, 5760665, 4499374, 8021166, 3530437, 6623180, 2140649, 237124, 5336701, 5654953, 6308525, 5700314, 5842901, 5102745, 4874723, 518909, 3539968, 7260833, 7568473, 6026966, 2706023, 4010497, 3699596, 876248, 2867647, 7300517, 5823537, 5260684, 3592148, 6262231, 6736599};

static uint32_t invgj[64] = {4193792, 1643818, 2118186, 4788269, 3119733, 2556880, 1079900, 5512770, 7504169, 4680821, 4369920, 5674394, 2353451, 811944, 1119584, 4840449, 7861508, 3505694, 3277672, 2537516, 2680103, 2071892, 2725464, 3043716, 8143293, 6239768, 1757237, 4849980, 359251, 3881043, 2619752, 539299, 8354570, 1699267, 4520680, 6718724, 2884855, 4479693, 3585928, 3574422, 7913949, 1600420, 280005, 95776, 1826347, 7849063, 549488, 8079950, 5771523, 3821735, 1399561, 3915439, 3111497, 2797779, 7356305, 3861115, 7602457, 3507263, 19422, 5303092, 6288512, 954230, 2108549, 6031717};

void dgt(uint32_t *poly)
{ 
  unsigned int i, index, j, m, window;
  uint32_t a, sub_re, sub_img;

  window = 1;
  for(m = 128; m >= 4; m >>= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];
      for(i = j; i < N; i += (m << 1)) 
      {
        sub_re = poly[i] + (2*Q - poly[i+m]);
        sub_img = poly[i+1] + (2*Q - poly[i+m+1]);
        
        poly[i] = reduce32(poly[i] + poly[i+m]);
        poly[i+1] = reduce32(poly[i+1] + poly[i+m+1]);
        
        poly[i+m] = montgomery_reduce((uint64_t)a * sub_re);
        poly[i+m+1] = montgomery_reduce((uint64_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;

    m >>= 1;

    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];
      for(i = j; i < N; i += (m << 1)) 
      {
        sub_re = poly[i] + (2*Q - poly[i+m]);
        sub_img = poly[i+1] + (2*Q - poly[i+m+1]);
        
        poly[i] = poly[i] + poly[i+m];
        poly[i+1] = poly[i+1] + poly[i+m+1];
        
        poly[i+m] = montgomery_reduce((uint64_t)a * sub_re);
        poly[i+m+1] = montgomery_reduce((uint64_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;

  }

  index = 0;
  for(j = 0; j < m; j += 2) 
  {
    a = gj[index];
    for(i = j; i < N; i += (m << 1)) 
    {
      sub_re = poly[i] + (2*Q - poly[i+m]);
      sub_img = poly[i+1] + (2*Q - poly[i+m+1]);
      
      poly[i] = poly[i] + poly[i+m];
      poly[i+1] = poly[i+1] + poly[i+m+1];
      
      poly[i+m] = montgomery_reduce((uint64_t)a * sub_re);
      poly[i+m+1] = montgomery_reduce((uint64_t)a * sub_img);        
    }
    index += window;
  }

}

void idgt(uint32_t *poly)
{
  unsigned int i, index, j, m, window;
  uint32_t a, mul_re, mul_img;

  window = 64;
  for(m = 2; m < N; m <<= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = invgj[index];
      for(i = j; i < N; i += (m << 1)) 
      {
        mul_re = montgomery_reduce((uint64_t)poly[i+m] * a);
        mul_img = montgomery_reduce((uint64_t)poly[i+m+1] * a);
        
        poly[i+m] = poly[i] + (2*Q - mul_re);
        poly[i+m+1] = poly[i+1] + (2*Q - mul_img);
        
        poly[i] = poly[i] + mul_re;
        poly[i+1] = poly[i+1] + mul_img;        
      }
      index += window;
    }
    window >>= 1;
  }
}
