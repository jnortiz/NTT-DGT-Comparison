#include <stdint.h>
#include "params.h"
#include "dgt.h"
#include "reduce.h"

/* Roots of unity in order needed by forward DGT */
static uint32_t gj[64] = {4193792, 25847, 2608894, 518909, 777960, 237124, 466468, 876248, 2091905, 8021166, 6554070, 6026966, 5268920, 5700314, 5495562, 5260684, 6271868, 5760665, 7830929, 7260833, 1024112, 5654953, 4794489, 7300517, 8360995, 6623180, 8100412, 4010497, 6980856, 5102745, 3859737, 6262231, 2348700, 7841118, 300467, 3539968, 4519302, 5336701, 4805995, 2867647, 3077325, 3530437, 8284641, 2706023, 4464978, 5842901, 1661693, 3592148, 7426187, 4499374, 531354, 7568473, 5582638, 6308525, 3900724, 5823537, 4873154, 2140649, 6779997, 3699596, 4558682, 4874723, 6681150, 6736599};

/* Roots of unity in order needed by inverse DGT */
static uint32_t invgj[64] = {4193792, 8354570, 7861508, 5771523, 7504169, 7913949, 8143293, 7602457, 3119733, 2884855, 2680103, 3111497, 2353451, 1826347, 359251, 6288512, 2118186, 4520680, 3277672, 1399561, 4369920, 280005, 1757237, 19422, 1079900, 3585928, 2725464, 7356305, 1119584, 549488, 2619752, 2108549, 1643818, 1699267, 3505694, 3821735, 4680821, 1600420, 6239768, 3507263, 2556880, 4479693, 2071892, 2797779, 811944, 7849063, 3881043, 954230, 4788269, 6718724, 2537516, 3915439, 5674394, 95776, 4849980, 5303092, 5512770, 3574422, 3043716, 3861115, 4840449, 8079950, 539299, 6031717};

/*************************************************
* Name:        dgt
*
* Description: Forward DGT, in-place. No modular reduction is performed after
*              additions or subtractions. If input coefficients are below 2*Q,
*              then output coefficients are below 18*Q.
*              Output vector is in bitreversed order.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
void dgt(uint32_t p[N]) {
  int m, distance;
  int j1, j2, j, k;
  uint32_t a, temp_re, temp_img;
  
  distance = 128;
  for(m = 1; m < 128; m <<= 1) {
    for(k = 0; k < m; ++k) {
      j1 = k*distance << 1;
      j2 = j1+distance-1;

      a = gj[k];
      for(j = j1; j <= j2; j = j+2) {
        temp_re = montgomery_reduce((uint64_t)a * p[j+distance]);
        temp_img = montgomery_reduce((uint64_t)a * p[j+distance+1]);

        p[j+distance] = p[j] + 2*Q - temp_re;
        p[j+distance+1] = p[j+1]+ 2*Q - temp_img;
        
        p[j] = p[j] + temp_re;
        p[j+1] = p[j+1] + temp_img;
      }
    }
    distance >>= 1;
  }
}

/*************************************************
* Name:        invdgt_tomont
*
* Description: Inverse DGT and multiplication by Montgomery factor 2^32.
*              In-place. No modular reductions after additions or
*              subtractions. Input coefficient need to be smaller than 2*Q.
*              Output coefficient are smaller than 2*Q.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
void invdgt_tomont(uint32_t p[N]) {
  int distance, j1, jtwiddle, j;
  uint32_t sub_re, sub_img;

  for(distance = 2; distance < N; distance <<= 1) {
    for(j1 = 0; j1 < distance; j1 += 2) {
      jtwiddle = 0;
      for(j = j1; j < N; j += distance << 1) {
        sub_re = p[j] + 256*Q - p[j+distance];
        sub_img = p[j+1] + 256*Q - p[j+distance+1];

        p[j] = p[j] + p[j+distance];
        p[j+1] = p[j+1] + p[j+distance+1];

        p[j+distance] = montgomery_reduce((uint64_t)invgj[jtwiddle] * sub_re);
        p[j+distance+1] = montgomery_reduce((uint64_t)invgj[jtwiddle++] * sub_img);
      }
    }   
  }
}
