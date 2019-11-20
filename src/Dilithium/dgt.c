#include "dgt.h"
#include "params.h"
#include "reduce.h"

static uint64_t gj[254] = {5754110, 2626307, 155746, 7823731, 8060855, 622214, 294099, 6645043, 7004294, 5994413, 2473040, 5761264, 8114715, 7733724, 3494836, 914272, 7186140, 7505333, 550081, 3586334, 3036846, 569503, 2421825, 5608021, 6632133, 303639, 4030378, 154128, 4948617, 5429939, 7888873, 5990337, 6290804, 5331993, 3965669, 3423912, 2884613, 6763448, 2279270, 4036898, 462476, 3233500, 2525328, 1115886, 6452587, 1993974, 3659106, 3388463, 3292687, 2015288, 8008212, 672143, 4202580, 3449530, 640503, 239478, 1901171, 4147766, 1432279, 3118556, 581040, 3032699, 4084600, 741917, 70815, 5524342, 8369167, 3386254, 3409235, 7083498, 6756926, 6784727, 5435651, 2005478, 7057124, 7586081, 6155651, 3176064, 4307002, 8346556, 2443186, 4484442, 1672797, 7554389, 2141617, 3014127, 8254185, 2772798, 3756217, 5094439, 2850584, 2927335, 5642630, 7677729, 1876861, 6669695, 2219673, 4974853, 3192319, 5896609, 2218864, 8229353, 3011725, 4663133, 8139083, 3916241, 2551932, 3375454, 4960675, 3652030, 6342907, 7018347, 6018145, 6334329, 5785896, 3992482, 2734871, 655207, 7122278, 6112161, 2479233, 1131217, 6165265, 6143324, 920237, 5405296, 7412702, 7925534, 5289613, 5799860, 2795146, 3656683, 4907177, 1818255, 3672408, 3597701, 2546731, 6611444, 7656307, 2640183, 8026561, 7235408, 6777715, 33683, 1366741, 2839977, 5328398, 670714, 5732921, 4482244, 4645763, 3509828, 4822283, 3283554, 2814359, 5439266, 7401908, 2259943, 277735, 5397567, 1072474, 6855062, 1659787, 1456518, 7350582, 8316604, 3568337, 4953036, 3688027, 7402230, 4702233, 869004, 6407080, 6369665, 2940159, 4507002, 4849299, 5153270, 2913025, 3207336, 5645159, 1195290, 2606445, 4775828, 8118105, 2336685, 3538690, 5850635, 6053679, 6577606, 4734974, 791048, 5353489, 1328943, 6644144, 1434429, 7723179, 5143979, 3367518, 7797031, 1863047, 1450437, 3206371, 1795007, 4039098, 1831568, 2194575, 2509366, 74927, 5918845, 5460841, 2730686, 6243867, 1683848, 5665147, 1506556, 5026908, 7336323, 3012619, 327407, 7510598, 180759, 330106, 6224023, 4610849, 165385, 1028434, 2712761, 3235261, 7395243, 3687926, 4339130, 5248672, 3873457, 1234205, 8331330, 6472232, 6245510, 1716497, 7173015, 6404393, 1879341, 3285498, 7309672, 7834770, 271497, 5346355, 3691555, 3036228, 6156504, 3946615, 4306529, 2274781, 7640848};
static uint64_t invgj[254] = {74927, 2461572, 2194575, 5871051, 5460841, 5649731, 2136550, 1683848, 3012619, 8053010, 869819, 180759, 5665147, 6873861, 3353509, 7336323, 5353489, 7051474, 4734974, 7589369, 6644144, 6945988, 657238, 5143979, 3206371, 6585410, 4341319, 1831568, 3367518, 583386, 6517370, 1450437, 4610849, 8215032, 330106, 2156394, 1028434, 5667656, 5145156, 7395243, 1234205, 49087, 1908185, 6245510, 3687926, 4041287, 3131745, 3873457, 3946615, 4073888, 6105636, 7640848, 5346355, 4688862, 5344189, 6156504, 3285498, 1070745, 545647, 271497, 1716497, 1207402, 1976024, 1879341, 1366741, 5540440, 6777715, 8346734, 5328398, 7709703, 2647496, 4482244, 2814359, 2941151, 978509, 2259943, 4645763, 4870589, 3558134, 3283554, 5289613, 2580557, 7412702, 454883, 2795146, 4723734, 3473240, 1818255, 7656307, 5740234, 353856, 7235408, 3672408, 4782716, 5833686, 6611444, 1072474, 1525355, 277735, 2982850, 1659787, 6923899, 1029835, 8316604, 4702233, 7511413, 1973337, 6369665, 3568337, 3427381, 4692390, 7402230, 3538690, 2529782, 2326738, 6577606, 2606445, 3604589, 262312, 2336685, 2913025, 5173081, 2735258, 1195290, 2940159, 3873415, 3531118, 5153270, 5368692, 4663133, 8139083, 4464176, 4960675, 4728387, 2551932, 5004963, 6503556, 6669695, 2219673, 3405564, 2218864, 151064, 3192319, 2483808, 2037510, 7018347, 6018145, 2046088, 2734871, 7725210, 5785896, 4387935, 920237, 2975121, 6165265, 2237093, 2479233, 7249200, 7122278, 2268256, 1623491, 6784727, 5435651, 6374939, 6155651, 5204353, 7057124, 794336, 4295817, 741917, 70815, 2856075, 3409235, 1296919, 8369167, 4994163, 4073415, 8346556, 2443186, 3895975, 2141617, 5366290, 1672797, 826028, 5642630, 702688, 2850584, 5453082, 3756217, 3285978, 8254185, 5607619, 8008212, 7708274, 4177837, 3449530, 3659106, 4991954, 5087730, 2015288, 640503, 8140939, 6479246, 4147766, 7799377, 3032699, 6948138, 3118556, 3965669, 4956505, 5495804, 6763448, 7888873, 2390080, 2089613, 5331993, 2279270, 4343519, 7917941, 3233500, 1927830, 1993974, 5855089, 1115886, 6632133, 8076778, 2421825, 2772396, 4030378, 8226289, 3431800, 5429939, 7186140, 875084, 3494836, 7466145, 550081, 4794083, 5343571, 569503, 5907377, 5761264, 8114715, 646693, 8086318, 6645043, 7004294, 2386004, 8060855, 7758203, 155746, 556686, 2626307, 2626307};

void dgt(uint32_t *poly)
{ // It is supposed that the input vector is already folded.
  unsigned int m, distance, jtwiddle;
  unsigned j1, j2;
  unsigned int j, k;
  uint64_t temp_re, temp_img;
  uint64_t a, b, c, d;
  uint64_t copy[N];

  for(j = 0; j < N; ++j)
    copy[j] = (uint64_t) poly[j];
  
  distance = 128;
  jtwiddle = 0;
  for(m = 1; m < 128; m <<= 1)
  {
    for(k = 0; k < m; ++k)
    {
      j1 = 2 * k * distance;
      j2 = j1 + distance - 1;
      for(j = j1; j <= j2; j = j+2)
      {
        a = montgomery_reduce(gj[jtwiddle] * copy[j + distance]);
        b = montgomery_reduce(gj[jtwiddle + 1] * copy[j + distance + 1]);
        c = montgomery_reduce(gj[jtwiddle] * copy[j + distance + 1]);
        d = montgomery_reduce(gj[jtwiddle + 1] * copy[j + distance]);

        temp_re = reduce32(a + (2*Q - b));
        temp_img = reduce32(c + d);

        copy[j + distance] = (copy[j] + (2*Q - temp_re));
        copy[j + distance + 1] = (copy[j + 1] + (2*Q - temp_img));
        
        copy[j] = (copy[j] + temp_re);
        copy[j + 1] = (copy[j + 1] + temp_img);
      }
      jtwiddle += 2;
    }
    distance >>= 1;
  }

  for(j = 0; j < N; ++j)
    poly[j] = (uint32_t) copy[j];
    
}

void idgt(uint32_t *poly) 
{
  unsigned int distance, j1, jtwiddle, m;
  unsigned int j, h;
  uint64_t temp_re, temp_img, sum_re, sum_img;
  uint64_t a, b, c, d;
  uint64_t copy[N];

  for(j = 0; j < N; ++j)
    copy[j] = (uint64_t) poly[j];

  distance = 2;
  h = 0;
  m = 128;
  //for(m = 128; m > 1; m >>= 1)
  //{
    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = (temp_re + copy[j + distance]);
        copy[j + 1] = (temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
    h += m;
    distance <<= 1;

    m >>= 1;

    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = reduce32(temp_re + copy[j + distance]);
        copy[j + 1] = reduce32(temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
    h += m;
    distance <<= 1;

    m >>= 1;

    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = reduce32(temp_re + copy[j + distance]);
        copy[j + 1] = reduce32(temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
    h += m;
    distance <<= 1;

    m >>= 1;

    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = reduce32(temp_re + copy[j + distance]);
        copy[j + 1] = reduce32(temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
    h += m;
    distance <<= 1;

    m >>= 1;

    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = reduce32(temp_re + copy[j + distance]);
        copy[j + 1] = reduce32(temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
    h += m;
    distance <<= 1;

    m >>= 1;

    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = reduce32(temp_re + copy[j + distance]);
        copy[j + 1] = reduce32(temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
    h += m;
    distance <<= 1;

    m >>= 1;

    for(j1 = 0; j1 < distance; j1 = j1 + 2)
    {
      jtwiddle = h;
      for(j = j1; j < N; j = j + 2*distance)
      {
        temp_re = copy[j];
        temp_img = copy[j + 1];
        
        copy[j] = (temp_re + copy[j + distance]);
        copy[j + 1] = (temp_img + copy[j + distance + 1]);
        
        sum_re = (temp_re + (2*Q - copy[j + distance]));
        sum_img = (temp_img + (2*Q - copy[j + distance + 1]));

        a = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_re);
        b = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_img);
        c = montgomery_reduce((uint64_t)invgj[jtwiddle] * sum_img);
        d = montgomery_reduce((uint64_t)invgj[jtwiddle + 1] * sum_re);

        copy[j + distance] = (a + (2*Q - b));
        copy[j + distance + 1] = (c + d);

        jtwiddle += 2;
      }
    }
  //}

  for(j = 0; j < N; ++j)
  {
    poly[j] = (uint32_t) (montgomery_reduce((uint64_t)copy[j] * 32764)); // invofkmodp * 2**32 mod p
  }

}