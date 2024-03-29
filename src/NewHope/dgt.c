#include "inttypes.h"
#include "dgt.h"
#include "params.h"
#include "reduce.h"

#if (NEWHOPE_N == 512)
static uint16_t gj[NEWHOPE_N/4] = {4075, 5315, 4324, 4916, 10120, 11767, 7210, 9027, 1973, 5574, 11011, 2344, 8775, 1041, 1018, 6364, 3789, 147, 5456, 7840, 7540, 5537, 4789, 4467, 9606, 1190, 8471, 6118, 5445, 3860, 7753, 1050, 7083, 5529, 9090, 12233, 8724, 11635, 10587, 1987, 6427, 6136, 6874, 3643, 400, 1728, 4948, 6137, 10256, 3998, 10367, 8410, 1254, 11316, 5435, 1359, 6950, 5446, 6093, 3710, 11907, 316, 8301, 468, 6415, 677, 6234, 3336, 12237, 9115, 1323, 2766, 12138, 10162, 8332, 9450, 2505, 5906, 10710, 11858, 5241, 9369, 9162, 8120, 787, 8807, 1010, 6821, 2049, 7377, 10968, 192, 3445, 7509, 7591, 7232, 11286, 3532, 12048, 12231, 7280, 1956, 11404, 6008, 8851, 2844, 975, 4212, 5681, 8812, 12147, 11184, 3600, 3263, 7665, 6077, 421, 8209, 6068, 3602, 8076, 11785, 605, 9987, 9260, 5594, 6403, 7507};
static uint16_t invgj[NEWHOPE_N/4] = {4075, 6974, 7373, 7965, 3262, 5079, 522, 2169, 5925, 11271, 11248, 3514, 9945, 1278, 6715, 10316, 11239, 4536, 8429, 6844, 6171, 3818, 11099, 2683, 7822, 7500, 6752, 4749, 4449, 6833, 12142, 8500, 11821, 3988, 11973, 382, 8579, 6196, 6843, 5339, 10930, 6854, 973, 11035, 3879, 1922, 8291, 2033, 6152, 7341, 10561, 11889, 8646, 5415, 6153, 5862, 10302, 1702, 654, 3565, 56, 3199, 6760, 5206, 4782, 5886, 6695, 3029, 2302, 11684, 504, 4213, 8687, 6221, 4080, 11868, 6212, 4624, 9026, 8689, 1105, 142, 3477, 6608, 8077, 11314, 9445, 3438, 6281, 885, 10333, 5009, 58, 241, 8757, 1003, 5057, 4698, 4780, 8844, 12097, 1321, 4912, 10240, 5468, 11279, 3482, 11502, 4169, 3127, 2920, 7048, 431, 1579, 6383, 9784, 2839, 3957, 2127, 151, 9523, 10966, 3174, 52, 8953, 6055, 11612, 5874};
static const uint16_t invnthroots[NEWHOPE_N] = {1024, 0, 1273, 5080, 1139, 5910, 11136, 5533, 10877, 7723, 5007, 1216, 12205, 525, 8572, 944, 2565, 10480, 11527, 11448, 4977, 4211, 1030, 10880, 3110, 1546, 733, 11854, 33, 8124, 4370, 8103, 253, 11411, 1922, 4844, 8097, 5620, 1471, 6688, 6201, 5531, 5892, 4288, 7894, 3454, 11616, 84, 3271, 270, 3507, 5738, 1492, 8762, 6727, 2045, 2550, 5228, 10573, 11133, 8834, 3155, 3275, 9596, 3090, 10079, 1916, 1469, 3879, 6855, 7598, 10976, 334, 70, 668, 2056, 6040, 7694, 3038, 3070, 5084, 8543, 1334, 1963, 3465, 6334, 9932, 2334, 9649, 8226, 457, 4594, 7528, 4978, 9337, 1651, 4461, 2445, 2477, 9461, 6680, 1716, 7568, 6614, 5783, 763, 3128, 9536, 4297, 7403, 4799, 2354, 7837, 2804, 5565, 9194, 6179, 11171, 2919, 7446, 7793, 4368, 5192, 10584, 3776, 4064, 9879, 2663, 8395, 7496, 7968, 3922, 6002, 2137, 9797, 5538, 1831, 6979, 6357, 11387, 5501, 1589, 2448, 10772, 10473, 4174, 8238, 9645, 12245, 9595, 205, 921, 10123, 2822, 1237, 7284, 10562, 4295, 3224, 10897, 3713, 7159, 11908, 5538, 531, 10827, 7781, 3865, 2344, 7967, 10242, 327, 7918, 8985, 10129, 7043, 682, 10581, 7857, 828, 7088, 10197, 9013, 2752, 8149, 9827, 5483, 379, 4036, 3346, 11604, 12253, 11292, 8798, 1354, 2271, 3258, 5424, 448, 2552, 4794, 4723, 4395, 2640, 6732, 9436, 3657, 1780, 10681, 8306, 479, 8085, 6006, 9103, 11103, 5313, 7235, 9830, 5988, 2821, 2522, 7903, 1156, 7467, 3997, 8261, 5498, 9877, 10232, 3107, 4411, 7639, 4514, 4029, 9890, 10445, 10858, 6496, 5538, 7217, 1108, 4655, 11935, 3327, 9097, 4408, 10935, 7406, 10084, 5394, 898, 11008, 3511, 9883, 6904, 4610, 11731, 11731, 8003, 7459, 4784, 10344, 11660, 6926, 6682, 8874, 7751, 9498, 2084, 3576, 4004, 671, 3617, 11397, 557, 6898, 4623, 4114, 11752, 6039, 2316, 11672, 5124, 4590, 4865, 6956, 11947, 3164, 3152, 5549, 10665, 7882, 8131, 1046, 4595, 8059, 835, 6184, 12087, 10486, 3149, 3321, 6485, 825, 2181, 5055, 4648, 4827, 9098, 2141, 1433, 8181, 10052, 706, 10674, 1373, 4418, 6596, 10189, 1459, 10457, 6049, 6885, 5044, 5462, 8792, 4961, 6824, 10417, 1604, 8725, 3028, 9494, 4277, 6930, 1376, 3301, 567, 3343, 5668, 9172, 3133, 12037, 5749, 4433, 8429, 10343, 2168, 5379, 2930, 11413, 1405, 5814, 11358, 7070, 9192, 9488, 9730, 9993, 5569, 3505, 9442, 5220, 4044, 7501, 825, 2772, 6255, 4530, 514, 7858, 10415, 11937, 7683, 4171, 8057, 2526, 3142, 250, 12141, 925, 10321, 10277, 10495, 10383, 7974, 7060, 7298, 11695, 9318, 7784, 8565, 11926, 4560, 3269, 3772, 6893, 2233, 5400, 3453, 2256, 10092, 4391, 552, 11373, 1540, 10663, 6539, 1854, 5403, 1543, 11078, 9186, 3017, 1901, 2650, 4326, 11165, 10690, 11471, 11563, 8836, 9099, 7911, 1888, 2854, 5662, 9962, 7566, 798, 6047, 7828, 3162, 11744, 11603, 3212, 11991, 2558, 5877, 5950, 8547, 8830, 675, 886, 11865, 3514, 4602, 1761, 11362, 7246, 4268, 2225, 9773, 7702, 7099, 7774, 8434, 7367, 3196, 2623, 8002, 7127, 9034, 6134, 1034, 7247, 8129, 4070, 6320, 6816, 10334, 8720, 526, 4118, 7707, 10153, 172, 4280, 3415, 5886, 7892, 1337, 12059, 6047, 3102, 1492, 1663, 9515, 3264, 10922, 7287, 1040, 9288, 9889, 1851, 11219, 3925, 1588, 3758, 11317, 11438, 7018, 10726, 7359, 5845, 7551, 3541, 3833, 496, 1102, 11663, 8847};
#elif (NEWHOPE_N == 1024)
static uint16_t gj[NEWHOPE_N/4] = {4075, 5315, 4324, 4916, 10120, 11767, 7210, 9027, 1973, 5574, 11011, 2344, 8775, 1041, 1018, 6364, 3789, 147, 5456, 7840, 7540, 5537, 4789, 4467, 9606, 1190, 8471, 6118, 5445, 3860, 7753, 1050, 7083, 5529, 9090, 12233, 8724, 11635, 10587, 1987, 6427, 6136, 6874, 3643, 400, 1728, 4948, 6137, 10256, 3998, 10367, 8410, 1254, 11316, 5435, 1359, 6950, 5446, 6093, 3710, 11907, 316, 8301, 468, 6415, 677, 6234, 3336, 12237, 9115, 1323, 2766, 12138, 10162, 8332, 9450, 2505, 5906, 10710, 11858, 5241, 9369, 9162, 8120, 787, 8807, 1010, 6821, 2049, 7377, 10968, 192, 3445, 7509, 7591, 7232, 11286, 3532, 12048, 12231, 7280, 1956, 11404, 6008, 8851, 2844, 975, 4212, 5681, 8812, 12147, 11184, 3600, 3263, 7665, 6077, 421, 8209, 6068, 3602, 8076, 11785, 605, 9987, 9260, 5594, 6403, 7507, 5297, 6170, 3956, 1360, 11089, 7105, 9734, 6167, 10695, 1962, 5106, 6328, 9597, 168, 7991, 8960, 3728, 8240, 6299, 1159, 1146, 11341, 11964, 10885, 8527, 2919, 8273, 8212, 5766, 11637, 295, 6190, 8049, 8719, 11454, 6224, 8243, 709, 1319, 9139, 1958, 7967, 10211, 11177, 8210, 1058, 11848, 11367, 6507, 1566, 2948, 9786, 11606, 9830, 8633, 12225, 10542, 9166, 9235, 5486, 3834, 5257, 7856, 5919, 10314, 3757, 9364, 11942, 7535, 10431, 426, 3315, 2738, 6421, 2655, 6554, 723, 174, 1693, 9280, 350, 1512, 10474, 6906, 9087, 7796, 5369, 2057, 11026, 12240, 6374, 1483, 1583, 6347, 2500, 10800, 6142, 2447, 3963, 11713, 1954, 2051, 1805, 2882, 9928, 10446, 9259, 4115, 9381, 218, 8760, 3434, 156, 9522, 8320, 3991, 5876, 2281, 10258, 6956, 4774, 6860, 4737, 1293, 11871, 8517, 6381, 11836};
static uint16_t invgj[NEWHOPE_N/4] = {4075, 6974, 7373, 7965, 3262, 5079, 522, 2169, 5925, 11271, 11248, 3514, 9945, 1278, 6715, 10316, 11239, 4536, 8429, 6844, 6171, 3818, 11099, 2683, 7822, 7500, 6752, 4749, 4449, 6833, 12142, 8500, 11821, 3988, 11973, 382, 8579, 6196, 6843, 5339, 10930, 6854, 973, 11035, 3879, 1922, 8291, 2033, 6152, 7341, 10561, 11889, 8646, 5415, 6153, 5862, 10302, 1702, 654, 3565, 56, 3199, 6760, 5206, 4782, 5886, 6695, 3029, 2302, 11684, 504, 4213, 8687, 6221, 4080, 11868, 6212, 4624, 9026, 8689, 1105, 142, 3477, 6608, 8077, 11314, 9445, 3438, 6281, 885, 10333, 5009, 58, 241, 8757, 1003, 5057, 4698, 4780, 8844, 12097, 1321, 4912, 10240, 5468, 11279, 3482, 11502, 4169, 3127, 2920, 7048, 431, 1579, 6383, 9784, 2839, 3957, 2127, 151, 9523, 10966, 3174, 52, 8953, 6055, 11612, 5874, 453, 5908, 3772, 418, 10996, 7552, 5429, 7515, 5333, 2031, 10008, 6413, 8298, 3969, 2767, 12133, 8855, 3529, 12071, 2908, 8174, 3030, 1843, 2361, 9407, 10484, 10238, 10335, 576, 8326, 9842, 6147, 1489, 9789, 5942, 10706, 10806, 5915, 49, 1263, 10232, 6920, 4493, 3202, 5383, 1815, 10777, 11939, 3009, 10596, 12115, 11566, 5735, 9634, 5868, 9551, 8974, 11863, 1858, 4754, 347, 2925, 8532, 1975, 6370, 4433, 7032, 8455, 6803, 3054, 3123, 1747, 64, 3656, 2459, 683, 2503, 9341, 10723, 5782, 922, 441, 11231, 4079, 1112, 2078, 4322, 10331, 3150, 10970, 11580, 4046, 6065, 835, 3570, 4240, 6099, 11994, 652, 6523, 4077, 4016, 9370, 3762, 1404, 325, 948, 11143, 11130, 5990, 4049, 8561, 3329, 4298, 12121, 2692, 5961, 7183, 10327, 1594, 6122, 2555, 5184, 1200, 10929, 8333, 6119, 6992};
static const uint16_t invnthroots[NEWHOPE_N] = {512, 0, 1726, 591, 1224, 8497, 4663, 4951, 4028, 303, 772, 6295, 5345, 7951, 464, 4915, 10412, 10720, 1028, 3129, 3670, 2278, 4510, 4283, 683, 3995, 9476, 3143, 9523, 3532, 4072, 9722, 2265, 12032, 5700, 380, 11384, 1620, 4920, 10153, 2346, 4815, 6695, 9147, 10547, 8921, 7520, 6653, 4302, 4082, 2350, 11478, 6338, 4395, 5756, 11187, 7763, 625, 8671, 1971, 6866, 260, 572, 10194, 10275, 8143, 6205, 12069, 1730, 9229, 7636, 11507, 10275, 5986, 10927, 10318, 1500, 11177, 3652, 575, 5335, 7738, 11069, 11386, 1898, 12013, 7677, 8221, 11326, 8325, 10306, 5231, 11927, 344, 3471, 3190, 9747, 10416, 9842, 6497, 3093, 1484, 10970, 9269, 5136, 3298, 402, 629, 11358, 9401, 10996, 6015, 3, 11128, 10543, 722, 9650, 6227, 10006, 10433, 4671, 3517, 2878, 6711, 7836, 7800, 515, 3369, 6608, 3623, 3717, 6736, 4611, 12141, 4050, 10464, 9351, 2747, 774, 5005, 9481, 8549, 6660, 10649, 5527, 1871, 12224, 11271, 1532, 6358, 4018, 1120, 9276, 1165, 11372, 7482, 9682, 11467, 11938, 10085, 3617, 5054, 1079, 1555, 8251, 12176, 11600, 7055, 4871, 1410, 2024, 271, 3078, 4930, 10302, 2747, 5084, 9343, 5802, 5730, 1088, 6380, 4464, 4906, 2233, 4890, 155, 4397, 6800, 7441, 10470, 1875, 8001, 2565, 7426, 865, 8362, 1599, 9566, 7554, 7255, 8905, 10842, 1359, 1690, 2311, 8334, 5325, 8867, 3185, 989, 7651, 815, 4348, 3153, 21, 548, 134, 5581, 9917, 8447, 12007, 2495, 231, 5240, 7427, 2107, 11068, 4936, 10629, 8091, 5094, 2578, 4790, 8106, 2946, 5156, 7335, 1642, 11909, 6358, 8727, 727, 5364, 6868, 10065, 3254, 7823, 8396, 3726, 1969, 7851, 6192, 873, 6209, 5530, 3987, 11336, 7220, 7174, 11498, 4868, 12020, 10145, 1712, 6359, 9016, 11844, 881, 8859, 4817, 711, 4929, 6781, 10829, 5003, 10737, 11388, 7065, 5300, 5362, 2620, 8331, 7581, 8797, 1882, 1369, 7306, 11735, 8472, 10339, 10399, 8041, 6691, 158, 11076, 4453, 5406, 4451, 8411, 5944, 8362, 3617, 7688, 3175, 4770, 2413, 7336, 5475, 7090, 768, 7731, 4346, 9427, 1105, 1465, 1866, 4366, 3699, 3191, 4922, 3914, 4634, 7691, 7344, 5546, 11203, 11908, 11780, 9175, 6871, 2716, 11243, 4006, 7547, 9753, 4511, 7723, 1756, 6784, 11722, 7135, 6054, 7821, 988, 3399, 9416, 2854, 5646, 8681, 12157, 3523, 2521, 11772, 5423, 11896, 5462, 6015, 1413, 3252, 6146, 2537, 1925, 6142, 3528, 7206, 11496, 10411, 5735, 8379, 5125, 6760, 11826, 11927, 10666, 7774, 7541, 6404, 2044, 5259, 12269, 8231, 8121, 10683, 11373, 3120, 9632, 3604, 180, 3490, 5171, 1028, 164, 5090, 966, 2947, 2759, 2889, 4934, 3155, 10903, 2794, 9816, 8059, 9315, 11184, 5195, 4363, 4676, 375, 10746, 10934, 2819, 4860, 1061, 1468, 11435, 341, 3232, 2948, 1252, 5412, 3446, 776, 10049, 10578, 9761, 10896, 7703, 4627, 921, 4832, 7176, 2207, 4002, 506, 90, 3013, 2946, 10357, 7937, 1736, 6175, 4285, 11670, 10268, 10975, 10786, 6018, 11850, 8769, 9787, 5759, 7184, 3009, 12024, 4185, 3324, 8495, 699, 8317, 10362, 9164, 776, 5803, 1577, 3845, 4646, 7695, 3227, 1934, 9512, 7541, 6600, 4794, 9184, 6040, 1787, 9970, 11796, 11161, 5046, 4822, 9660, 6257, 3904, 9722, 5814, 5277, 10996, 9849, 9101, 2463, 5254, 10495, 11218, 1213, 9680, 10821, 5158, 9115, 8397, 10114, 8226, 70, 378, 4072, 443, 4119, 9986, 5719, 10264, 279, 279, 6931, 6375, 10990, 3017, 11364, 6727, 4070, 1664, 4407, 4307, 11685, 3789, 3071, 1275, 10105, 7963, 9680, 11602, 6519, 841, 3844, 10960, 7796, 773, 3571, 2340, 7849, 7666, 6858, 11789, 10063, 10887, 3091, 1777, 10025, 5070, 4229, 749, 3719, 1622, 488, 616, 358, 3984, 7265, 5539, 1560, 9177, 7891, 5327, 3723, 10193, 2825, 1720, 8162, 6971, 3939, 5079, 6624, 2491, 9974, 10427, 3946, 5476, 1989, 11158, 450, 8444, 2715, 5031, 249, 9173, 2084, 944, 9584, 7220, 5061, 9864, 4187, 10028, 8924, 9524, 5361, 2564, 5080, 6599, 10108, 7516, 9174, 2922, 8400, 9927, 10930, 7014, 6428, 7795, 4727, 3023, 3589, 1174, 2007, 2748, 6330, 6180, 8445, 2842, 8291, 10904, 11859, 1085, 9371, 11898, 9312, 11731, 10866, 12036, 8144, 3049, 2597, 6670, 6336, 2081, 4724, 3624, 989, 8549, 6643, 4399, 12108, 4856, 11402, 2600, 10378, 7957, 7727, 880, 4631, 1877, 2140, 10545, 6347, 8784, 11401, 12240, 8992, 1258, 7115, 11260, 11612, 5152, 6124, 11474, 5240, 10802, 7452, 4732, 4202, 1992, 12154, 4509, 4445, 9740, 7918, 9667, 6533, 3805, 2198, 7671, 10196, 6219, 1391, 7628, 5533, 10975, 7784, 637, 9448, 5804, 6045, 10502, 10272, 6114, 9089, 10290, 4553, 8817, 3275, 4904, 9556, 5983, 562, 10126, 719, 5166, 12206, 5836, 8681, 8561, 589, 9574, 11384, 11665, 4534, 11613, 5696, 6363, 72, 1671, 11971, 3172, 643, 9798, 5451, 6866, 5170, 4836, 7238, 5997, 9389, 921, 4882, 8158, 7953, 446, 7958, 2931, 12139, 12058, 25, 9705, 3307, 4207, 5308, 166, 3541, 11487, 3166, 4480, 10014, 6420, 7338, 1855, 11963, 7235, 6703, 4812, 5185, 12222, 11916, 3719, 9787, 5410, 346, 8389, 3748, 10342, 11498, 403, 5677, 6566, 4358, 9702, 564, 3606, 11228, 10503, 405, 10972, 9006, 11269, 1391, 5949, 6679, 6131, 8790, 7184, 3578, 7242, 11623, 9678, 4201, 11167, 2112, 3107, 11526, 2303, 794, 8107, 8512, 4052, 5104, 10620, 11956, 10394, 3297, 644, 8787, 9649, 6315, 12068, 2414, 1960, 979, 11650, 5310, 12153, 12009, 7543, 12071, 1055, 6424, 4937, 1916, 8169, 2382, 3972, 2485, 9755, 9382, 5679, 8487, 4724, 3812, 6496, 1704, 4121, 7300, 1122, 416, 3184, 8864, 1997, 12287, 7603, 1802, 11899, 7917, 9166, 7948, 3003, 9814, 10129, 8599, 7647, 6312, 4334, 8643, 6631, 6481, 3936, 8136, 804, 12194, 3269, 3899, 5414, 366, 1798, 6599, 10180, 9199, 5812, 732, 3497, 9928, 4521, 8064, 7355, 10750, 1604, 401, 6439, 2632, 5632, 5828, 8391, 9601, 9332, 9641, 130, 11997, 1510, 1545, 1105, 5877, 12, 6693, 3680, 7034, 3738, 820, 2623, 4945, 10989, 10250, 6534, 8146, 5968, 7227, 9984, 7990, 3860, 4622, 6826, 8134, 744, 5920, 10121, 5370, 677, 1168, 10881, 9403, 58, 1581, 6417, 3131, 8072, 8486, 11504, 2175, 6477, 2952, 9968, 5358, 5712, 11757, 8111, 197, 3927, 6188, 5089, 1329, 8505, 8320, 3011, 1074, 6073, 2635, 4911, 2230, 7140, 1004, 5858, 223, 11306, 8919, 1576, 8422, 3103, 5536, 7989, 872, 10520, 10574, 6324, 1656, 4866, 3374, 2858, 6971, 11273, 4271, 9878, 5492, 5467, 5987, 10416, 1007, 4077, 3321, 10898, 152, 8625, 7814, 2993, 11630, 2404, 908, 2087, 1444, 1027, 11099, 6137, 7602, 9618, 1900, 2555, 5304, 7014, 9256, 2213, 11631, 911, 5875, 7976, 4406, 4939, 9560, 1286, 3333, 11626, 10105, 10997, 7858, 10789, 9788, 6174, 5756, 2733};
#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif

#if(NEWHOPE_N == 512)

void dgt(uint16_t *x)
{
  int m, distance;
  int j1, j2, j, k;
  uint16_t a, temp_re, temp_img;
  
  distance = NEWHOPE_N/2;
  for(m = 1; m < NEWHOPE_N/2; m <<= 1)
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


void idgt(uint16_t *x)
{
  int distance, i, j1, j2, k, j, m;
  uint16_t a, a_re, a_img, sub_re, sub_img;
  uint16_t p[NEWHOPE_N];

  m = NEWHOPE_N/4;
  for(distance = 2; distance < NEWHOPE_N/2; distance <<= 1)
  {
    for(k = 0; k < m; ++k)
    {
      j1 = 2 * k * distance;
      j2 = j1 + distance - 1;

      a = invgj[k];
      for(j = j1; j <= j2; j += 2)
      {
        sub_re = x[j] + 3*NEWHOPE_Q - x[j+distance];
        sub_img = x[j+1] +  3*NEWHOPE_Q - x[j+distance+1];
        
        x[j] = x[j] + x[j+distance];
        x[j+1] = x[j+1] + x[j+distance+1];

        x[j+distance] = montgomery_reduce((uint32_t)a * sub_re);
        x[j+distance+1] = montgomery_reduce((uint32_t)a * sub_img);
      }
    }
    m >>= 1;
    distance <<= 1;
    if(distance < NEWHOPE_N/2) {
      for(k = 0; k < m; ++k)
      {
        j1 = 2 * k * distance;
        j2 = j1 + distance - 1;

        a = invgj[k];
        for(j = j1; j <= j2; j += 2)
        {
          sub_re = x[j] + 3*NEWHOPE_Q - x[j+distance];
          sub_img = x[j+1] +  3*NEWHOPE_Q - x[j+distance+1];
          
          x[j] = (x[j] + x[j+distance]) % NEWHOPE_Q;
          x[j+1] = (x[j+1] + x[j+distance+1]) % NEWHOPE_Q;

          x[j+distance] = montgomery_reduce((uint32_t)a * sub_re);
          x[j+distance+1] = montgomery_reduce((uint32_t)a * sub_img);
        }
      }
      m >>= 1;      
    }
  }

  for(j = 0; j < NEWHOPE_N; ++j)
    p[j] = x[j];

  i = 0;
  for(j = 0; j < NEWHOPE_N/2; j += 2)
  {
    sub_re = p[j] + 3*NEWHOPE_Q - p[j+NEWHOPE_N/2];
    sub_img = p[j+1] +  3*NEWHOPE_Q - p[j+NEWHOPE_N/2+1];
    
    a_re = p[j] + p[j+NEWHOPE_N/2];
    a_img = p[j+1] + p[j+NEWHOPE_N/2+1];

    x[i] = montgomery_reduce((uint32_t)a_re * invnthroots[j]) + 3*NEWHOPE_Q - montgomery_reduce((uint32_t)a_img * invnthroots[j+1]);
    x[i+NEWHOPE_N/2] = montgomery_reduce((uint32_t)a_re * invnthroots[j+1] + (uint32_t)a_img * invnthroots[j]);    
    x[i+NEWHOPE_N/4] = montgomery_reduce((uint32_t)sub_re * invnthroots[j+NEWHOPE_N/2]) + 3*NEWHOPE_Q - montgomery_reduce((uint32_t)sub_img * invnthroots[j+NEWHOPE_N/2+1]);
    x[i+NEWHOPE_N/2+NEWHOPE_N/4] = montgomery_reduce((uint32_t)sub_re * invnthroots[j+NEWHOPE_N/2+1] + (uint32_t)sub_img * invnthroots[j+NEWHOPE_N/2]);

    i++;
  }
}

#elif (NEWHOPE_N == 1024)

void dgt(uint16_t *x)
{
  int m, distance;
  int j1, j2, j, k;
  uint16_t a, temp_re, temp_img;
  
  distance = NEWHOPE_N/2;
  for(m = 1; m < NEWHOPE_N/2; m <<= 1)
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

void idgt(uint16_t *x)
{
  int distance, i, j1, j2, k, j, m;
  uint16_t a, a_re, a_img, sub_re, sub_img;
  uint16_t p[NEWHOPE_N];

  m = NEWHOPE_N/4;
  for(distance = 2; distance < NEWHOPE_N/2; distance <<= 1)
  {
    for(k = 0; k < m; ++k)
    {
      j1 = 2 * k * distance;
      j2 = j1 + distance - 1;

      a = invgj[k];
      for(j = j1; j <= j2; j += 2)
      {
        sub_re = x[j] + 3*NEWHOPE_Q - x[j+distance];
        sub_img = x[j+1] +  3*NEWHOPE_Q - x[j+distance+1];
        
        x[j] = x[j] + x[j+distance];
        x[j+1] = x[j+1] + x[j+distance+1];

        x[j+distance] = montgomery_reduce((uint32_t)a * sub_re);
        x[j+distance+1] = montgomery_reduce((uint32_t)a * sub_img);
      }
    }
    m >>= 1;
    distance <<= 1;
    if(distance < NEWHOPE_N/2)     {
      for(k = 0; k < m; ++k)
      {
        j1 = 2 * k * distance;
        j2 = j1 + distance - 1;

        a = invgj[k];
        for(j = j1; j <= j2; j += 2)
        {
          sub_re = x[j] + 3*NEWHOPE_Q - x[j+distance];
          sub_img = x[j+1] +  3*NEWHOPE_Q - x[j+distance+1];
          
          x[j] = (x[j] + x[j+distance]) % NEWHOPE_Q;
          x[j+1] = (x[j+1] + x[j+distance+1]) % NEWHOPE_Q;

          x[j+distance] = montgomery_reduce((uint32_t)a * sub_re);
          x[j+distance+1] = montgomery_reduce((uint32_t)a * sub_img);
        }
      }
      m >>= 1;      
    }
  }

  for(j = 0; j < NEWHOPE_N; ++j)
    p[j] = x[j];

  i = 0;
  for(j = 0; j < NEWHOPE_N/2; j += 2)
  {
    sub_re = p[j] + 3*NEWHOPE_Q - p[j+NEWHOPE_N/2];
    sub_img = p[j+1] +  3*NEWHOPE_Q - p[j+NEWHOPE_N/2+1];
    
    a_re = p[j] + p[j+NEWHOPE_N/2];
    a_img = p[j+1] + p[j+NEWHOPE_N/2+1];

    x[i] = montgomery_reduce((uint32_t)a_re * invnthroots[j]) + (3*NEWHOPE_Q - 
           montgomery_reduce((uint32_t)a_img * invnthroots[j+1]));
    x[i+NEWHOPE_N/2] = montgomery_reduce((uint32_t)a_re * invnthroots[j+1] + (uint32_t)a_img * invnthroots[j]);    
    x[i+NEWHOPE_N/4] = montgomery_reduce((uint32_t)sub_re * invnthroots[j+NEWHOPE_N/2]) + (3*NEWHOPE_Q - 
                       montgomery_reduce((uint32_t)sub_img * invnthroots[j+NEWHOPE_N/2+1]));
    x[i+NEWHOPE_N/2+NEWHOPE_N/4] = montgomery_reduce((uint32_t)sub_re * invnthroots[j+NEWHOPE_N/2+1] + (uint32_t)sub_img * invnthroots[j+NEWHOPE_N/2]);
    i++;
  }
}

#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif
