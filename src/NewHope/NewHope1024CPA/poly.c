#include "poly.h"
#include "ntt.h"
#include "dgt.h"
#include "reduce.h"
#include "fips202.h"

static uint16_t nthroots[NEWHOPE_N] = {4075, 0, 3993, 1027, 1000, 9181, 6531, 10539, 3457, 8728, 4113, 9689, 8814, 45, 6264, 12166, 418, 852, 29, 7430, 9893, 10207, 3759, 3876, 11633, 451, 7205, 11875, 2680, 1570, 2155, 3782, 6635, 407, 8528, 9945, 3736, 8218, 9331, 508, 2447, 2385, 3202, 1473, 2510, 5833, 9355, 10341, 8756, 9303, 9649, 3349, 7688, 8518, 9889, 3725, 10054, 4930, 11450, 8776, 10193, 140, 6572, 4805, 7192, 1446, 11128, 8709, 6752, 6100, 9153, 7118, 8345, 3416, 2485, 11325, 576, 8768, 7061, 2452, 10535, 4719, 1510, 5144, 3793, 2532, 93, 1905, 1719, 1486, 1521, 9507, 7083, 1822, 5275, 8673, 7775, 8694, 4574, 9414, 6449, 11659, 3981, 11572, 7616, 11296, 9248, 7590, 421, 2473, 10682, 11489, 7852, 12234, 6007, 9989, 4780, 11320, 2811, 8286, 6972, 2932, 3641, 4992, 11022, 8849, 9906, 11497, 12159, 2010, 10494, 1017, 9734, 11268, 4334, 9565, 3531, 9735, 3217, 7932, 6874, 3562, 10069, 8537, 2175, 9588, 11654, 2989, 11206, 4916, 3683, 102, 11997, 5799, 6051, 5663, 8338, 3609, 8862, 2420, 4277, 7150, 6237, 8654, 11328, 2998, 7564, 2988, 4596, 1342, 5836, 6448, 8412, 10657, 9800, 4194, 9351, 8606, 3363, 5316, 695, 5752, 6457, 356, 9627, 5978, 9525, 6182, 11171, 6320, 10810, 1867, 1108, 3245, 11619, 9614, 9110, 4704, 11743, 3202, 9026, 11456, 9187, 5388, 6257, 11449, 11759, 9813, 5735, 8529, 4607, 445, 10135, 6374, 2631, 12063, 328, 3394, 4837, 6026, 10908, 2365, 4435, 8236, 9743, 5156, 3926, 7369, 7207, 306, 1252, 2273, 4695, 2582, 5717, 3656, 8990, 3054, 8920, 4827, 4647, 6800, 3853, 5660, 9644, 4183, 10178, 7998, 11510, 3074, 9032, 355, 3426, 4325, 4945, 10400, 10397, 1455, 8844, 3274, 9587, 11001, 10162, 10548, 10149, 6211, 8485, 4455, 1145, 10002, 97, 3503, 3781, 12106, 3552, 2773, 3807, 11794, 179, 11147, 1959, 7165, 4372, 2035, 7966, 5611, 9350, 11592, 11159, 434, 8123, 6193, 9710, 7045, 11177, 12155, 1604, 12176, 8666, 7106, 416, 10097, 10164, 6223, 4911, 3599, 3302, 1640, 473, 3748, 8756, 9401, 4766, 7664, 10308, 852, 2618, 7902, 10209, 8635, 8723, 3323, 1567, 358, 6391, 3680, 353, 1830, 1480, 11264, 11781, 9426, 10061, 987, 816, 285, 10568, 2421, 5484, 11818, 5233, 11370, 4340, 5346, 12248, 8980, 11590, 731, 7035, 11653, 11523, 1382, 991, 9391, 5076, 6421, 4007, 9094, 5579, 3802, 1909, 3108, 1799, 10137, 12260, 10308, 3381, 8890, 3479, 9524, 11817, 5957, 8561, 7413, 7582, 6846, 9531, 11098, 9760, 4881, 7091, 10403, 7303, 11833, 3821, 11626, 2144, 6203, 10417, 9399, 10411, 4691, 6474, 993, 4468, 7285, 11933, 11826, 6517, 9903, 422, 2097, 5259, 2131, 5747, 9267, 4478, 3442, 11120, 4251, 9080, 1021, 9222, 2936, 2443, 7680, 9427, 47, 3185, 1291, 5187, 7882, 12089, 9483, 3352, 1413, 11101, 7154, 11065, 2748, 7312, 1573, 4190, 91, 11043, 732, 9234, 4004, 919, 6899, 9955, 3961, 10982, 4816, 10841, 8865, 3271, 8575, 1415, 4845, 3724, 4697, 2791, 10529, 9822, 10535, 4228, 12147, 8651, 9295, 5299, 8173, 9167, 8189, 5631, 11770, 178, 6173, 12123, 3666, 10162, 6931, 3654, 11349, 10850, 3492, 2532, 3430, 8797, 5323, 3735, 5331, 3411, 8270, 9108, 9301, 2425, 987, 11449, 1521, 5980, 5178, 10568, 7558, 2202, 7898, 2625, 2528, 7300, 10494, 10968, 9701, 10721, 174, 7310, 11042, 6675, 12258, 10294, 8202, 2691, 782, 7669, 4620, 8784, 1312, 7396, 4386, 11558, 1608, 8200, 11681, 7377, 9375, 11428, 382, 11147, 1632, 68, 6429, 9321, 4441, 11817, 11164, 10569, 3064, 8601, 10084, 9053, 5185, 3324, 10004, 763, 5951, 11950, 6111, 5866, 11165, 6562, 10555, 11674, 5585, 1865, 9615, 8572, 1146, 1066, 4571, 1132, 9226, 696, 11699, 8379, 8952, 8783, 9237, 458, 7682, 6821, 11705, 3773, 31, 9591, 8035, 11408, 9814, 4275, 10646, 3249, 722, 2875, 3197, 3658, 10148, 8545, 8619, 2965, 3473, 8164, 4464, 7128, 4185, 5384, 95, 4730, 11245, 7983, 1113, 4113, 8380, 8203, 11754, 1113, 7472, 10637, 7382, 592, 6220, 5062, 3165, 3110, 9742, 4395, 2311, 6387, 10257, 3393, 10454, 4327, 5001, 1529, 11199, 6393, 3488, 7249, 4260, 1708, 3658, 5797, 6156, 6071, 9776, 6549, 3751, 888, 10697, 9076, 232, 12119, 1227, 7745, 11823, 11241, 9885, 5924, 1940, 9182, 9163, 8045, 7125, 1971, 11090, 8087, 7796, 106, 9246, 7442, 10112, 4403, 4471, 5953, 3437, 7847, 2537, 4417, 1855, 3743, 5069, 4278, 8920, 8500, 1242, 8411, 6966, 2810, 381, 1038, 10198, 5558, 4419, 6356, 387, 8902, 6574, 8975, 9216, 12090, 7170, 7793, 11475, 97, 10497, 2145, 6812, 9209, 3660, 2287, 3172, 7479, 6459, 6723, 2469, 7220, 10320, 1820, 7333, 7562, 5892, 759, 8976, 9896, 9819, 10057, 11645, 5038, 12163, 4896, 12238, 11083, 4076, 10653, 7523, 3425, 9589, 7580, 3495, 6999, 8788, 11869, 4983, 12232, 1728, 10329, 10922, 10529, 11312, 10687, 4341, 7475, 12125, 4118, 2456, 3769, 8357, 11382, 4009, 10604, 12198, 6083, 4830, 6987, 2258, 9670, 10798, 7128, 4993, 12023, 7081, 3760, 10237, 6219, 7212, 2883, 11089, 8040, 1906, 10498, 8517, 1370, 3624, 2139, 11894, 754, 6675, 7205, 8323, 9875, 8943, 3499, 10525, 3074, 242, 5418, 1808, 6521, 10993, 7983, 5648, 4493, 776, 8812, 9625, 6604, 4517, 11916, 4655, 2430, 2293, 2147, 4726, 2885, 9973, 6553, 7541, 7090, 1212, 11274, 1349, 7221, 6817, 2354, 110, 8633, 499, 1129, 4453, 1133, 10806, 2645, 7207, 884, 6465, 8094, 6808, 2117, 11612, 2049, 4448, 10917, 10656, 4028, 7831, 4980, 4162, 1361, 9302, 6616, 2003, 6542, 2653, 992, 11284, 9612, 2293, 2372, 2688, 1104, 4055, 7198, 1441, 6500, 9164, 4895, 1064, 3886, 10485, 4175, 8108, 6930, 10140, 4039, 9285, 4368, 7411, 1001, 353, 6673, 1717, 6875, 2556, 6192, 3040, 4294, 8502, 4200, 9968, 2505, 7385, 9927, 3314, 11766, 2646, 4261, 12072, 1552, 6383, 6230, 367, 8865, 2174, 5506, 5622, 11213, 4493, 2747, 8661, 3401, 2837, 11149, 4204, 8686, 8422, 3189, 8596, 4909, 5822, 668, 983, 3130, 3710, 1363, 549, 5573, 4057, 4363, 4632, 2565, 382, 7004, 6060, 8448, 9386, 2360, 2065, 4657, 1231, 2516, 5854, 10203, 11690, 1575, 12265, 12234, 6821, 11083, 11685, 8093, 6360, 3203, 7762, 7455, 9173, 9378, 641, 4661, 5567, 824, 2854, 11198, 8719, 1812, 7771, 5951, 12239, 304, 8785, 7906, 8332, 4679, 11915, 819, 861, 7205, 10251, 2842, 9065, 11405, 6173, 7747, 1411, 7715, 1215, 5037, 270, 1297, 10940, 7202, 12250, 5933, 5643, 11412, 360, 5018, 6136, 2487, 2057, 4503, 2567, 5409, 8602, 4046, 11116, 1724, 4614, 3509, 5408, 4295, 10306, 10788, 6561, 11532, 6687, 3152, 12188, 6757, 1000, 9056, 4640, 3024, 10484, 5520, 8125, 11383, 4422, 6873, 10118, 5008, 6354, 10245, 11659, 11615, 10861, 7869, 2104};
static uint16_t invnthroots[NEWHOPE_N] = {1024, 0, 9609, 9042, 7099, 4515, 5662, 200, 2233, 6893, 5896, 10390, 8946, 5668, 3248, 6891, 617, 9973, 10771, 4573, 1156, 4822, 4848, 3828, 10827, 531, 4586, 1686, 9161, 9536, 10774, 8675, 12019, 9018, 7148, 9873, 4930, 1563, 11963, 7231, 11865, 8775, 408, 7373, 7298, 7060, 11783, 2137, 1872, 1604, 3307, 811, 11618, 8285, 3125, 1650, 6006, 3186, 10921, 6724, 7284, 1237, 11801, 1748, 7828, 2445, 5663, 3390, 5601, 10818, 8435, 2303, 2400, 3001, 7951, 953, 11603, 9077, 10476, 2685, 3142, 2526, 6099, 7740, 2100, 1459, 8383, 6276, 1945, 7505, 4414, 3557, 4395, 9649, 1152, 1866, 4174, 10473, 1655, 10167, 2357, 2334, 11481, 2679, 4165, 12256, 4056, 4389, 6242, 230, 11155, 511, 1888, 9435, 10515, 3979, 6255, 2772, 1972, 7486, 3191, 2141, 5218, 5604, 1281, 11391, 7121, 6725, 11292, 3491, 2098, 8966, 5538, 9797, 502, 6219, 6249, 7694, 6962, 4104, 841, 762, 10416, 7904, 8171, 11763, 685, 9514, 1901, 9639, 807, 4976, 9730, 9488, 1548, 2280, 202, 10486, 3607, 790, 7634, 11181, 653, 751, 9013, 9537, 4840, 837, 4064, 3776, 569, 2193, 10373, 1469, 2532, 4775, 4566, 1412, 1984, 3836, 6155, 3255, 537, 6683, 11373, 10749, 7835, 6791, 2168, 10343, 8197, 3883, 9137, 5549, 5071, 9009, 4650, 7878, 7698, 40, 7918, 3304, 5019, 7553, 9194, 5565, 10445, 4070, 9739, 5228, 7458, 5036, 626, 8847, 7009, 2391, 4587, 2516, 8982, 436, 3269, 8517, 6702, 7319, 567, 3301, 6643, 4744, 537, 6039, 10051, 6674, 4386, 9767, 4982, 5150, 11908, 6751, 11543, 2225, 763, 5783, 11528, 1995, 673, 84, 10938, 8195, 851, 7018, 6033, 1748, 11403, 11614, 6671, 4360, 10383, 4315, 9010, 55, 6824, 4961, 8444, 3205, 10205, 3576, 184, 8434, 4204, 11810, 6116, 10341, 10123, 9467, 2217, 3573, 1651, 9337, 4368, 5143, 4192, 5620, 4401, 8866, 5002, 1040, 7872, 2419, 545, 9127, 10029, 1543, 4171, 4232, 3079, 9553, 6596, 4418, 8291, 4149, 4286, 7459, 6710, 4643, 7566, 7495, 976, 3557, 2448, 1517, 6048, 8333, 6334, 3465, 7959, 8026, 11556, 11854, 12120, 6595, 4397, 1337, 8489, 2468, 4378, 3190, 9482, 2611, 7501, 11464, 8655, 11917, 4827, 4648, 12282, 4656, 2205, 5394, 12221, 2737, 36, 685, 3962, 9470, 6002, 10152, 9090, 8715, 2056, 668, 6290, 11990, 9724, 10480, 1829, 2463, 1955, 8720, 2468, 7618, 9272, 3103, 1114, 7986, 7070, 3097, 6594, 2431, 6184, 835, 10885, 5088, 6751, 7217, 8396, 8585, 2092, 5201, 2580, 10621, 5192, 1705, 4791, 11141, 10079, 3090, 2549, 4129, 1153, 5533, 2682, 5757, 4287, 7127, 10023, 8761, 11737, 7898, 8751, 8250, 4433, 3860, 4622, 11928, 3164, 11947, 5732, 11072, 2057, 3107, 11513, 11640, 11962, 2047, 11655, 1812, 7837, 9485, 8665, 3838, 2045, 6727, 10014, 4268, 1102, 496, 3297, 800, 8021, 2225, 956, 5838, 7729, 363, 8033, 10341, 6930, 10913, 5342, 7141, 4114, 4623, 4456, 3538, 6301, 2821, 2211, 8484, 5130, 8576, 8199, 6868, 7568, 5675, 10801, 9058, 3454, 7894, 12223, 9695, 11317, 3758, 419, 8944, 3742, 8830, 11845, 5180, 1794, 2012, 6437, 5919, 5462, 3497, 2278, 7347, 9498, 7751, 4796, 11210, 1608, 8306, 9288, 12046, 11368, 12084, 9189, 4429, 7528, 7311, 396, 3463, 4844, 1922, 8034, 5004, 10922, 3264, 7687, 6629, 6242, 7828, 4693, 10004, 4606, 352, 8839, 1779, 10674, 10916, 10758, 3806, 11731, 11731, 11620, 6278, 11841, 2552, 1488, 886, 10700, 6788, 3937, 10759, 1334, 10326, 137, 6308, 1546, 3110, 1855, 7384, 5886, 3415, 9525, 4393, 726, 8836, 1780, 5799, 8245, 7069, 8577, 3971, 2181, 7234, 2565, 4430, 7406, 10935, 4430, 4662, 8253, 3346, 11404, 8505, 8367, 4321, 10921, 8544, 334, 12219, 176, 9901, 944, 8572, 8533, 9711, 6816, 6320, 10458, 4474, 10746, 11078, 8035, 6706, 931, 6475, 6226, 9540, 4595, 4230, 2323, 10998, 6496, 10858, 651, 7245, 4432, 828, 1731, 1637, 7921, 4496, 10519, 6884, 3275, 2693, 700, 4421, 5910, 1139, 11798, 2689, 2623, 3196, 11766, 5412, 10033, 10092, 30, 6663, 6540, 252, 148, 5725, 4865, 5333, 6830, 6031, 9877, 5498, 8313, 6652, 9945, 7967, 3927, 5591, 9935, 7490, 6701, 5551, 1492, 3527, 8323, 8168, 3541, 8456, 9481, 4600, 7246, 11362, 9599, 6452, 4505, 8565, 11492, 1498, 8012, 2795, 3195, 7555, 557, 5391, 9000, 1922, 9830, 7235, 3395, 2531, 9065, 10897, 4984, 329, 10573, 5609, 1839, 9723, 5892, 8001, 11306, 10157, 3925, 10701, 7398, 99, 5950, 5877, 1779, 3315, 11364, 10321, 1461, 11696, 7245, 5404, 10319, 5163, 6682, 3415, 10898, 2301, 1780, 3657, 2245, 194, 44, 9595, 10472, 4655, 7695, 11832, 7394, 1482, 253, 878, 3346, 10761, 1663, 2774, 8260, 3224, 798, 7566, 3445, 2462, 4431, 10415, 7303, 11028, 11583, 2237, 8521, 8601, 6904, 7679, 3004, 10431, 5424, 3258, 6432, 8069, 5932, 11387, 9961, 8945, 3746, 7205, 8878, 6396, 1030, 1409, 632, 943, 172, 8009, 5660, 8592, 11471, 10690, 2075, 10739, 8784, 9442, 3902, 4656, 11464, 5804, 7362, 1743, 9097, 7881, 9834, 8505, 379, 5483, 2605, 6686, 3894, 7496, 4649, 3869, 1313, 4691, 7538, 3421, 12205, 11764, 9344, 9273, 8129, 8219, 10967, 609, 5403, 1854, 1795, 10334, 876, 1405, 8847, 4058, 11243, 4158, 1497, 10940, 9890, 1844, 2669, 4113, 10581, 682, 10209, 3967, 9370, 7446, 5688, 10974, 9134, 3455, 10848, 3398, 1273, 7209, 3218, 12223, 8434, 4922, 8247, 11427, 3453, 5400, 11734, 10861, 3117, 3133, 6061, 10659, 7699, 7165, 5510, 4707, 3997, 4028, 8696, 2742, 3865, 7781, 10489, 2507, 7992, 7403, 7490, 2711, 6551, 8782, 2901, 1667, 4738, 6444, 7601, 3418, 4602, 10528, 10010, 9770, 9318, 11695, 4044, 6691, 3564, 3028, 2434, 4440, 892, 8672, 6730, 7550, 11103, 6976, 10823, 7801, 4295, 10562, 7393, 4738, 9812, 9461, 3996, 6387, 6758, 6088, 8050, 4041, 1070, 10438, 223, 6948, 11991, 9731, 7862, 11408, 12141, 250, 340, 12152, 1832, 6049, 1149, 8506, 5363, 629, 6164, 2722, 6732, 2853, 7809, 12038, 9645, 8238, 12088, 5689, 2640, 8226, 1533, 8887, 4186, 7919, 4350, 1500, 10797, 9187, 9838, 6753, 5662, 2327, 2821, 10639, 514, 4530, 1522, 7420, 10856, 8181, 10646, 11563, 2406, 8778, 6872, 4464, 1354, 10018, 11157, 11035, 6979, 1831, 5194, 7125, 9251, 3070, 9405, 8679, 8078, 7312, 2428, 9229, 2136, 4582, 12083, 5717, 4326, 1124, 2851, 8231, 5569, 9993, 8867, 8884, 9140, 3321, 5535, 7690, 8962, 354, 11294, 8459, 8149, 2462, 7481, 6636, 2663, 9879, 3030, 6228, 8410, 6855, 194, 5484, 11073, 7282, 1911, 3807, 5042, 11255, 2799, 5128, 10663, 5750, 7583, 1748, 2930, 5379, 10576, 8386, 1624, 7882, 12030, 1392, 8260, 7775, 5904, 5736, 10129, 5246, 1910, 800, 11171, 6179, 10266, 6263, 1716, 11133, 12149, 7329};

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
void poly_sample_dgt(poly *r, const unsigned char *seed, unsigned char nonce)
{
#if NEWHOPE_K != 8
#error "poly_sample in poly.c only supports k=8"
#endif
  unsigned char buf[128], a, b, c, d;
  uint16_t x, y;
  int i, j, k;

  unsigned char extseed[NEWHOPE_SYMBYTES+2];

  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[NEWHOPE_SYMBYTES] = nonce;

  k = 0;
  for(i=0;i<NEWHOPE_N/64;i++) /* Generate noise in blocks of 64 coefficients */
  {
    extseed[NEWHOPE_SYMBYTES+1] = i;
    shake256(buf,128,extseed,NEWHOPE_SYMBYTES+2);
    for(j=0;j<64;j+=2)
    {
      a = buf[2*j];
      b = buf[2*j+1];
      c = buf[2*j+2];
      d = buf[2*j+3];
      x = hw(a) + NEWHOPE_Q - hw(b);
      y = hw(c) + NEWHOPE_Q - hw(d);

      r->coeffs[64*i+j] = montgomery_reduce((uint32_t)x * nthroots[k]) +
      (NEWHOPE_3Q - montgomery_reduce((uint32_t)y * nthroots[k+1]));
      r->coeffs[64*i+j+1] = montgomery_reduce((uint32_t)x * nthroots[k+1]) +
      montgomery_reduce((uint32_t)y * nthroots[k]);

      k += 2;
    }
  } 

  dgt((uint16_t *)r->coeffs);

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
  uint16_t t, s;

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
    t = montgomery_reduce(3186*a->coeffs[i]);             
    s = montgomery_reduce(3186*a->coeffs[i+1]);

    r->coeffs[i]   = montgomery_reduce((uint32_t)t * b->coeffs[i]) +
                    (NEWHOPE_3Q - montgomery_reduce((uint32_t)s * b->coeffs[i+1]));

    r->coeffs[i+1] = montgomery_reduce((uint32_t)t * b->coeffs[i+1]) + 
                    montgomery_reduce((uint32_t)s * b->coeffs[i]);
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

  int i;
  
  idgt((uint16_t *)r->coeffs);

  uint16_t a, b, c, d;

  for(i = 0; i < NEWHOPE_N; i += 2) 
  {
    a = montgomery_reduce((uint32_t)r->coeffs[i] * invnthroots[i]);
    b = montgomery_reduce((uint32_t)r->coeffs[i+1] * invnthroots[i+1]);
    c = montgomery_reduce((uint32_t)r->coeffs[i] * invnthroots[i+1]);
    d = montgomery_reduce((uint32_t)r->coeffs[i+1] * invnthroots[i]);

    r->coeffs[i] = a + NEWHOPE_3Q - b;
    r->coeffs[i+1] = c + d;
  }

}
