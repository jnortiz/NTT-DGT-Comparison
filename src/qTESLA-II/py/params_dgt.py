from gauss import GaussianInteger

# DGT_2
invofkmodp = 8273665
gj_powers = [[1, 1, 1, 1, 1, 1], [4704222, 5240745, 684268, 5750773, 3585646, 8404992], [5240745, 684268, 5750773, 3585646, 8404992, 1], [1597874, 8191273, 3506438, 6322675, 4819347, 8404992], [684268, 5750773, 3585646, 8404992, 1, 1], [4360356, 5121240, 7690526, 2654220, 3585646, 8404992], [8191273, 3506438, 6322675, 4819347, 8404992, 1], [2126834, 6920830, 1270094, 2082318, 4819347, 8404992], [5750773, 3585646, 8404992, 1, 1, 1], [5639303, 1651513, 7720725, 5750773, 3585646, 8404992], [5121240, 7690526, 2654220, 3585646, 8404992, 1], [8314555, 973655, 4898555, 6322675, 4819347, 8404992], [3506438, 6322675, 4819347, 8404992, 1, 1], [3463953, 2379409, 714467, 2654220, 3585646, 8404992], [6920830, 1270094, 2082318, 4819347, 8404992, 1], [2134075, 5433596, 7134899, 2082318, 4819347, 8404992], [3585646, 8404992, 1, 1, 1, 1], [5330453, 3164248, 684268, 5750773, 3585646, 8404992], [1651513, 7720725, 5750773, 3585646, 8404992, 1], [4153273, 213720, 3506438, 6322675, 4819347, 8404992], [7690526, 2654220, 3585646, 8404992, 1, 1], [2436145, 3283753, 7690526, 2654220, 3585646, 8404992], [973655, 4898555, 6322675, 4819347, 8404992, 1], [5146046, 1484163, 1270094, 2082318, 4819347, 8404992], [6322675, 4819347, 8404992, 1, 1, 1], [5400177, 6753480, 7720725, 5750773, 3585646, 8404992], [2379409, 714467, 2654220, 3585646, 8404992, 1], [2786978, 7431338, 4898555, 6322675, 4819347, 8404992], [1270094, 2082318, 4819347, 8404992, 1, 1], [5597909, 6025584, 714467, 2654220, 3585646, 8404992], [5433596, 7134899, 2082318, 4819347, 8404992, 1], [5785355, 2971397, 7134899, 2082318, 4819347, 8404992], [8404992, 1, 1, 1, 1, 1], [3700771, 5240745, 684268, 5750773, 3585646, 8404992], [3164248, 684268, 5750773, 3585646, 8404992, 1], [6807119, 8191273, 3506438, 6322675, 4819347, 8404992], [7720725, 5750773, 3585646, 8404992, 1, 1], [4044637, 5121240, 7690526, 2654220, 3585646, 8404992], [213720, 3506438, 6322675, 4819347, 8404992, 1], [6278159, 6920830, 1270094, 2082318, 4819347, 8404992], [2654220, 3585646, 8404992, 1, 1, 1], [2765690, 1651513, 7720725, 5750773, 3585646, 8404992], [3283753, 7690526, 2654220, 3585646, 8404992, 1], [90438, 973655, 4898555, 6322675, 4819347, 8404992], [4898555, 6322675, 4819347, 8404992, 1, 1], [4941040, 2379409, 714467, 2654220, 3585646, 8404992], [1484163, 1270094, 2082318, 4819347, 8404992, 1], [6270918, 5433596, 7134899, 2082318, 4819347, 8404992], [4819347, 8404992, 1, 1, 1, 1], [3074540, 3164248, 684268, 5750773, 3585646, 8404992], [6753480, 7720725, 5750773, 3585646, 8404992, 1], [4251720, 213720, 3506438, 6322675, 4819347, 8404992], [714467, 2654220, 3585646, 8404992, 1, 1], [5968848, 3283753, 7690526, 2654220, 3585646, 8404992], [7431338, 4898555, 6322675, 4819347, 8404992, 1], [3258947, 1484163, 1270094, 2082318, 4819347, 8404992], [2082318, 4819347, 8404992, 1, 1, 1], [3004816, 6753480, 7720725, 5750773, 3585646, 8404992], [6025584, 714467, 2654220, 3585646, 8404992, 1], [5618015, 7431338, 4898555, 6322675, 4819347, 8404992], [7134899, 2082318, 4819347, 8404992, 1, 1], [2807084, 6025584, 714467, 2654220, 3585646, 8404992], [2971397, 7134899, 2082318, 4819347, 8404992, 1], [2619638, 2971397, 7134899, 2082318, 4819347, 8404992]]
invgj_powers = [[1, 1, 1, 1, 1, 1], [8404992, 4819347, 2082318, 7134899, 2971397, 2619638], [1, 8404992, 4819347, 2082318, 7134899, 2971397], [8404992, 3585646, 2654220, 714467, 6025584, 2807084], [1, 1, 8404992, 4819347, 2082318, 7134899], [8404992, 4819347, 6322675, 4898555, 7431338, 5618015], [1, 8404992, 3585646, 2654220, 714467, 6025584], [8404992, 3585646, 5750773, 7720725, 6753480, 3004816], [1, 1, 1, 8404992, 4819347, 2082318], [8404992, 4819347, 2082318, 1270094, 1484163, 3258947], [1, 8404992, 4819347, 6322675, 4898555, 7431338], [8404992, 3585646, 2654220, 7690526, 3283753, 5968848], [1, 1, 8404992, 3585646, 2654220, 714467], [8404992, 4819347, 6322675, 3506438, 213720, 4251720], [1, 8404992, 3585646, 5750773, 7720725, 6753480], [8404992, 3585646, 5750773, 684268, 3164248, 3074540], [1, 1, 1, 1, 8404992, 4819347], [8404992, 4819347, 2082318, 7134899, 5433596, 6270918], [1, 8404992, 4819347, 2082318, 1270094, 1484163], [8404992, 3585646, 2654220, 714467, 2379409, 4941040], [1, 1, 8404992, 4819347, 6322675, 4898555], [8404992, 4819347, 6322675, 4898555, 973655, 90438], [1, 8404992, 3585646, 2654220, 7690526, 3283753], [8404992, 3585646, 5750773, 7720725, 1651513, 2765690], [1, 1, 1, 8404992, 3585646, 2654220], [8404992, 4819347, 2082318, 1270094, 6920830, 6278159], [1, 8404992, 4819347, 6322675, 3506438, 213720], [8404992, 3585646, 2654220, 7690526, 5121240, 4044637], [1, 1, 8404992, 3585646, 5750773, 7720725], [8404992, 4819347, 6322675, 3506438, 8191273, 6807119], [1, 8404992, 3585646, 5750773, 684268, 3164248], [8404992, 3585646, 5750773, 684268, 5240745, 3700771], [1, 1, 1, 1, 1, 8404992], [8404992, 4819347, 2082318, 7134899, 2971397, 5785355], [1, 8404992, 4819347, 2082318, 7134899, 5433596], [8404992, 3585646, 2654220, 714467, 6025584, 5597909], [1, 1, 8404992, 4819347, 2082318, 1270094], [8404992, 4819347, 6322675, 4898555, 7431338, 2786978], [1, 8404992, 3585646, 2654220, 714467, 2379409], [8404992, 3585646, 5750773, 7720725, 6753480, 5400177], [1, 1, 1, 8404992, 4819347, 6322675], [8404992, 4819347, 2082318, 1270094, 1484163, 5146046], [1, 8404992, 4819347, 6322675, 4898555, 973655], [8404992, 3585646, 2654220, 7690526, 3283753, 2436145], [1, 1, 8404992, 3585646, 2654220, 7690526], [8404992, 4819347, 6322675, 3506438, 213720, 4153273], [1, 8404992, 3585646, 5750773, 7720725, 1651513], [8404992, 3585646, 5750773, 684268, 3164248, 5330453], [1, 1, 1, 1, 8404992, 3585646], [8404992, 4819347, 2082318, 7134899, 5433596, 2134075], [1, 8404992, 4819347, 2082318, 1270094, 6920830], [8404992, 3585646, 2654220, 714467, 2379409, 3463953], [1, 1, 8404992, 4819347, 6322675, 3506438], [8404992, 4819347, 6322675, 4898555, 973655, 8314555], [1, 8404992, 3585646, 2654220, 7690526, 5121240], [8404992, 3585646, 5750773, 7720725, 1651513, 5639303], [1, 1, 1, 8404992, 3585646, 5750773], [8404992, 4819347, 2082318, 1270094, 6920830, 2126834], [1, 8404992, 4819347, 6322675, 3506438, 8191273], [8404992, 3585646, 2654220, 7690526, 5121240, 4360356], [1, 1, 8404992, 3585646, 5750773, 684268], [8404992, 4819347, 6322675, 3506438, 8191273, 1597874], [1, 8404992, 3585646, 5750773, 684268, 5240745], [8404992, 3585646, 5750773, 684268, 5240745, 4704222]]
nthroots = [GaussianInteger(1,0), GaussianInteger(458475,2159086), GaussianInteger(7767882,3021529), GaussianInteger(3735285,1895821), GaussianInteger(73026,1633051), GaussianInteger(520345,4836327), GaussianInteger(4253900,4818341), GaussianInteger(4071267,1224421), GaussianInteger(2311455,1878291), GaussianInteger(1299208,4939937), GaussianInteger(3328455,5617398), GaussianInteger(4422782,4095225), GaussianInteger(5607962,2515739), GaussianInteger(2460988,2080420), GaussianInteger(5875920,4596116), GaussianInteger(6581058,7803074), GaussianInteger(7427812,5598482), GaussianInteger(1238430,368946), GaussianInteger(3770340,6890108), GaussianInteger(8145273,1730137), GaussianInteger(7327962,1448761), GaussianInteger(4153459,6177308), GaussianInteger(2937626,4766109), GaussianInteger(7435388,6285825), GaussianInteger(6728461,2667480), GaussianInteger(6421274,5952093), GaussianInteger(4254147,7957992), GaussianInteger(1270078,585163), GaussianInteger(5592766,6652379), GaussianInteger(5794663,7080751), GaussianInteger(8180550,1891210), GaussianInteger(4275095,3154194), GaussianInteger(8119042,8119042), GaussianInteger(3736060,4061111), GaussianInteger(5477944,3130628), GaussianInteger(4701662,8167134), GaussianInteger(4110293,4989465), GaussianInteger(3816734,5424227), GaussianInteger(1387812,1513245), GaussianInteger(7136469,3499536), GaussianInteger(702877,2453760), GaussianInteger(2450320,4492243), GaussianInteger(3519904,818333), GaussianInteger(8395278,731771), GaussianInteger(6219506,8350675), GaussianInteger(4049196,1689450), GaussianInteger(566609,8119832), GaussianInteger(7518234,2899871), GaussianInteger(2806511,977181), GaussianInteger(6683642,4672622), GaussianInteger(3227341,7525150), GaussianInteger(382893,968251), GaussianInteger(3294709,1381241), GaussianInteger(4257377,1330300), GaussianInteger(7142782,7508357), GaussianInteger(586817,7443711), GaussianInteger(659942,3362454), GaussianInteger(6127128,5996249), GaussianInteger(5013665,7917472), GaussianInteger(1017721,6733976), GaussianInteger(5804306,220912), GaussianInteger(6085966,2434992), GaussianInteger(6909842,8247634), GaussianInteger(1646204,8196862)]
invNthroots = [GaussianInteger(1,0), GaussianInteger(8196862,6758789), GaussianInteger(8247634,1495151), GaussianInteger(2434992,2319027), GaussianInteger(220912,2600687), GaussianInteger(6733976,7387272), GaussianInteger(7917472,3391328), GaussianInteger(5996249,2277865), GaussianInteger(3362454,7745051), GaussianInteger(7443711,7818176), GaussianInteger(7508357,1262211), GaussianInteger(1330300,4147616), GaussianInteger(1381241,5110284), GaussianInteger(968251,8022100), GaussianInteger(7525150,5177652), GaussianInteger(4672622,1721351), GaussianInteger(977181,5598482), GaussianInteger(2899871,886759), GaussianInteger(8119832,7838384), GaussianInteger(1689450,4355797), GaussianInteger(8350675,2185487), GaussianInteger(731771,9715), GaussianInteger(818333,4885089), GaussianInteger(4492243,5954673), GaussianInteger(2453760,7702116), GaussianInteger(3499536,1268524), GaussianInteger(1513245,7017181), GaussianInteger(5424227,4588259), GaussianInteger(4989465,4294700), GaussianInteger(8167134,3703331), GaussianInteger(3130628,2927049), GaussianInteger(4061111,4668933), GaussianInteger(8119042,285951), GaussianInteger(3154194,4129898), GaussianInteger(1891210,224443), GaussianInteger(7080751,2610330), GaussianInteger(6652379,2812227), GaussianInteger(585163,7134915), GaussianInteger(7957992,4150846), GaussianInteger(5952093,1983719), GaussianInteger(2667480,1676532), GaussianInteger(6285825,969605), GaussianInteger(4766109,5467367), GaussianInteger(6177308,4251534), GaussianInteger(1448761,1077031), GaussianInteger(1730137,259720), GaussianInteger(6890108,4634653), GaussianInteger(368946,7166563), GaussianInteger(5598482,977181), GaussianInteger(7803074,1823935), GaussianInteger(4596116,2529073), GaussianInteger(2080420,5944005), GaussianInteger(2515739,2797031), GaussianInteger(4095225,3982211), GaussianInteger(5617398,5076538), GaussianInteger(4939937,7105785), GaussianInteger(1878291,6093538), GaussianInteger(1224421,4333726), GaussianInteger(4818341,4151093), GaussianInteger(4836327,7884648), GaussianInteger(1633051,8331967), GaussianInteger(1895821,4669708), GaussianInteger(3021529,637111), GaussianInteger(2159086,7946518)]

# NTT_6
theta_powers = [1, 5763158, 4509843, 1405441, 4458480, 3769554, 6999551, 6588348, 125596]
zeta = 1405441
zeta_squared = 6999551
one_minus_zeta_by_nine = 8248833
zeta_minus_one_by_nine = 156160
theta_inv = 125596
theta_squared_inv = 6588348
