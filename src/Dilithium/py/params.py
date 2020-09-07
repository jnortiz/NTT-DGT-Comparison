from gaussian import GaussianInteger

modinv = lambda y,p:pow(y,p-2,p)

p = 8380417
N = 256
k = N//2

# k-th root of i mod p
if p == 12289:
	kthroots = {
	2: GaussianInteger(1550, 1550),
	4: GaussianInteger(4730, 9782),
	8: GaussianInteger(10885, 9419),
	16: GaussianInteger(3512,1024),
	32: GaussianInteger(959, 2358),
	64: GaussianInteger(11255, 11828),
    512:GaussianInteger(8458, 7974)
	}

	invkthroots = {
	2: GaussianInteger(1550, 10739),
	4: GaussianInteger(4730, 2507),
	8: GaussianInteger(11341, 2358),
	16: GaussianInteger(8777,1024),
	32: GaussianInteger(8622, 5044),
	64: GaussianInteger(6207, 9177),
    512: GaussianInteger(1330, 3875)
	}

elif p == 0xFFFFFFFF00000001:
	kthroots = {
	    4: GaussianInteger(17872535116924321793, 18446744035054843905),
	    8: GaussianInteger(18446741870391328801, 18293621682117541889),
	    16: GaussianInteger(13835058050988244993, 68719460336),
	    32: GaussianInteger(2711615452560896026, 1702263320480097066),
	    64: GaussianInteger(5006145799600815370, 13182758589097755684),
	    128: GaussianInteger(1139268601170473317, 15214299942841048380),
	    256: GaussianInteger(4169533218321950981, 11340503865284752770),
	    512: GaussianInteger(1237460157098848423, 590072184190675415),
	    1024: GaussianInteger(13631489165933064639, 9250462654091849156),
	    2048: GaussianInteger(12452373509881738698, 10493048841432036102),
	    4096: GaussianInteger(12694354791372265231, 372075061476485181),
	    8192: GaussianInteger(9535633482534078963, 8920239275564743446),
	    16384: GaussianInteger(9868966171728904500, 6566969707398269596),
	    32768: GaussianInteger(10574165493955537594, 3637150097037633813),
	    65536: GaussianInteger(2132094355196445245, 12930307277289315455)
	}

	invkthroots = {
	    4: GaussianInteger(34359736320, 17868031517297999873),
	    8: GaussianInteger(18311636080627023873, 18446741870391328737),
	    16: GaussianInteger(18446739675663041537, 18446462594505048065),
	    32: GaussianInteger(9223372049739972605, 9223372049739382781),
	    64: GaussianInteger(3985917792403544159, 10871216858344511621),
	    128: GaussianInteger(697250266387245748, 7269985899340929883),
	    256: GaussianInteger(16440350318199496391, 8259263625103887671),
	    512: GaussianInteger(11254465366323603399, 282547220712277683),
	    1024: GaussianInteger(4772545667722300316, 8077569763565898552),
	    2048: GaussianInteger(13028894352332048345, 9995848711590186809),
	    4096: GaussianInteger(11525613835860693, 17335883825168514904),
	    8192: GaussianInteger(17414606149056687587, 3916527805974289959),
	    16384: GaussianInteger(9801605401783064476, 2749242888134484347),
	    32768: GaussianInteger(10469048769509713349, 8715957816394874924),
	    65536: GaussianInteger(15132804493885713016, 7997468840100395569)
	}
elif p == 8380417:
	kthroots = {
		128:GaussianInteger(2089168, 5159577)
	}

	invkthroots = {
		128:GaussianInteger(4799789, 4692735)
	}
else:
	raise Exception()

# Primitive root of p
PROOTS = {
    0xFFFFFFFF00000001:7,
    12289:11,
	8380417:10
}

# Computing the k-th primitive root of p
r = PROOTS[p] # Primitive root of p
assert (p-1)%k == 0
n = (p-1)//k

g = int(pow(r, n, p))
g_inv = modinv(g, p)

assert pow(g, k, p) == 1 # k-th primitive root of unity
assert pow(g_inv, k, p) == 1 # Inverse of the k-th primitive root of unity

# Scalar k^-1 mod p used in to scale the output signal in backward DGT
invkmodp = modinv(k, p)

w = [GaussianInteger(4616857, 3763560), GaussianInteger(3042578, 27908), GaussianInteger(14629305284132, 134187078152), GaussianInteger(4644566, 1062415), GaussianInteger(22331974373804, 5108297428510), GaussianInteger(21433789622460, 4902843366150), GaussianInteger(21452804475664, 4907192893160), GaussianInteger(3120684, 3716310), GaussianInteger(15004854084696, 17868739444140), GaussianInteger(14401363730040, 17150064551100), GaussianInteger(14414139810336, 17165279124240), GaussianInteger(8999197588584, 10716819771060), GaussianInteger(16161782143332, 19246483334130), GaussianInteger(16161981867108, 19246721177970), GaussianInteger(9816667003752, 11690314608180), GaussianInteger(2227391, 3639383), GaussianInteger(10709728041854, 17498859504302), GaussianInteger(10278986260710, 16795061062230), GaussianInteger(10288105199464, 16809960696232), GaussianInteger(6423185338866, 10494983381058), GaussianInteger(11535486479893, 18848084324509), GaussianInteger(11535629032917, 18848317245021), GaussianInteger(7006654866098, 11448327036674), GaussianInteger(13922363130275, 22748054426075), GaussianInteger(5974504150608, 9761873348304), GaussianInteger(1340183299053, 2189754881589), GaussianInteger(10775966195412, 17607087475956), GaussianInteger(17424788469969, 28470743994297), GaussianInteger(2720038659207, 4444330813791), GaussianInteger(16427405100598, 26841097435174), GaussianInteger(10280634530050, 16797754205650), GaussianInteger(567084, 1577601), GaussianInteger(2726649886296, 7585411662594), GaussianInteger(2616984914040, 7280328870810), GaussianInteger(2619306555936, 7286787569304), GaussianInteger(1635314874984, 4549369021326), GaussianInteger(2936884370532, 8170274103723), GaussianInteger(2936920663908, 8170375070187), GaussianInteger(1783863662952, 4962624758478), GaussianInteger(3544572719100, 9860834490525), GaussianInteger(1521082608192, 4231580231088), GaussianInteger(341204802372, 949215702483), GaussianInteger(2743513830288, 7632326361132), GaussianInteger(4436274881556, 12341507941359), GaussianInteger(692509937868, 1926530056377), GaussianInteger(4182345440952, 11635088187978), GaussianInteger(2617404556200, 7281496295550), GaussianInteger(3994812463404, 11113380270081), GaussianInteger(4368717030468, 12153565179027), GaussianInteger(2753745726900, 7660791015975), GaussianInteger(4506215052528, 12536078205492), GaussianInteger(2448308046576, 6811077763764), GaussianInteger(1911270992316, 5317065952749), GaussianInteger(949647372660, 2641874298615), GaussianInteger(1510359049752, 4201747796178), GaussianInteger(2719348679796, 7565100049719), GaussianInteger(364424056752, 1013810575428), GaussianInteger(3457331949456, 9618134775084), GaussianInteger(1596703826676, 4441954902039), GaussianInteger(1391614495572, 3871406034783), GaussianInteger(828349806312, 2304430177518), GaussianInteger(3780251128248, 10516480733322), GaussianInteger(2921923560444, 8128653834141), GaussianInteger(2719536, 8139413), GaussianInteger(13076056677984, 39135876750122), GaussianInteger(12550141928160, 37561844506530), GaussianInteger(12561275708544, 37595167263352), GaussianInteger(7842396671136, 23471836892838), GaussianInteger(14084267539728, 42153393192199), GaussianInteger(14084441590032, 42153914114631), GaussianInteger(8554784565408, 25603972407014), GaussianInteger(16998527756400, 50875604441825), GaussianInteger(7294578778368, 21832249816944), GaussianInteger(1636298579088, 4897346432079), GaussianInteger(13156930239552, 39377926613916), GaussianInteger(21274818627024, 63674294183067), GaussianInteger(3321034813872, 9939663949101), GaussianInteger(20057062077408, 60029619690514), GaussianInteger(12552154384800, 37567867672150), GaussianInteger(19157719680816, 57337940229653), GaussianInteger(20950834864272, 62704629633551), GaussianInteger(13205998827600, 39524786042675), GaussianInteger(21610227160512, 64678152406596), GaussianInteger(11741226823104, 35140808667332), GaussianInteger(9165785438064, 27432662465137), GaussianInteger(4554175778640, 13630383100995), GaussianInteger(7243152352608, 21678333517114), GaussianInteger(13041042651984, 39031081807747), GaussianInteger(1747649980608, 5230614697364), GaussianInteger(16580151618624, 49623429006492), GaussianInteger(7657231623504, 22917648679907), GaussianInteger(6673695111888, 19973981131979), GaussianInteger(3972475186848, 11889387078534), GaussianInteger(18128758759392, 54258320066386), GaussianInteger(14012520740976, 41938659198433), GaussianInteger(8816659564992, 26387749042436), GaussianInteger(5863455592800, 17548981398650), GaussianInteger(10929785269104, 32712211313457), GaussianInteger(5634413551344, 16863471896377), GaussianInteger(1071899675328, 3208133355124), GaussianInteger(2525766340464, 7559471684337), GaussianInteger(19811681063664, 59295208594937), GaussianInteger(9535726639680, 28539874954940), GaussianInteger(7453005348048, 22306411321259), GaussianInteger(13037376717456, 39020109879023), GaussianInteger(14134894421904, 42304916504607), GaussianInteger(9796361530848, 29319940018034), GaussianInteger(19636523908512, 58770973421846), GaussianInteger(10075399522128, 30155084488899), GaussianInteger(8459981540448, 25320232469834), GaussianInteger(17075985580752, 51107431202891), GaussianInteger(21276455787696, 63679194109693), GaussianInteger(12265379313600, 36709566571300), GaussianInteger(8681721627744, 25983887648202), GaussianInteger(7879349726304, 23582435236682), GaussianInteger(5054029214976, 15126415349808), GaussianInteger(13691378893344, 40977500335502), GaussianInteger(7933805715168, 23745418842594), GaussianInteger(5041492154016, 15088892655878), GaussianInteger(12574120077072, 37633609710951), GaussianInteger(5264475069264, 15756267545987), GaussianInteger(18121435048944, 54236400627177), GaussianInteger(12048534391104, 36060562336332), GaussianInteger(13069956758736, 39117620046763), GaussianInteger(2223318583296, 6654263146368), GaussianInteger(9287403087984, 27796657014497), GaussianInteger(4784607502992, 14320051843311)]


w_inv = [GaussianInteger(3892641, 8377596), GaussianInteger(7195353, 4400568), GaussianInteger(5425935, 6760166), GaussianInteger(2345055, 3591409), GaussianInteger(861450, 7476452), GaussianInteger(2380950, 2855496), GaussianInteger(97464, 2674952), GaussianInteger(7098624, 5986556), GaussianInteger(4654163, 5678949), GaussianInteger(8338474, 5182144), GaussianInteger(3515619, 588930), GaussianInteger(590683, 2929378), GaussianInteger(6537316, 4487296), GaussianInteger(6267523, 4503509), GaussianInteger(2438943, 7598314), GaussianInteger(7539166, 333657), GaussianInteger(2122408, 3675083), GaussianInteger(1695586, 7854831), GaussianInteger(3495343, 5900731), GaussianInteger(2141266, 461601), GaussianInteger(6278349, 6662720), GaussianInteger(7485261, 2129880), GaussianInteger(852284, 2683516), GaussianInteger(1674151, 1601861), GaussianInteger(7077742, 4535706), GaussianInteger(7213984, 5914144), GaussianInteger(872625, 4361651), GaussianInteger(6046804, 2537943), GaussianInteger(6827326, 7254674), GaussianInteger(7890630, 5104097), GaussianInteger(8323600, 2074686), GaussianInteger(2325132, 2512194), GaussianInteger(4028641, 1911023), GaussianInteger(8369446, 4810516), GaussianInteger(5399004, 1654381), GaussianInteger(5994602, 6453482), GaussianInteger(560883, 2855097), GaussianInteger(676132, 257878), GaussianInteger(7096646, 5237512), GaussianInteger(7117390, 3322915), GaussianInteger(8228714, 206518), GaussianInteger(2719536, 241004), GaussianInteger(6547725, 5768336), GaussianInteger(5626501, 4830911), GaussianInteger(2704403, 2416219), GaussianInteger(3771111, 5062776), GaussianInteger(2442536, 1519094), GaussianInteger(3737561, 4627620), GaussianInteger(1874697, 5257626), GaussianInteger(2614646, 1320979), GaussianInteger(3632328, 925684), GaussianInteger(2179874, 4735672), GaussianInteger(5768330, 4791574), GaussianInteger(158407, 4908186), GaussianInteger(2813177, 6179570), GaussianInteger(2341091, 4327329), GaussianInteger(6053316, 2605920), GaussianInteger(5786376, 2058645), GaussianInteger(316440, 5310825), GaussianInteger(1699075, 5676464), GaussianInteger(7789561, 3478960), GaussianInteger(5671898, 7244185), GaussianInteger(7584940, 6414725), GaussianInteger(2180572, 2568065), 

GaussianInteger(5638723, 2792435), GaussianInteger(7950045, 6938739), GaussianInteger(6776344, 538031), GaussianInteger(3659471, 1878133), GaussianInteger(6614241, 340960), GaussianInteger(5648368, 488551), GaussianInteger(6168709, 3215634), GaussianInteger(7744036, 7956652), GaussianInteger(4452117, 5439232), GaussianInteger(1623507, 1990730), GaussianInteger(8116917, 1324641), GaussianInteger(7676940, 2772480), GaussianInteger(7369224, 7058762), GaussianInteger(4247471, 3994974), GaussianInteger(4248781, 128308), GaussianInteger(4350637, 3022120), GaussianInteger(6176728, 989957), GaussianInteger(754550, 3940002), GaussianInteger(8100332, 6776412), GaussianInteger(3526258, 3538542), GaussianInteger(7813333, 1577601), GaussianInteger(5791593, 5459118), GaussianInteger(7222586, 1773030), GaussianInteger(6956199, 7550434), GaussianInteger(74803, 66641), GaussianInteger(3401024, 2387641), GaussianInteger(8336581, 4334946), GaussianInteger(4504634, 609439), GaussianInteger(7577509, 5003657), GaussianInteger(6152515, 6221061), GaussianInteger(2559490, 2105668), GaussianInteger(577602, 5338278), 

GaussianInteger(3773247, 6178426), GaussianInteger(1384455, 7347713), GaussianInteger(1509426, 891235), GaussianInteger(4074113, 3268773), GaussianInteger(7552727, 4777718), GaussianInteger(1666900, 2929100), GaussianInteger(122615, 4468765), GaussianInteger(5628640, 4552979), GaussianInteger(2550194, 6941215), GaussianInteger(4785165, 691778), GaussianInteger(2227391, 4741034), GaussianInteger(722045, 4081197), GaussianInteger(6348799, 5577633), GaussianInteger(6673101, 6905655), GaussianInteger(2277177, 6981314), GaussianInteger(6587166, 5174491), 

GaussianInteger(5491743, 7489321), GaussianInteger(4118972, 7206038), GaussianInteger(4085511, 4297157), GaussianInteger(2622542, 1631528), GaussianInteger(2380374, 2444094), GaussianInteger(3120684, 4664107), GaussianInteger(2614722, 6296882), GaussianInteger(1440907, 6199820), 

GaussianInteger(6798825, 6106555), GaussianInteger(1268623, 2183275), GaussianInteger(3735851, 1062415), GaussianInteger(6760544, 2514674), 

GaussianInteger(5206831, 158852), GaussianInteger(3042578, 8352509), 

GaussianInteger(3763560, 3763560)]