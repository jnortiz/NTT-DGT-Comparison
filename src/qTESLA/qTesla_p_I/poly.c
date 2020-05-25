/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: DGT, modular reduction and polynomial functions
**************************************************************************************/

#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"


/* The k-th root of i and its inverse are multiplied by the Montgomery constant R = 2^32 mod p. */
static const int32_t nthroots[PARAM_N] = {172048372, 0, 317977323, 337163154, 265441122, 76896845, 204070189, 98753319, 173436411, 142249661, 48289095, 31690006, 168154760, 115504205, 85341862, 207021055, 300948020, 234654008, 274344386, 333694015, 55562192, 315132040, 155549124, 70569583, 280982638, 39696014, 39958875, 330723128, 309004679, 325769056, 103548902, 262123333, 135490431, 104460452, 315925655, 304595989, 138025277, 111496939, 236561099, 257296478, 238187521, 328527502, 8606419, 270750620, 233439059, 246526343, 129835268, 89487289, 71753438, 75103137, 184123196, 4841271, 290795988, 206477384, 93609091, 142549013, 30420351, 259598127, 272706097, 89889337, 197602390, 198147530, 41349798, 208102942, 276947471, 212904634, 64128521, 204053129, 300577464, 103019415, 156008152, 277638475, 94764841, 221216756, 61791730, 233480537, 240400233, 235691747, 167451083, 119190363, 124028559, 95350249, 141545700, 323051450, 66477766, 271364016, 144823751, 124192166, 330816348, 297731631, 210813716, 177596718, 28110246, 160407028, 161482893, 3184410, 309222922, 255673257, 268383012, 35142603, 309459503, 42488211, 43925281, 94549066, 213717722, 1517282, 335753719, 228154944, 288972923, 266836964, 8075056, 64623876, 153028338, 84977571, 332420586, 33773389, 135695830, 46625988, 106281800, 101909431, 7253775, 93503043, 284368345, 219071892, 299896934, 253140660, 324504721, 248987403, 243069542, 308135050, 178868, 292966729, 169309271, 176743269, 66718820, 192394657, 82115121, 271831030, 227109347, 104238958, 67854916, 228166895, 82982195, 68083191, 56595181, 261737414, 190600807, 193627187, 147853830, 256041539, 155297494, 11609903, 343415422, 179336289, 101970814, 86460626, 258654802, 168918972, 26813551, 216066481, 17277594, 321162525, 42695862, 248468402, 85054974, 30269591, 57808436, 220971825, 221688099, 3697712, 241810146, 247647950, 14935846, 34495046, 202625515, 238946676, 84938389, 145000049, 116335219, 28328801, 4325286, 84887642, 217951388, 226634933, 167174960, 11538274, 321522020, 283851712, 68116813, 136022754, 301711822, 106971414, 229061660, 305233635, 83623516, 200402907, 263065853, 67637195, 136350879, 181203944, 107267729, 230218976, 266536939, 24012695, 96179371, 3927930, 199652456, 228908062, 42417057, 263833527, 29796697, 342399674, 266693627, 227840873, 134527122, 302697183, 45908500, 3416778, 89482263, 166928435, 282786016, 199798000, 133276089, 116011422, 156964060, 61742845, 308137782, 132295548, 86892150, 269964373, 105852138, 139489212, 46203987, 171944813, 242742603, 210216010, 157876616, 161070465, 180791527, 124463535, 139795889, 215432472, 197951505, 290043569, 274630233, 310912488, 102600989, 160855644, 65036177, 48363408, 156545577, 243457722, 224997181, 89100652, 165442463, 145712114, 312193573, 18174287, 226131565, 306270915, 160542112, 238802755, 159331653, 142629279, 312543555, 38874948, 5228512, 56883520, 271133466, 55716246, 38822243, 273926000, 208539278, 278974226, 169041604, 269922432, 108998200, 160650051, 240207849, 80783147, 203353715, 27800405, 333749405, 161413484, 133515974, 150675107, 30399932, 248870171, 286511333, 343484057, 125776071, 335916371, 211401750, 146912568, 161616182, 226481386, 168779456, 54327512, 293680271, 339959723, 62847388, 52336156, 234112811, 31414543, 269002188, 204196911, 125967393, 83828181, 63939080, 230422429, 48636152, 68531512, 26637208, 175191081, 196957072, 30111891, 40728699, 32792426, 22284080, 67029858, 141461549, 267814432, 115557732, 299148003, 99447361, 18203990, 227760953, 222852895, 338969550, 288078750, 177448079, 59120435, 285662014, 294555538, 204480372, 217163185, 166022566, 24340454, 31069144, 267006181, 74423906, 137701606, 134418290, 69807788, 62243660, 341716174, 164533240, 275813578, 183775820, 118247387, 101973929, 271665061, 79445846, 331742648, 256353272, 219173729, 137266479, 237916938, 263497105, 106914242, 94153433, 46380871, 23419237, 27380569, 104936206, 266643298, 11415414, 89232306, 91176027, 144871749, 225584207, 182345818, 94406327, 312122428, 222733932, 341692373, 98194032, 32437955, 338476162, 30531729, 45780758, 157011125, 146970484, 170364917, 282748758, 330442944, 252001412, 334213998, 160024716, 188248941, 280347330, 52521614, 13316568, 232092012, 139871452, 67548107, 64326038, 146041310, 72040520, 94987773, 133532664, 338505894, 202156358, 126509152, 45244050, 68985178, 61124071, 11070163, 4030559, 329255563, 76607764, 57019610, 43083837, 336841830, 173611634, 222785699, 119225879, 9448901, 223102514, 20476390, 313335986, 106407, 231728308, 163889921, 46689205, 24942895, 40441970, 129794145, 52408164, 222309586, 138537916, 295855913, 271605595, 219790977, 341226172, 45152110, 314760234, 130376028, 17343362, 25020606, 199989291, 291040221, 341156339, 162237396, 246799407, 76712571, 16117992, 342596201, 62365488, 319982580, 321016439, 305336736, 34481930, 28724132, 166638686, 289196767, 315569638, 21684185, 297600534, 325735133, 43911902, 74599599, 112186011, 37805887, 100813847, 128930178, 88279206, 77156786, 306650028, 103280324, 182010936, 133182162, 115690681, 92632688, 282909261, 43481354, 321432787, 41069750, 55370485, 22474417, 95520362, 231113820, 290332258, 102582977, 245484388, 76864394, 156740554, 271202042, 186287315, 89102724, 315235438, 243049873, 241690149, 250470760, 12218043, 151210962, 267972452, 196292367, 218174059, 46161349, 103516582, 26582427, 65618741, 155056578, 166419013, 153414236, 99672984, 219435283, 268093071, 326902272, 43034682, 191040688, 258652380, 84924197, 276469363, 282904658, 52819556, 227690368, 234692762, 158061995, 228066720, 67302255, 116815348, 254415434, 260103374, 90962482, 261528800, 93677043, 237653145, 303917590, 54906072, 69127916, 271863862, 302164160, 138700242, 40508470, 204033325, 3168991, 75738315, 287982743, 141531868, 148110025, 186150144, 136787397, 177855349, 84338616, 1003738, 104309841, 256162323, 298912938, 311849720, 81188914, 307092122, 209481473, 259089858, 168736599, 140491628, 61562455, 291439513, 15060091, 1038892, 311507063, 109653926, 277588017, 232662043, 304210161, 214853136, 52932452, 126168333, 147711970, 92823313, 91073409, 312335334, 139577399, 326503374, 242479214, 298455455, 86335693, 223140343, 167204823, 250079316, 161333767, 18488084, 244444252, 126924547, 174195866, 94050490, 218355419, 118409369, 281576985, 235032396, 85276325, 179818327, 64308219, 205204712, 274180384, 287051822, 84580043, 273819399, 199405222, 249764876, 162712060, 196815048, 83128691, 269700337, 314892636, 79973082, 324726945, 75590890, 295510164, 268386619, 300343844, 332346174, 214154901, 262903158, 11969403, 170701507, 6676159, 242680830, 137487347, 137490235, 212318257, 249185512, 241438758, 337231183, 5256856, 142924806, 74372519, 328325639, 65510128, 116606959, 134115208, 203117978, 163312741, 307046743, 178827105, 18732571, 251585580, 84372223, 234179784, 143510264, 33615069, 275496748, 159203797, 172572797, 253742099, 44262436, 83152948, 260221580, 147959514, 189107245, 230670902, 7965992, 334838309, 86807173, 70281503, 259506592, 181117068, 178292440, 191320938, 75707651, 272061324, 271530146, 283512113, 313422018, 107535666, 182856453, 183610938, 59331568, 275233086, 286636399, 9460332, 214417832, 3818356, 166675762, 73830309, 9046055, 217458059, 331601142, 169012508, 168545439, 96508064, 108949165, 263225470, 168259908, 59213403, 76886693, 186922429, 307672689, 182639393, 864480, 148606694, 177696284, 47374831, 220228586, 120492156, 236224691, 29220771, 111140949, 313499326, 29434284, 1352417, 52170072, 238095442, 92513845, 125527693, 277884843, 119843787, 55534047, 24621919, 322819025, 175844749, 194776089, 144079558, 81678223, 287047284, 230958701, 131922387, 299266096, 102767987, 167855468, 313912845, 221448855, 214984484, 264856617, 181294336, 312362586, 254303480, 243746981, 44233948, 45191201, 60776271, 189254037, 285449463, 309435470, 91315870, 86604762, 302262587, 200231566, 203364036, 212580482, 9036368, 242682650, 27709740, 247537429, 286787125, 306778421, 96793133, 329636128, 325143363, 9374916, 255210719, 13556380, 162121436, 96749849, 254809161, 296714509, 20797919, 79145552, 149943272, 99764890, 205638593, 17890177, 255111514, 235751791, 240586572, 123465304, 228682775, 291757248, 55681685, 147240129, 220957158, 291985757, 157608758, 45510709, 276900146, 284185915, 242839135, 97262770, 41447375, 117450698, 198518246, 311067983, 116963577, 207556672, 10337180, 342126799, 223241259, 68033082, 94852863, 214622357, 250211057, 63521643, 57716712, 59746153, 36121568, 314733284, 221482256, 158674286, 116350752, 332221775, 146259259, 150500774, 272672280, 107404815, 219412596, 89382433, 310394093, 306649757, 322551788, 335414745, 149032672, 200944142, 73483759, 54817677, 338162693, 3589109, 329739756, 225565627, 265148372, 241197186, 8265065, 232472823, 181031655, 46812431, 237773843, 333327594, 25157527, 155532001, 131466927, 170842362, 257540671, 261103408, 199174350, 209290436, 128608200, 75406016, 6589908, 295231045, 204464303, 124704213, 184008413, 143671236, 337264105, 317110889, 305818682, 82183594, 266370590, 258222890, 237005321, 318579388, 258996687, 159267551, 286600249, 264914731, 233665129, 82081910, 20706757, 86331240, 24064167, 40502783, 134090755, 209970994, 68527024, 49895356, 232306481, 99281949, 271661466, 59297367, 151578868, 74728939, 33857048, 243815934, 225209918, 55319038, 237379042, 342109586, 242752727, 233827174, 315149473, 333300167, 175715761, 10161331, 105265715, 17690627, 85996922, 248754700, 107090821, 177273476, 307995444, 271708720, 116149346, 287229325, 39679567, 329231345, 10674953, 126485734, 329143637, 342241374, 205778498, 129794946, 42660821, 302303357, 103577843, 97914650, 212990101, 323011035, 257180954, 65861960, 65323301, 270762367, 308062676, 64458518, 332695523, 306450490, 222727581, 245643478, 139397526, 90748078, 269509877, 242764444, 174842333, 153388031, 229578666, 143033511, 201412194, 124519330, 299525061, 283267670, 244850926, 236841026, 228347502, 246143805, 210422994, 307985168, 316803311, 157789452, 8917219, 55895108, 246419072, 33521851, 235859389, 191641212, 296212518, 94649066, 156596475, 130578970, 287604226, 46672810, 119756093, 236227847, 73071185, 132276796, 282854881, 129682373, 262693114, 146599465, 17902158, 29318594, 233176987, 239606852, 44571981, 196327301, 53274331, 33504920, 281238405, 204092013, 191927788, 214380818, 66458035, 176624631, 251446012, 85681972, 264054689, 109857654, 74716188, 25968793, 54526135, 73781214, 329557101, 4338060, 157167826, 163508655, 113705204, 59742204, 82790542, 79138592, 336007060, 128508129, 304559401, 43414910, 336016320, 38978615, 298264654, 336820096, 167610646, 40707749, 109571358, 189314684, 262027071, 266799759, 64779062, 218622044, 24121209, 97215210, 336310543, 336794933, 236158331, 99364986, 277615850, 172886890, 231050835, 133851883, 116089994, 56898903, 286872571, 34291768, 132712867, 237051548, 60124269, 307918898, 129750801, 237145757, 95789606};


static const int32_t invnthroots[1024] = {8388608, 0, 318263842, 104891361, 160633892, 32950996, 270549775, 104220497, 62666669, 37511712, 3915538, 230058334, 251198691, 157434851, 288330864, 242581707, 187093563, 81673784, 89710631, 206696030, 162379425, 156835359, 90638592, 61309420, 135678218, 95438771, 130024036, 82840197, 151870855, 110643413, 157352597, 257024581, 291817381, 36831510, 43603354, 256255542, 144870163, 150734808, 73129448, 193107257, 248449461, 340775701, 180062945, 288902336, 202963465, 265726536, 271076011, 320616841, 220881289, 102619624, 324933163, 250757387, 64265292, 206515437, 174963586, 337192175, 51800497, 183448440, 13795819, 72745615, 83088198, 101933857, 196721116, 260654221, 120204647, 45163282, 68231272, 301914338, 291269798, 27226641, 176327649, 303731459, 21354939, 40004527, 289364406, 176024242, 131088259, 285775291, 342796208, 17192211, 175449381, 156840371, 318655293, 82835653, 45420879, 125420504, 258163720, 45522094, 191266097, 179532682, 183814861, 181923523, 203067479, 212241467, 270207288, 43826588, 184345377, 57156871, 40176842, 280254863, 278207256, 111785654, 58158524, 342605943, 23828169, 104209339, 179696203, 116585110, 71403349, 143124500, 66868767, 38322248, 20110188, 186425449, 30127797, 256482552, 147087097, 219975108, 154843347, 318116920, 244006421, 189044297, 180714213, 83290564, 28938386, 258771023, 299689320, 62410072, 115392067, 154093998, 165769706, 332867827, 290641284, 265845062, 278040727, 31679623, 137495346, 88232099, 82748066, 7566727, 90088395, 1978592, 309558726, 120097746, 53355987, 304675864, 236153375, 239778491, 199775383, 269093114, 20595071, 342126436, 302411464, 341758278, 327537552, 329338614, 180136917, 58265361, 68977485, 19266482, 139360660, 297847861, 53817682, 204930589, 341825329, 342826422, 231558561, 241408665, 39632275, 79023348, 273572914, 200125941, 305886607, 149332635, 33387205, 125534799, 338000046, 199467963, 323965392, 285705934, 305253095, 69840683, 16092823, 303033089, 195634360, 314477950, 226542521, 45725689, 29538983, 128693939, 37829876, 131116638, 205729701, 203488625, 301803562, 167428324, 226399948, 150682027, 70509176, 156374202, 331291067, 181762578, 188918066, 165294808, 101344391, 86765152, 30043980, 296833707, 14064983, 336188039, 340881811, 266298991, 11551339, 9002203, 193552904, 316750598, 14050944, 104755611, 185144439, 21970009, 266163549, 9855945, 328003988, 282888310, 132482119, 224152210, 150542000, 245293658, 161484103, 152384231, 21544086, 242131637, 16217880, 243466359, 182342702, 99566968, 172644596, 38787907, 315157531, 138238720, 67124990, 42541688, 176714069, 224193526, 17834979, 307110589, 32962304, 119256578, 280301312, 210825070, 138776709, 35476656, 135859526, 233625468, 17878805, 129224686, 243699178, 128271378, 328589119, 58811081, 121258536, 245814164, 111892233, 258318539, 42677662, 103146539, 16398010, 53529259, 22185205, 179932392, 295087747, 188375524, 196262660, 104657011, 194431331, 131507098, 127463118, 192618004, 313568466, 206754656, 308571164, 185396825, 142316297, 188761547, 289910386, 86821044, 212448364, 94897739, 132786813, 317907603, 237729344, 51737378, 178385239, 271404804, 118894200, 22056320, 196032411, 110246852, 177653359, 59784237, 172142379, 26995670, 148050452, 100895731, 335466064, 243933630, 94818487, 31625800, 30454821, 72693141, 266295648, 267588625, 210990478, 16395778, 77513967, 311406815, 108086818, 321323528, 35799614, 285994753, 249203980, 35384853, 289015671, 230738618, 192593418, 287151056, 216689760, 56822008, 62464536, 119656216, 260601961, 328383387, 11500345, 18442281, 189525784, 150313064, 150000422, 75227501, 45325299, 110572751, 282626860, 217090920, 259538641, 115878513, 236397389, 106367442, 153329047, 164430150, 78937339, 26153204, 61209568, 205015150, 93954178, 284105568, 271121870, 245714780, 130720877, 332723927, 305014408, 174786387, 138445919, 129571160, 139460671, 108179444, 78372938, 173653566, 101030874, 144598141, 249983601, 278649120, 151123069, 89079839, 319401782, 198614650, 307790514, 29827810, 197577096, 93438475, 308844489, 24071278, 331322255, 83543944, 29166009, 200776321, 195340623, 166139611, 209824359, 170952451, 314541843, 18081709, 107045904, 229569764, 276119699, 274732966, 80787704, 331941013, 182653006, 4056075, 59868531, 3747138, 241587548, 224142426, 318548313, 155196446, 311109909, 211111585, 8992154, 341760495, 215419448, 130520959, 89272762, 250458470, 244008696, 212744151, 320005458, 210855928, 244838641, 156877594, 236172081, 146132267, 110666897, 161195810, 94780129, 133825196, 3001989, 7564754, 232572072, 175950821, 207519028, 30307560, 150850261, 155282346, 3480842, 229818257, 233020212, 222328889, 327350332, 102439076, 197043540, 124631236, 240575411, 173553533, 282317591, 185173139, 339194346, 55879596, 27839540, 337101324, 206180361, 277902000, 270234889, 309386475, 306941549, 325519299, 85400974, 183014809, 20419940, 180265490, 104786872, 10317133, 182448169, 17663878, 311237389, 282968089, 264330277, 29524083, 261067087, 273889415, 285315640, 243987075, 113736676, 301465566, 171526384, 102741603, 310182749, 165810823, 230082231, 12249595, 289425419, 146286504, 48480180, 249953531, 39858995, 171424714, 208985206, 147354131, 242810795, 136074818, 193939062, 201586946, 307419103, 312437470, 215325525, 35705608, 185344264, 144839130, 329407110, 102206178, 8906587, 149803954, 310201838, 137727874, 294415930, 335295846, 97433410, 278710526, 57347795, 102883008, 86338852, 67001637, 289103189, 76630541, 147125380, 147125380, 225845255, 38836732, 638481, 95436245, 159466961, 273592911, 282139798, 45977275, 41907821, 207225671, 67827767, 282980076, 39681991, 71376015, 162105952, 238369705, 31163542, 125462113, 102488497, 24356725, 74961035, 73199586, 83383981, 141227286, 176344268, 307033851, 79333790, 238413627, 214264670, 65195649, 282291555, 70944525, 224844976, 164298615, 221526055, 140963330, 165162733, 51118139, 268600124, 123917922, 202916618, 37223198, 212252888, 286609619, 52492441, 111221548, 84803865, 15237202, 129586104, 103793327, 226288879, 148886891, 195911170, 14181805, 26213224, 325512982, 194497710, 62753046, 61792518, 178431421, 65016969, 79227727, 8677542, 203876785, 273785670, 155651656, 261187502, 31057226, 29842982, 334186632, 24055118, 28464460, 216126325, 258990655, 110306513, 27569250, 98061196, 125490567, 257440664, 232323177, 144853164, 297674731, 203090695, 238790729, 278738425, 76420484, 140297749, 293156787, 171437338, 120336046, 59723480, 118834560, 17487241, 33116653, 39610287, 15201241, 234630881, 76160387, 275787575, 40849780, 197399479, 185059624, 250943982, 63741688, 202007070, 326680994, 111528705, 97884641, 279403057, 218366812, 232843753, 337947387, 2198667, 48174752, 23100869, 272319852, 35697474, 147357374, 99768409, 316708648, 248390342, 129635759, 99011730, 93634172, 98625770, 88086146, 215380757, 228946173, 7714272, 249342805, 222423551, 186461929, 246334249, 86575154, 212785572, 74965591, 29522432, 72038157, 303923311, 122617397, 283538401, 52572198, 125768929, 60887290, 52516026, 250949657, 106546372, 180978007, 112118494, 239518396, 124234468, 145433523, 41813793, 268575606, 165542490, 203059446, 107124706, 337707503, 48963391, 229343249, 61595964, 284996233, 25059728, 50640711, 208563579, 337886839, 206008103, 222666367, 228292664, 248025225, 189504485, 65617345, 190428081, 316673975, 17494788, 282857994, 53436939, 249230481, 74390584, 212835331, 137680310, 95613286, 66325358, 309362127, 151421060, 209593179, 114113715, 43423887, 191832954, 238667390, 322626116, 29920868, 106156502, 203955069, 100721248, 82459356, 244991333, 268034519, 316405779, 273735559, 134343451, 166324912, 66883796, 5243503, 29018791, 21898554, 312436143, 308156686, 161783925, 200186101, 324889452, 276349028, 57703064, 95386270, 198736317, 128511569, 186993691, 207709224, 134496539, 14350162, 116747391, 133292896, 241577100, 153110400, 154156067, 297885938, 234490040, 132606731, 99630365, 104702682, 114803506, 76773345, 143091004, 327673317, 127141842, 37365801, 258209625, 131195249, 160254296, 94881513, 97165924, 238146216, 71239909, 188705979, 128952317, 322092829, 126232952, 216138068, 65370229, 217779406, 127294484, 278842410, 213320402, 72702569, 118810993, 24219023, 9679266, 278161791, 99489129, 212282769, 219237152, 27207215, 75251836, 268963225, 282154332, 191048288, 16712403, 59186885, 225367572, 11021192, 156103903, 290290746, 54597981, 154658980, 85537687, 262742464, 294329603, 221642782, 98979886, 55606742, 45903704, 156818490, 287064771, 79013953, 88836726, 207423049, 173922024, 147323990, 65989289, 295671864, 182915287, 149762436, 309679161, 15259334, 208031554, 13331295, 194524078, 97710260, 101773250, 339700451, 219430398, 231453364, 49501804, 279744171, 323892223, 69399046, 89928104, 114561358, 70506940, 328963989, 193711471, 97092453, 16459066, 106430322, 35697648, 169261344, 25220188, 61573089, 191844840, 18341958, 157234161, 309434896, 192185398, 223325944, 43501470, 238265116, 106048120, 203672079, 222559535, 104257802, 283348053, 283173809, 106751962, 76943303, 268031356, 99820257, 94413413, 329088816, 296670590, 197258941, 205824376, 151184563, 311373495, 194841984, 139338522, 34781589, 242307449, 342068358, 52827034, 122047346, 152955167, 103307647, 53434797, 74433958, 45961184, 54520751, 202825365, 256141175, 278835187, 83881315, 140271708, 143300957, 299116442, 3737511, 309731307, 68745677, 39431990, 158927839, 6172407, 205849662, 183641743, 45498735, 98176600, 323672612, 7241400, 32720971, 64125331, 66974490, 104357642, 58050494, 123373985, 120788291, 251573779, 67972151, 164893064, 96668162, 223282352, 67861160, 72901059, 316850298, 298798984, 181839841, 49185546, 219896163, 212142235, 315127530, 120183556, 328164262, 62573436, 162094733, 169230315, 203982821, 19310510, 81974504, 301075120, 336881376, 191922692, 145871047, 294103691, 193847081, 73227217, 30934843, 35634182, 195421831, 58880538, 6106529, 154347260, 94302372, 8365871, 283127354, 141266949, 184797503, 312618828, 232878579, 60636883, 305043982, 54213848, 132737665, 47604320, 152722488, 15620335, 95717621, 110955714, 307683980, 191709016, 241778797, 112521033, 334061218, 8484640, 199116169, 122001951, 144641665, 70661250, 122885767, 252041542, 5914181, 219177478, 9524807, 323851588, 111313207, 36623599, 57324186, 305502406, 32348723, 140756060, 256951970, 329762985, 87724458, 295664397, 313482504, 134890104, 211691553, 63896246, 148832513, 47819188, 174218896, 81678305, 227029337, 67633610, 141574319, 248929416, 205838, 281000597, 135760711, 186098064, 212452636, 330079483, 100711206, 233728729, 256746755, 125997936, 325927085, 108003948, 261776244, 134855544, 232775612, 248365292, 73266486, 269228080, 303009890, 252929541, 311257753, 301952299, 301435772, 134667909, 34306707, 1075386, 111227285, 292131475, 90934101, 157087127, 219338382, 217026336, 339882594, 145810294, 201586874, 292056070, 237032553, 74473802, 157075231};


static const int32_t gj[256] = {172048372, 309870606, 200431766, 324059659, 14516260, 71200620, 112658784, 331539029, 298230187, 199465001, 95054535, 46153749, 68257019, 34453836, 264967053, 333985274, 2157599, 336914091, 97363940, 244886526, 273585413, 313646880, 304604206, 154887261, 259345227, 70957183, 256585624, 125430881, 104908181, 38821441, 203436199, 263147140, 306069179, 219236223, 120513729, 144095014, 262687761, 102505389, 80715937, 16753350, 97015745, 309468163, 207480694, 283016648, 256009818, 175603941, 309005608, 106291750, 92762708, 129369847, 225846253, 182195945, 120488377, 100183661, 144843084, 54014060, 328793523, 80692520, 188489502, 311644197, 44711423, 101476595, 97132319, 66787068, 341297933, 336126437, 215466666, 304156831, 192691578, 235006678, 266568758, 267029920, 230730845, 143640880, 58600117, 62618402, 44164717, 89460722, 328674715, 328313139, 61987110, 277923611, 167341181, 209325648, 195893203, 163176989, 92082345, 101291223, 141586883, 116097341, 119707476, 270409387, 52968398, 301302769, 23858345, 325103190, 269720605, 305570320, 236766448, 179240984, 12957952, 59734849, 225980378, 226444843, 38891102, 153376031, 292044590, 234984074, 211967562, 220775631, 182996330, 275058303, 139026960, 56954705, 137934911, 252295904, 44896463, 92821217, 329737605, 104245071, 263574185, 328595292, 156039906, 248376507, 245797804, 253708871, 70747422, 308520011, 184007114, 56390934, 92905498, 114780906, 87424978, 200016287, 89254003, 48454192, 37043356, 119331054, 7072757, 334740528, 20594664, 303758926, 227847682, 291682234, 36206659, 59506067, 224003799, 202581571, 151922543, 30359439, 327443646, 326283837, 170077366, 237697179, 133577520, 172821269, 203788155, 217317485, 13052895, 15599810, 260374244, 293351741, 274875394, 263545324, 333967048, 183332133, 180664712, 76840946, 156688276, 172483421, 268253444, 212224416, 172063915, 76494101, 315790590, 162553334, 278504038, 308113789, 290150331, 217965216, 169571794, 322975554, 309484036, 104099734, 289946488, 310642530, 184778640, 50537885, 329131026, 270750726, 57066271, 4755224, 172127038, 189147643, 227373088, 273708648, 190961717, 21821056, 290599527, 225963915, 107820937, 17590914, 264117725, 301165510, 36060712, 147514506, 208663554, 29607387, 112759694, 217306643, 333540791, 117231636, 261889302, 50717815, 296093912, 204194230, 100118513, 62196815, 58169614, 156316970, 123297488, 328727631, 173478900, 137254480, 20436325, 18378834, 84216704, 187811389, 326574509, 283127491, 103512701, 108106635, 170244636, 58864560, 166300682, 107363471, 135533785, 304963431, 333212133, 64179189, 183642748, 190254906, 101127339, 332003328, 128576039, 283119421, 141656867, 260949411, 132014089, 12121494, 286669262, 154072039, 220629853, 176201368, 202106993, 181897417};


static const int32_t invgj[256] = {172048372, 161679160, 141469584, 167375209, 122946724, 189504538, 56907315, 331455083, 211562488, 82627166, 201919710, 60457156, 215000538, 11573249, 242449238, 153321671, 159933829, 279397388, 10364444, 38613146, 208042792, 236213106, 177275895, 284712017, 173331941, 235469942, 240063876, 60449086, 17002068, 155765188, 259359873, 325197743, 323140252, 206322097, 170097677, 14848946, 220279089, 187259607, 285406963, 281379762, 243458064, 139382347, 47482665, 292858762, 81687275, 226344941, 10035786, 126269934, 230816883, 313969190, 134913023, 196062071, 307515865, 42411067, 79458852, 325985663, 235755640, 117612662, 52977050, 321755521, 152614860, 69867929, 116203489, 154428934, 171449539, 338821353, 286510306, 72825851, 14445551, 293038692, 158797937, 32934047, 53630089, 239476843, 34092541, 20601023, 174004783, 125611361, 53426246, 35462788, 65072539, 181023243, 27785987, 267082476, 171512662, 131352161, 75323133, 171093156, 186888301, 266735631, 162911865, 160244444, 9609529, 80031253, 68701183, 50224836, 83202333, 327976767, 330523682, 126259092, 139788422, 170755308, 209999057, 105879398, 173499211, 17292740, 16132931, 313217138, 191654034, 140995006, 119572778, 284070510, 307369918, 51894343, 115728895, 39817651, 322981913, 8836049, 336503820, 224245523, 306533221, 295122385, 254322574, 143560290, 256151599, 228795671, 250671079, 287185643, 159569463, 35056566, 272829155, 89867706, 97778773, 95200070, 187536671, 14981285, 80002392, 239331506, 13838972, 250755360, 298680114, 91280673, 205641666, 286621872, 204549617, 68518274, 160580247, 122800946, 131609015, 108592503, 51531987, 190200546, 304685475, 117131734, 117596199, 283841728, 330618625, 164335593, 106810129, 38006257, 73855972, 18473387, 319718232, 42273808, 290608179, 73167190, 223869101, 227479236, 201989694, 242285354, 251494232, 180399588, 147683374, 134250929, 176235396, 65652966, 281589467, 15263438, 14901862, 254115855, 299411860, 280958175, 284976460, 199935697, 112845732, 76546657, 77007819, 108569899, 150884999, 39419746, 128109911, 7450140, 2278644, 276789509, 246444258, 242099982, 298865154, 31932380, 155087075, 262884057, 14783054, 289562517, 198733493, 243392916, 223088200, 161380632, 117730324, 214206730, 250813869, 237284827, 34570969, 167972636, 87566759, 60559929, 136095883, 34108414, 246560832, 326823227, 262860640, 241071188, 80888816, 199481563, 223062848, 124340354, 37507398, 80429437, 140140378, 304755136, 238668396, 218145696, 86990953, 272619394, 84231350, 188689316, 38972371, 29929697, 69991164, 98690051, 246212637, 6662486, 341418978, 9591303, 78609524, 309122741, 275319558, 297422828, 248522042, 144111576, 45346390, 12037548, 230917793, 272375957, 329060317, 19516918, 143144811, 33705971};


void poly_uniform(poly_k a, const unsigned char *seed)         
{ // Generation of polynomials "a_i"
  unsigned int pos=0, i=0, nbytes = (PARAM_Q_LOG+7)/8;
  unsigned int nblocks=PARAM_GEN_A;
  uint32_t val1, val2, val3, val4, mask = (uint32_t)(1<<PARAM_Q_LOG)-1;
  unsigned char buf[SHAKE128_RATE*PARAM_GEN_A];
  uint16_t dmsp=0;

  cshake128_simple(buf, SHAKE128_RATE*PARAM_GEN_A, dmsp++, seed, CRYPTO_RANDOMBYTES);    
     
  while (i < PARAM_K*PARAM_N) {   
    if (pos > SHAKE128_RATE*nblocks - 4*nbytes) {
      nblocks = 1;
      cshake128_simple(buf, SHAKE128_RATE*nblocks, dmsp++, seed, CRYPTO_RANDOMBYTES);    
      pos = 0;
    } 
    val1  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    val2  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    val3  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    val4  = (*(uint32_t*)(buf+pos)) & mask;
    pos += nbytes;
    if (val1 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val1*PARAM_R2_INVN);
    if (val2 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val2*PARAM_R2_INVN);
    if (val3 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val3*PARAM_R2_INVN);
    if (val4 < PARAM_Q && i < PARAM_K*PARAM_N)
      a[i++] = reduce((int64_t)val4*PARAM_R2_INVN);
  }
}


int32_t reduce(int64_t a)
{ // Montgomery reduction
  int64_t u;

  u = (a*PARAM_QINV) & 0xFFFFFFFF;
  u *= PARAM_Q;
  a += u;
  return (int32_t)(a>>32);
}

int64_t barr_reduce64(int64_t a)
{ // Barrett reduction
  int64_t u = (a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - u*PARAM_Q;
}


sdigit_t barr_reduce(sdigit_t a)
{ // Barrett reduction
  digit_t u = ((int64_t)a*PARAM_BARR_MULT)>>PARAM_BARR_DIV;
  return a - (digit_t)u*PARAM_Q;
}


void dgt(poly x)
{
  int32_t a, sub_re, sub_img;
  int i, index, j, m, window;

  window = 1;
  for(m = 512; m >= 2; m >>= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];
      for(i = j; i < PARAM_N; i += (m << 1)) 
      {
        sub_re = x[i] - x[i+m];
        sub_img = x[i+1] - x[i+m+1];
        
        x[i] = x[i] + x[i+m] - PARAM_Q;
        x[i+1] = x[i+1] + x[i+m+1] - PARAM_Q;

        x[i] += (x[i] >> (RADIX32-1)) & PARAM_Q;
        x[i+1] += (x[i+1] >> (RADIX32-1)) & PARAM_Q;
        
        x[i+m] = reduce((int64_t)a * sub_re);
        x[i+m+1] = reduce((int64_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;
  }

}


void idgt(poly x)
{
  int32_t a, mul_re, mul_img;
  int i, index, j, m, window;

  window = 256;
  for(m = 2; m <= 512; m <<= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = invgj[index];
      for(i = j; i < PARAM_N; i += (m << 1)) 
      {
        mul_re = reduce((int64_t)x[i+m] * a);
        mul_img = reduce((int64_t)x[i+m+1] * a);        
	
	x[i+m] = x[i] - mul_re;
        x[i+m] += (x[i+m] >> (RADIX32-1)) & PARAM_Q;

        x[i+m+1] = x[i+1] - mul_img;
        x[i+m+1] += (x[i+m+1] >> (RADIX32-1)) & PARAM_Q;

        x[i] = x[i] + mul_re - PARAM_Q;
        x[i] += (x[i] >> (RADIX32-1)) & PARAM_Q;

        x[i+1] = x[i+1] + mul_img - PARAM_Q;       
        x[i+1] += (x[i+1] >> (RADIX32-1)) & PARAM_Q;
      }
      index += window;
    }
    window >>= 1;
  }
}


static void poly_pointwise(poly result, const poly x, const poly y)
{ // Pointwise polynomial multiplication result = x.y

  for (int i=0; i<PARAM_N; i++)
    result[i] = reduce((int64_t)x[i]*y[i]);
}


void poly_dgt(poly x_dgt, const poly x)
{
  int i, j;

  j = 0;
  for(i = 0; i < PARAM_N; i+=2) {             
      x_dgt[i] = reduce(
        (int64_t)x[j] * nthroots[i] - 
        (int64_t)x[512+j] * nthroots[i+1]
      );
      x_dgt[i+1] = reduce(
        (int64_t)x[j] * nthroots[i+1] + 
        (int64_t)x[512+j] * nthroots[i]
      );
      j++;
  } 

  dgt(x_dgt);
}


void poly_mul(poly result, const poly a, const poly b)
{ /* It is assumed that both signals are already in the DGT domain. 
     The DGT counterpart of poly_b was computed in sign.c. */  
  poly mul;
  int i, j;

  /* Calculating the point-wise multiplication of input signals */
  for(i = 0; i < PARAM_N; i+=2) {             
    mul[i] = reduce(
      (int64_t)a[i] * b[i] -
      (int64_t)a[i+1] * b[i+1]
    );

    mul[i+1] = reduce(
      (int64_t)a[i] * b[i+1] + 
      (int64_t)a[i+1] * b[i]
    );
  }

  /* Recovering the multiplication result in Z[x]/<x^n+1> */
  idgt(mul);

  /* Removing the twisting factors and writing the result from the Gaussian integer to the polynomial form */
  j = 0;
  for(i = 0; i < PARAM_N; i+=2) {
      result[j] = reduce(
               (int64_t)mul[i] * invnthroots[i] -
               (int64_t)mul[i+1] * invnthroots[i+1]);

      result[j+512] = reduce(
               (int64_t)mul[i] * invnthroots[i+1] + 
               (int64_t)mul[i+1] * invnthroots[i]);
      j++;
  }
}


void poly_add(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y

    for(int i=0; i<PARAM_N; i++)
      result[i] = x[i] + y[i];
}


void poly_add_correct(poly result, const poly x, const poly y)
{ // Polynomial addition result = x+y with correction

    for (int i=0; i<PARAM_N; i++) {
      result[i] = x[i] + y[i];
      result[i] += (result[i] >> (RADIX32-1)) & PARAM_Q;    // If result[i] < 0 then add q
      result[i] -= PARAM_Q;
      result[i] += (result[i] >> (RADIX32-1)) & PARAM_Q;    // If result[i] >= q then subtract q
    }
}


void poly_sub(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y

    for (int i=0; i<PARAM_N; i++)
      result[i] = x[i] - y[i];
}    


void poly_sub_reduce(poly result, const poly x, const poly y)
{ // Polynomial subtraction result = x-y

    for (int i=0; i<PARAM_N; i++)
      result[i] = (int32_t)barr_reduce(x[i] - y[i]);
}


/********************************************************************************************
* Name:        sparsemul8
* Description: performs sparse polynomial multiplication
* Parameters:  inputs:
*              - const unsigned char* s: part of the secret key
*              - const uint32_t pos_list[PARAM_H]: list of indices of nonzero elements in c
*              - const int16_t sign_list[PARAM_H]: list of signs of nonzero elements in c
*              outputs:
*              - poly prod: product of 2 polynomials
*
* Note: pos_list[] and sign_list[] contain public information since c is public
*********************************************************************************************/
void sparse_mul8(poly prod, const unsigned char *s, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H])
{
  int i, j, pos;
  int8_t *t = (int8_t*)s;

  for (i=0; i<PARAM_N; i++)
    prod[i] = 0;

  for (i=0; i<PARAM_H; i++) {
    pos = pos_list[i];
    for (j=0; j<pos; j++) {
        prod[j] = prod[j] - sign_list[i]*t[j+PARAM_N-pos];
    }
    for (j=pos; j<PARAM_N; j++) {
        prod[j] = prod[j] + sign_list[i]*t[j-pos];
    }
  }
}


/********************************************************************************************
* Name:        sparsemul32
* Description: performs sparse polynomial multiplication 
* Parameters:  inputs:
*              - const int32_t* pk: part of the public key
*              - const uint32_t pos_list[PARAM_H]: list of indices of nonzero elements in c
*              - const int16_t sign_list[PARAM_H]: list of signs of nonzero elements in c
*              outputs:
*              - poly prod: product of 2 polynomials
*********************************************************************************************/
void sparse_mul32(poly prod, const int32_t *pk, const uint32_t pos_list[PARAM_H], const int16_t sign_list[PARAM_H])
{
  int i, j, pos;
  int64_t temp[PARAM_N] = {0};
  
  for (i=0; i<PARAM_H; i++) {
    pos = pos_list[i];
    for (j=0; j<pos; j++) {
        temp[j] = temp[j] - sign_list[i]*pk[j+PARAM_N-pos];
    }
    for (j=pos; j<PARAM_N; j++) {
        temp[j] = temp[j] + sign_list[i]*pk[j-pos];
    }
  }
  for (i=0; i<PARAM_N; i++)
    prod[i] = (int32_t)barr_reduce64(temp[i]);
}
