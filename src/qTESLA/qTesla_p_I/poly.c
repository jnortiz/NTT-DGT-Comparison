/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: NTT, modular reduction and polynomial functions
**************************************************************************************/

#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"


static const int32_t nthroots[1024] = {172048372, 0, 92194885, 238672960, 178990711, 235148722, 202462144, 30820657, 62607340, 123830012, 71502837, 110369095, 60754183, 653736, 59667725, 221258940, 310813721, 163069688, 173872410, 212038276, 208950413, 44602975, 218366454, 103711520, 86870415, 318271054, 125648848, 23677553, 243865732, 53319855, 270943024, 197419224, 26653162, 69342818, 322015396, 192307168, 17602852, 88718744, 205663327, 264965791, 84823617, 123437202, 48606141, 337720740, 145607952, 317119115, 316148004, 287409325, 72418041, 204991267, 226580841, 187809305, 318722538, 307820569, 129496053, 273720106, 43852040, 23162641, 277393187, 81512031, 216779435, 127046784, 20706942, 273448811, 75322917, 257840912, 228962572, 25514756, 48246709, 200169505, 50192135, 176295900, 55712376, 99282101, 237436074, 75058658, 328831755, 106392182, 204011178, 125863293, 153281095, 117045479, 117621711, 267604239, 44597889, 329475460, 207604895, 42669499, 5870126, 284102237, 245656751, 108906081, 220059646, 338060582, 11126353, 122667794, 98316983, 214006646, 300285822, 179730074, 9323813, 173065299, 244542679, 246930179, 146857685, 45905266, 97494097, 145223082, 285261465, 235740255, 54000181, 25988226, 33068496, 323612353, 9634917, 16996576, 223453636, 239872774, 90401473, 6220050, 334539518, 37672557, 54840268, 336902284, 115395750, 75477828, 81289702, 313564020, 305056652, 195375356, 334278543, 238222024, 213966350, 129562512, 313195793, 277471729, 214956310, 154806056, 99430114, 177619, 101246708, 253238410, 272253465, 329229178, 316690865, 154632654, 24821491, 142497786, 157120342, 147363710, 95933529, 40232482, 319375035, 118340768, 117945756, 331317634, 154603355, 53544841, 201543935, 111332218, 200349623, 26531390, 287105141, 93686205, 160409377, 32624426, 283768666, 306140647, 305866757, 323902756, 306761469, 289489522, 122016411, 332417425, 69710899, 120496603, 183887485, 254514048, 185695644, 214727464, 90756304, 172197960, 45292683, 117064516, 221275405, 320885849, 255660987, 260282603, 45109952, 20698650, 259039557, 249406145, 83279220, 27437106, 186487063, 164062151, 62837951, 322479346, 15554110, 20997771, 136091129, 186848653, 315065315, 69981509, 330157762, 242409851, 360366, 242766489, 226307503, 252517986, 199458211, 26159487, 164009086, 311212324, 311386571, 226733266, 294046420, 38711587, 299837677, 76624292, 73040479, 39339227, 75672839, 194416998, 17555462, 249629287, 172966662, 87144239, 95397910, 50073905, 340596827, 66560634, 185575132, 193725988, 115114872, 203039082, 221765212, 235138695, 66550019, 80534855, 60900990, 31077484, 160788341, 326570387, 276477497, 66553955, 228121840, 311664680, 172652345, 251557936, 43038274, 257301899, 134039934, 248716737, 325761806, 125861014, 312193573, 18174287, 103640206, 337263697, 12316245, 265559954, 304594523, 171454417, 14222436, 231687353, 150296348, 255766948, 125203111, 40104387, 211405712, 291193392, 141577653, 342985147, 194105611, 238761029, 237563067, 72438883, 149906244, 3109882, 153056812, 88746185, 137609156, 213213369, 177780628, 174011010, 323932312, 269841546, 263796961, 82599900, 6659649, 156662837, 16339185, 234781906, 201181433, 50376243, 332896178, 166494212, 125461623, 121454930, 226078906, 21372613, 161860486, 24701648, 304061771, 10904593, 130817533, 229804804, 339297355, 202878385, 158515427, 128633512, 289639484, 333056121, 278924093, 150612099, 131757018, 85213563, 100001310, 61271934, 310737440, 329269386, 14458607, 314589702, 237335474, 83899989, 194957569, 118426187, 260349235, 52362080, 337943670, 12336690, 92116988, 23565018, 214354523, 133803756, 282229292, 39043589, 231563912, 70826808, 130817729, 212440826, 175720867, 193996119, 223147106, 87058537, 209206576, 83417632, 278565793, 326678621, 260587352, 102194828, 210988801, 133984468, 175894211, 102752423, 5653901, 316014192, 155384134, 165456424, 246653786, 294293484, 74680512, 31024891, 4338282, 236022965, 326860158, 288089535, 181436287, 231373413, 134085727, 36303245, 33776334, 228777567, 331698660, 78733116, 89160327, 64552338, 293651020, 155850726, 143614775, 250991092, 187095098, 291134966, 122814929, 125893327, 58895223, 62140605, 48774583, 162044224, 98845692, 222921607, 337575095, 224217671, 302022582, 101097913, 150549075, 135966907, 62497182, 184453363, 79478, 201242082, 218325926, 215804061, 240234420, 209033960, 106357364, 235287401, 232508419, 107849397, 103043371, 189651213, 34934025, 25496861, 212385347, 316174078, 208313538, 129468820, 316692653, 152363759, 127385366, 315545717, 136865954, 60996752, 326011873, 137548121, 283916344, 268376723, 104181869, 228023144, 185810279, 75633365, 316092579, 242216653, 199545795, 154404478, 312411127, 141455018, 208303693, 135362599, 284268281, 240252591, 204021935, 1007481, 147039228, 26916438, 86433665, 225969058, 19442480, 103560361, 17723229, 106797072, 250265209, 265603412, 298639023, 7759262, 185771102, 282824643, 57342620, 93976343, 151774919, 326249997, 22322550, 146020642, 104611938, 36477842, 6824490, 192625101, 230598684, 176654639, 181470933, 160427822, 201658664, 122406791, 106868199, 235000617, 69725945, 171218832, 61288763, 131908460, 207070634, 251341270, 87098370, 219888496, 99791825, 47850584, 63037447, 297374369, 216176213, 67646260, 9288372, 5224785, 108950361, 39309089, 28887895, 275725717, 15696156, 159759911, 312047632, 15429619, 196769568, 293706646, 279979711, 81892956, 299877240, 125001443, 318562857, 191096656, 332859045, 287185205, 289106315, 187047012, 258652380, 84924197, 27997942, 117267166, 246825063, 267807009, 10974878, 143677303, 267386821, 234821352, 131583445, 296789740, 13921969, 300571086, 83284261, 224893959, 296991034, 118165869, 266948204, 304170607, 145838126, 127373883, 310318184, 78634514, 209287432, 156367439, 62489022, 244185642, 30585813, 51971445, 294964611, 47849832, 3947476, 250278849, 247498975, 151834061, 5896162, 342723184, 101692277, 102464751, 172355436, 98223872, 156019438, 189051088, 230605812, 286586969, 296303461, 251794994, 284120137, 254991360, 300283653, 77849938, 22525699, 279074300, 222986688, 148263934, 257681155, 342908450, 276156615, 273975206, 50245498, 295677552, 270002292, 182011412, 20555645, 190903562, 114268113, 223414521, 143172036, 300911556, 249349415, 2612193, 323700383, 33839302, 297976721, 91905769, 149322163, 208685538, 177192579, 103682978, 222852771, 335342147, 230191275, 15662092, 227961340, 117435719, 141408875, 276699335, 310864762, 239133448, 77409273, 99166229, 316838478, 316417170, 43560153, 191139196, 146835988, 69471289, 100139690, 195472462, 304691166, 182583084, 83234575, 193650412, 139174797, 337352688, 343102256, 264482136, 169148178, 341284227, 89877383, 243212441, 240847831, 179513852, 23072115, 21320495, 161631508, 284245008, 177824916, 59408571, 182267060, 327990345, 274967318, 327632041, 191395385, 147180785, 15868397, 10260485, 107303605, 31457470, 242226014, 237565249, 248547447, 107528001, 88387049, 13252895, 215174436, 304488108, 301564369, 143753349, 199065167, 128770921, 4296993, 247296056, 98473081, 145213576, 137390416, 131312971, 126657908, 280380315, 254781980, 44164132, 138629476, 121803213, 80506549, 248733384, 332999638, 141283295, 304395624, 35132036, 38933274, 50741740, 316181647, 165680378, 17264270, 278961935, 260434930, 63650698, 65836705, 153908531, 196180477, 141334398, 140004325, 320846752, 105718527, 131066150, 263701235, 332408404, 294422677, 202593463, 211282334, 216795966, 290145207, 233744678, 216862575, 187793681, 167833031, 183133882, 97562307, 167801401, 88768806, 147357626, 283725651, 331522553, 267985760, 104446743, 20153236, 281023459, 187457897, 217356275, 76627644, 113250086, 102459222, 65132575, 135722773, 190801576, 181819150, 64659911, 12060572, 87796288, 182206001, 170847100, 206154834, 240982990, 106858215, 82509102, 250503445, 227288535, 308200724, 165503188, 325350949, 117232128, 295372927, 157474764, 58277074, 246051719, 184389167, 181165133, 238709889, 73358518, 95835007, 326367189, 157887236, 328803025, 260587098, 67355544, 305487009, 13924628, 245497389, 311725941, 136664459, 99033993, 203242278, 143478369, 98180248, 304798787, 210362363, 252742128, 170317819, 66015875, 303059161, 221633540, 318270896, 24828740, 282014090, 315610861, 291757248, 55681685, 249554273, 128097931, 311168535, 116938740, 255549499, 309808074, 132321658, 263063541, 81233496, 239302052, 306469375, 325511215, 134396363, 258561040, 299730998, 83630441, 17732164, 32648079, 13672719, 167779621, 136247310, 224838051, 22195355, 23348609, 92179684, 40948671, 336963660, 73440398, 31646372, 64752579, 328177944, 210813419, 154210294, 58716679, 969547, 9413604, 119953223, 65356386, 228533222, 98588615, 234705234, 10852212, 111179910, 21585097, 224062211, 43942408, 331161028, 110648584, 175847726, 129777574, 180763471, 304485635, 30689711, 88051469, 271647941, 261073706, 298111408, 168333907, 186720164, 100251774, 284680987, 78900841, 219908535, 203325867, 97255654, 174971677, 9704966, 191994637, 125726132, 208381117, 163068035, 24704477, 37877891, 298273404, 259510693, 272098853, 256509928, 311721754, 22193445, 270078587, 242909137, 128164244, 75994813, 275764467, 24815373, 312534046, 337332603, 53580223, 247595163, 170675151, 201488330, 248974466, 62484938, 299556500, 335786479, 328684891, 12833410, 7921817, 314435974, 99379943, 185788760, 172589121, 163914363, 143218971, 296357242, 202176199, 202133783, 308538693, 295081061, 50433836, 74708560, 141712797, 83965918, 197778252, 181719782, 57748679, 10838905, 165483713, 282198205, 123724751, 39149115, 238401931, 185919913, 328505208, 270833097, 110964448, 219248442, 50639599, 126821503, 55649465, 45650633, 25267923, 92687633, 165822344, 11818502, 123282469, 309397428, 114087646, 50467894, 229724433, 22300358, 204469304, 129618462, 180452465, 106814489, 227377711, 151937741, 240786138, 301518898, 4215675, 194881311, 119318846, 12496698, 94131086, 325704136, 310737960, 175519905, 70990612, 335900970, 313900226, 245801624, 95108333, 117376443, 258823189, 250932529, 178109708, 343349510, 203008458, 329350764, 301262876, 201888897, 218648897, 309395965, 180639834, 276271908, 210159993, 41612084, 39186589, 12651137, 117879529, 308611081, 72039806, 246698079, 281629516, 221520931, 174253823, 331734685, 133230137, 253559267, 172049628, 16882973, 65306166, 109634619, 340435360, 164913047, 46599604, 35777129, 183202004, 59702676, 333472767, 321448297, 148441512, 247888690, 219931587, 185364534, 87110604, 190804766, 112633666, 289550909, 297306265, 330497276, 203697178, 142173714, 241254243, 223610237, 98428398, 279159595, 129889324, 25087960, 236044819, 83761783, 210437080, 272967305, 144465804, 72459811, 59760941, 326341715, 16995642, 263820575, 99845623, 278424126, 285032719, 339572135, 306608222, 8911127, 308657936, 161506867, 315776465, 287975121, 160746244, 182469984, 235001875, 281325293, 97636889, 16074869, 215419686, 224156102, 85956812, 216906754, 335395900, 199406071, 336215415, 246636475, 262137099};

static const int32_t invnthroots[1024] = {8388608, 0, 164918747, 296792553, 92590247, 337147680, 303968767, 918450, 206850669, 304217988, 146709206, 78481220, 326991074, 158488913, 331285653, 235852510, 169418054, 139686581, 32156007, 33908005, 161654368, 187204988, 281096909, 283190078, 162279273, 41061179, 6234443, 192075502, 132900699, 55059598, 141707849, 23345157, 78123728, 91400437, 27252930, 79691115, 331287690, 316685657, 314304155, 200098119, 12271107, 255232550, 105825736, 11801181, 326527175, 168458594, 241486909, 40368399, 127719108, 191547063, 35064634, 35874552, 41363482, 204856531, 59342149, 241620500, 1322362, 270986786, 201672217, 70390163, 51090663, 272794441, 65085521, 39377702, 307467536, 191215706, 196282051, 323620952, 305587056, 104706617, 172799676, 23054022, 121338698, 63267728, 87376943, 5436676, 187452626, 86540483, 238298577, 206601511, 91001949, 195406423, 283535069, 340959144, 128597217, 86170878, 325375637, 71829921, 36633093, 175143972, 164083583, 204179538, 329990083, 297045014, 184723959, 101519215, 214006351, 199316248, 158505982, 107695916, 317341567, 305361748, 76683322, 210684665, 130416357, 192210149, 89257618, 32963496, 197758398, 137268087, 312481418, 16567578, 268100599, 19878279, 306397243, 132823949, 160829153, 36138022, 194826748, 291972637, 318988587, 4003205, 80849632, 183015074, 202034800, 134791487, 219541387, 84975398, 183295010, 210280853, 193478552, 306139966, 91904140, 284832276, 78978245, 39515369, 33123002, 126276904, 214387523, 252963927, 210150815, 153986119, 121174926, 320596949, 66710535, 53537925, 142360680, 67199518, 126088591, 186827600, 210432900, 253077322, 153949717, 82218759, 300295545, 230477644, 188758590, 261094588, 240921705, 87211175, 111364883, 331512928, 244846544, 307217944, 256826615, 306946453, 33214702, 103528856, 43722769, 253668340, 324055555, 180463445, 181050515, 126679645, 72723505, 311562889, 261565170, 196573718, 154949872, 326970429, 319279244, 282004352, 260949190, 87833308, 23534931, 87588796, 217155499, 292331370, 163439654, 4007333, 152669638, 154151088, 57436201, 294160564, 273270640, 189350567, 259220333, 281475473, 288208370, 289310488, 144114182, 131665894, 163236640, 289161748, 170369843, 224448027, 60647792, 202313044, 166636015, 302666898, 338294018, 173363811, 208738087, 89703284, 93967916, 226355816, 38442292, 325682975, 277941522, 219198413, 18807730, 222114995, 338993922, 336564905, 14503753, 187252468, 41060398, 109990063, 248431199, 50341516, 43698098, 152818906, 85268699, 103969090, 62846601, 9128564, 147287207, 181827304, 75892190, 45596631, 101491589, 35651180, 333344811, 307748540, 312673085, 342306956, 62203803, 58893565, 7895324, 253397704, 36841686, 39092712, 8280972, 14826353, 250551096, 236392534, 243699178, 128271378, 185154628, 263842104, 126205518, 289271114, 341325263, 145025504, 255798226, 339888685, 11230405, 339810473, 299211671, 101807538, 278765152, 196891155, 80048138, 264797038, 93884512, 29717673, 158394525, 107442074, 69920546, 59885361, 203969738, 174164107, 29492500, 256824207, 221589119, 86098962, 292259719, 31179121, 81677378, 140806258, 35202064, 84646196, 228969, 24193325, 201637649, 184607294, 27956891, 185391032, 141081231, 154803381, 34023070, 54623290, 88912022, 32525481, 300800982, 276448221, 38376025, 180827845, 231213172, 185615213, 322901301, 229298301, 146509656, 126007361, 9148148, 241882200, 20009266, 271064031, 156558182, 235685487, 5344841, 226931121, 47261168, 197114736, 305654578, 130663809, 219119332, 304999042, 160747257, 245851057, 146744997, 294023383, 130606742, 277401212, 221170483, 99411108, 287857780, 162549623, 146544453, 170910759, 65047262, 325184835, 87512284, 341180268, 137865443, 107910141, 79307982, 333673243, 162267416, 95255097, 176138171, 95342322, 13520065, 189159495, 255066858, 240982758, 22420528, 314742171, 209852784, 121302449, 34461345, 238622328, 319505107, 18291723, 68323466, 249382478, 121716159, 225203788, 252597669, 81004478, 306480889, 193924480, 101579754, 52623993, 254607961, 311448328, 228751024, 195525754, 151011685, 327969842, 128380184, 79369268, 128634163, 234393703, 216138897, 121250111, 340241377, 330795672, 182141471, 37875917, 230809371, 186014299, 295230679, 120432650, 218206633, 98967790, 333394975, 184893562, 140290674, 248913746, 3034805, 315593203, 241381257, 262204227, 84547571, 183536785, 113252631, 203328471, 11395660, 266134231, 185587472, 181691441, 142618784, 320165845, 291616615, 113882575, 305462527, 13134171, 88951655, 145532338, 246221571, 51722919, 329678253, 338728051, 166886962, 253048987, 210578453, 71525947, 124373247, 340447148, 335554591, 114299616, 299271326, 108945565, 63281018, 87561209, 20539029, 291614236, 187401896, 268508263, 210775165, 106064405, 20807590, 197472150, 40850598, 303034064, 176921981, 311814141, 165450667, 255629141, 330511108, 77314220, 268996696, 253558009, 60929428, 132999184, 321430688, 173298149, 1631674, 300193984, 3229259, 173757437, 116914403, 174557037, 498030, 48431582, 310515964, 238985419, 284413045, 335073599, 204367833, 159404699, 171980132, 243924841, 183396231, 250102287, 269088582, 323433621, 34519999, 149160308, 42764850, 314042755, 102092753, 259119474, 93377179, 314661775, 315869486, 299165360, 290869190, 176747908, 331651296, 69182902, 304233522, 30583368, 269684332, 71951801, 98874849, 210129012, 339318487, 196454400, 34139453, 118748305, 185788916, 314464514, 188352076, 303462505, 92214196, 106004149, 150166765, 68635864, 194832958, 165023125, 147125380, 147125380, 276837103, 264499300, 136783653, 110072805, 118477684, 26890775, 191492824, 251728351, 282000108, 127623333, 71704734, 192877509, 8753760, 10127300, 317717736, 190546976, 72340663, 230113043, 149720480, 231298767, 117443605, 120770499, 137025914, 56616861, 236118658, 4574216, 284617811, 311842415, 97060382, 1171982, 28674917, 285462014, 99572738, 211260416, 75491788, 166954769, 144734306, 326591650, 81435884, 198236344, 141233416, 142578789, 140594061, 18338956, 34599669, 28170687, 73886526, 65558384, 149257853, 250928354, 335490160, 133242117, 156537732, 104571491, 41486320, 62715678, 65777859, 116850121, 115939013, 80708009, 333048397, 234161137, 230371731, 203960619, 63519858, 258184665, 285919020, 340605199, 90593448, 117034920, 226612420, 166535741, 317670085, 51263853, 229774696, 336926967, 259326100, 302252911, 206484816, 108763456, 200791074, 240543322, 102444654, 72940752, 245456694, 37024164, 112333665, 322137347, 247064799, 108442460, 265680261, 186302545, 183493690, 115472764, 72726053, 129776450, 173747912, 44545401, 152377695, 177759490, 163435077, 200442096, 222327532, 1559026, 101787794, 77633839, 188301710, 292107721, 238643533, 282755839, 20524491, 79183509, 180872173, 277691808, 218356161, 55402943, 47841865, 293329137, 296370094, 122142458, 81632202, 340699327, 129157709, 294494810, 216869873, 251527971, 204915526, 311797448, 7279103, 38555364, 8542793, 208415431, 17751644, 50426111, 73941359, 90417339, 131679184, 324139384, 280273848, 138169919, 77241425, 63487674, 276252630, 256657017, 44180795, 255701938, 222577871, 71122615, 175204124, 128695356, 13995752, 231700863, 317057814, 218458164, 269036416, 266395014, 231712249, 130510817, 201576089, 172047249, 249829456, 230331553, 281807157, 279282943, 322265966, 203589986, 273286573, 65326870, 113786011, 280825908, 176229500, 129256761, 134347934, 90810255, 340297594, 28974883, 185470584, 60646707, 197334138, 340712469, 310048272, 78523618, 107469950, 205503242, 293479278, 172078560, 286701364, 290771287, 340164722, 160352233, 251615057, 20866627, 87355912, 19936125, 89415817, 317819415, 255963453, 212177443, 263030269, 39697178, 231091750, 152018296, 53409039, 49665910, 169552934, 341307930, 332190059, 221523018, 204046838, 261392587, 211421864, 124373368, 116999569, 79609670, 341217569, 290584645, 309451520, 166698020, 203115055, 161690656, 308317015, 43605113, 24319056, 150470572, 251498987, 101366616, 256680202, 270756165, 219849129, 303044929, 208869261, 29227173, 4032362, 217126767, 277284307, 125021986, 297740595, 178790704, 46972205, 292971458, 311935010, 268006299, 41012257, 111820479, 62235961, 190284085, 219885211, 67077022, 31874128, 60470457, 74333952, 57015025, 289880407, 95086394, 118810993, 24219023, 243165199, 180546707, 42761799, 256078540, 79015159, 44205109, 311857597, 209700813, 317343377, 160605970, 277272813, 37709741, 73782065, 249986865, 271164091, 84433101, 124301175, 173671452, 253444352, 61303281, 101053761, 252089215, 319797220, 147268109, 261838722, 60400140, 139004737, 14576732, 140419236, 175476751, 274946190, 3992000, 103721113, 176337826, 196693899, 63606903, 61886073, 158452756, 148377217, 141687867, 203099335, 307402856, 148238397, 84231718, 86616285, 280779546, 278307072, 288779682, 241380386, 288549937, 89722841, 129538473, 125622659, 323500823, 77535459, 338613438, 250341916, 41574597, 181812803, 128047438, 38570170, 186858160, 310077765, 249467202, 214551434, 217584664, 326169756, 128753111, 186388661, 296774927, 284480035, 180079734, 126385666, 93187210, 295597445, 139400726, 145365758, 276109090, 86391242, 83521845, 25064121, 34087294, 75135581, 103774127, 267642514, 169176001, 148263183, 2757844, 17374131, 231628373, 143667992, 193619573, 45143197, 77951864, 302023420, 294198765, 90137878, 170723600, 79959292, 232551697, 85199153, 276241414, 236440031, 124862197, 320839524, 59535902, 87524060, 229191541, 4304604, 163016185, 33854417, 289945248, 90563458, 16244455, 251466559, 163537965, 30197507, 151462649, 145248724, 186131508, 10607657, 182584394, 75410428, 8976769, 209832254, 267095265, 174854073, 265139195, 94559150, 325970560, 126304442, 111168586, 250287869, 308574971, 270505924, 169792795, 159721573, 129335699, 168230502, 302877264, 193295021, 67757030, 214696368, 311301686, 256391095, 35460075, 280287447, 102113193, 57993767, 54164470, 95378475, 142646393, 170257430, 143797698, 288217611, 196598854, 240586220, 255753841, 93022606, 122609759, 159949010, 54333157, 146948739, 341804677, 278697628, 116953602, 281724000, 30856743, 46385651, 106291154, 83853459, 86478087, 162245233, 310465495, 188793093, 47345031, 88153115, 113679701, 264600710, 178527567, 20278039, 113614418, 222310799, 80416947, 24502056, 176387593, 321822949, 293153745, 169153930, 179393672, 161555116, 24681661, 100520232, 84511605, 257930571, 114325811, 279986220, 194733185, 160425712, 177071024, 144138882, 335942127, 327401589, 328191016, 327167192, 241805788, 148701982, 166949511, 228046619, 195328540, 250920278, 182240666, 62396027, 298521427, 245844657, 43452449, 237397456, 63347874, 69962271, 195911636, 21849136, 109422936, 277949308, 328761464, 199015793, 203469407, 226247317, 88102036, 267794398, 311120865, 291856456, 267578484, 150517314, 336439600, 279914132, 94209662, 255412378, 188895942, 177475168, 16169142, 217851699, 180395374, 59053501, 176366964, 124359445, 336055395, 174714336, 330033337, 138967133, 300234071, 52801019, 79505121, 301095663, 46122244};

static int32_t gj[256] = {172048372, 184007114, 341297933, 172127038, 306069179, 260374244, 269720605, 20436325, 2157599, 36206659, 61987110, 112759694, 92762708, 278504038, 139026960, 183642748, 298230187, 37043356, 230730845, 107820937, 97015745, 156688276, 38891102, 170244636, 259345227, 170077366, 141586883, 100118513, 328793523, 289946488, 263574185, 132014089, 14516260, 87424978, 192691578, 190961717, 262687761, 333967048, 12957952, 326574509, 273585413, 151922543, 195893203, 261889302, 120488377, 169571794, 44896463, 128576039, 68257019, 20594664, 44164717, 36060712, 256009818, 172063915, 211967562, 135533785, 104908181, 203788155, 52968398, 123297488, 44711423, 329131026, 245797804, 220629853, 200431766, 92905498, 215466666, 227373088, 120513729, 274875394, 236766448, 84216704, 97363940, 224003799, 167341181, 333540791, 225846253, 290150331, 137934911, 101127339, 95054535, 7072757, 58600117, 264117725, 207480694, 268253444, 292044590, 166300682, 256585624, 133577520, 119707476, 58169614, 188489502, 184778640, 156039906, 286669262, 112658784, 89254003, 266568758, 290599527, 80715937, 180664712, 225980378, 103512701, 304604206, 327443646, 92082345, 296093912, 144843084, 309484036, 329737605, 141656867, 264967053, 227847682, 328674715, 208663554, 309005608, 315790590, 182996330, 333212133, 203436199, 13052895, 23858345, 173478900, 97132319, 57066271, 70747422, 202106993, 309870606, 56390934, 336126437, 189147643, 219236223, 293351741, 305570320, 18378834, 336914091, 59506067, 277923611, 217306643, 129369847, 308113789, 56954705, 190254906, 199465001, 119331054, 143640880, 17590914, 309468163, 172483421, 153376031, 58864560, 70957183, 237697179, 116097341, 62196815, 80692520, 310642530, 328595292, 12121494, 71200620, 200016287, 235006678, 21821056, 102505389, 183332133, 59734849, 283127491, 313646880, 30359439, 163176989, 50717815, 100183661, 322975554, 92821217, 283119421, 34453836, 303758926, 89460722, 147514506, 175603941, 76494101, 220775631, 304963431, 38821441, 217317485, 301302769, 328727631, 101476595, 270750726, 253708871, 176201368, 324059659, 114780906, 304156831, 273708648, 144095014, 263545324, 179240984, 187811389, 244886526, 202581571, 209325648, 117231636, 182195945, 217965216, 252295904, 332003328, 46153749, 334740528, 62618402, 301165510, 283016648, 212224416, 234984074, 107363471, 125430881, 172821269, 270409387, 156316970, 311644197, 50537885, 248376507, 154072039, 331539029, 48454192, 267029920, 225963915, 16753350, 76840946, 226444843, 108106635, 154887261, 326283837, 101291223, 204194230, 54014060, 104099734, 104245071, 260949411, 333985274, 291682234, 328313139, 29607387, 106291750, 162553334, 275058303, 64179189, 263147140, 15599810, 325103190, 137254480, 66787068, 4755224, 308520011, 181897417};

static int32_t invgj[256] = {172048372, 159569463, 171449539, 2278644, 323140252, 73855972, 83202333, 37507398, 159933829, 204549617, 65072539, 250813869, 230816883, 281589467, 307369918, 341418978, 211562488, 80002392, 53630089, 14783054, 243458064, 201989694, 173499211, 84231350, 173331941, 304685475, 186888301, 246560832, 235755640, 112845732, 306533221, 45346390, 122946724, 97778773, 14445551, 298865154, 220279089, 290608179, 139788422, 238668396, 208042792, 131609015, 171512662, 87566759, 307515865, 299411860, 322981913, 275319558, 215000538, 298680114, 174004783, 223088200, 81687275, 147683374, 191654034, 69991164, 17002068, 330618625, 9609529, 80888816, 152614860, 150884999, 256151599, 329060317, 141469584, 272829155, 286510306, 246444258, 170097677, 319718232, 330523682, 140140378, 10364444, 160580247, 27785987, 34570969, 134913023, 14901862, 115728895, 78609524, 201919710, 13838972, 34092541, 198733493, 47482665, 251494232, 16132931, 38972371, 240063876, 117596199, 162911865, 262860640, 52977050, 77007819, 254322574, 230917793, 56907315, 187536671, 158797937, 155087075, 285406963, 223869101, 209999057, 86990953, 177275895, 51531987, 75323133, 136095883, 79458852, 284976460, 336503820, 248522042, 242449238, 205641666, 53426246, 117730324, 10035786, 176235396, 119572778, 246212637, 259359873, 106810129, 68701183, 223062848, 116203489, 128109911, 250671079, 143144811, 161679160, 35056566, 338821353, 276789509, 206322097, 18473387, 327976767, 80429437, 279397388, 68518274, 181023243, 237284827, 313969190, 15263438, 51894343, 9591303, 82627166, 239331506, 239476843, 289562517, 139382347, 242285354, 17292740, 188689316, 235469942, 117131734, 266735631, 326823227, 117612662, 76546657, 295122385, 12037548, 189504538, 95200070, 293038692, 31932380, 187259607, 73167190, 170755308, 218145696, 236213106, 108592503, 131352161, 60559929, 42411067, 280958175, 8836049, 297422828, 11573249, 91280673, 125611361, 161380632, 226344941, 134250929, 140995006, 98690051, 155765188, 164335593, 80031253, 199481563, 69867929, 39419746, 228795671, 19516918, 167375209, 89867706, 72825851, 242099982, 14848946, 42273808, 126259092, 304755136, 38613146, 122800946, 267082476, 167972636, 196062071, 254115855, 39817651, 309122741, 60457156, 250755360, 20601023, 243392916, 292858762, 180399588, 313217138, 29929697, 60449086, 283841728, 160244444, 241071188, 321755521, 108569899, 143560290, 272375957, 331455083, 14981285, 32934047, 262884057, 281379762, 227479236, 105879398, 272619394, 284712017, 190200546, 171093156, 34108414, 325985663, 199935697, 224245523, 144111576, 153321671, 286621872, 35462788, 214206730, 126269934, 65652966, 284070510, 6662486, 325197743, 38006257, 50224836, 124340354, 154428934, 7450140, 287185643, 33705971};


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
  int m, distance;
  int j1, j2, j, k;
  int32_t a, temp_re, temp_img;
  
  distance = 512;
  for(m = 1; m < 512; m <<= 1)
  {
    for(k = 0; k < m; ++k)
    {
      j1 = k*distance << 1;
      j2 = j1 + distance;

      a = gj[k];
      for(j = j1; j < j2; j = j+2)
      {
        temp_re = reduce((int64_t)a * x[j+distance]);
        temp_img = reduce((int64_t)a * x[j+distance+1]);

        x[j+distance] = x[j] - temp_re;
        x[j+distance] += (x[j+distance] >> (RADIX32-1)) & PARAM_Q;

        x[j+distance+1] = x[j + 1] - temp_img;
        x[j+distance+1] += (x[j+distance+1] >> (RADIX32-1)) & PARAM_Q;
        
        x[j] = x[j] + temp_re - PARAM_Q;
        x[j] += (x[j] >> (RADIX32-1)) & PARAM_Q;

        x[j+1] = x[j + 1] + temp_img - PARAM_Q;        
        x[j+1] += (x[j+1] >> (RADIX32-1)) & PARAM_Q;
      }
    }
    distance >>= 1;
  }
}


void idgt(poly poly)
{
  int distance, j1, jtwiddle, j;
  int32_t sub_re, sub_img;

  for(distance = 2; distance < PARAM_N; distance <<= 1)
  {
    for(j1 = 0; j1 < distance; j1 += 2)
    {
      jtwiddle = 0;
      for(j = j1; j < PARAM_N; j += distance << 1)
      {
        sub_re = poly[j];
        sub_img = poly[j+1];

        poly[j] = barr_reduce(sub_re + poly[j+distance]);
        poly[j+1] = barr_reduce(sub_img + poly[j+distance+1]);

        poly[j+distance] = reduce((int64_t)invgj[jtwiddle] * (sub_re - poly[j+distance]));
        poly[j+distance+1] = reduce((int64_t)invgj[jtwiddle++] * (sub_img - poly[j+distance+1]));
      }
    }
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

      ++j;
  } 

  dgt(x_dgt);

}


void poly_mul(poly result, const poly x, const poly y)
{ /* It is assumed that both signals are already in the DGT domain. 
     The DGT counterpart of poly_b was computed in sign.c. */
  
  poly aux_mul;
  int i;

  /* Calculating the point-wise multiplication of input signals */
  for(i = 0; i < PARAM_N; i+=2) {             
    aux_mul[i] = reduce(
      (int64_t)x[i] * y[i] -
      (int64_t)x[i+1] * y[i+1]
    );

    aux_mul[i+1] = reduce(
      (int64_t)x[i] * y[i+1] + 
      (int64_t)x[i+1] * y[i]
    );
  }

  /* Recovering the multiplication result in Z[x]/<x^n+1> */
  idgt(aux_mul);

  /* Removing the twisting factors and writing the result from the Gaussian integer to the polynomial form */
  int j = 0;
  for(i = 0; i < PARAM_N; i+=2) {
      result[j] = reduce(
               (int64_t)aux_mul[i] * invnthroots[i] -
               (int64_t)aux_mul[i+1] * invnthroots[i+1]);

      result[j+512] = reduce(
               (int64_t)aux_mul[i] * invnthroots[i+1] + 
               (int64_t)aux_mul[i+1] * invnthroots[i]);
      ++j;
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
* Name:        sparse_mul8
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
* Name:        sparse_mul32
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
