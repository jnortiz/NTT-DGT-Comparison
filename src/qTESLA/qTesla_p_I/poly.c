/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: DGT, modular reduction and polynomial functions
**************************************************************************************/

#include "poly.h"
#include "sha3/fips202.h"
#include "api.h"

static const int32_t nthroots[1024] = {1, 0, 95598679, 230765503, 169905438, 233756460, 260688610, 144613436, 213445814, 341735197, 291047220, 289796272, 134919063, 331758667, 15353298, 58343183, 134378070, 175281391, 270815073, 114094683, 326599814, 15612353, 277281968, 215710474, 128975649, 147041113, 50813004, 141073720, 206107103, 271292966, 214800757, 24059662, 320099168, 326187203, 306934857, 69315579, 35674473, 178225281, 166122319, 201397764, 265895951, 247002466, 177994675, 304854525, 92329853, 321637827, 182825192, 117534032, 306071865, 177303797, 132618415, 89570472, 265219570, 58749261, 39865720, 208781989, 197768825, 312804364, 172698570, 269680462, 206970995, 183276596, 221571946, 42888295, 175619471, 2414777, 55723482, 326679389, 69647694, 101944834, 112546975, 291093134, 118200888, 76745925, 218871295, 265738498, 249571277, 163779570, 149922452, 320794310, 342341930, 89157113, 126445810, 213190377, 165238955, 128706839, 326803989, 211052972, 113967758, 29345559, 66363552, 111564452, 171869886, 187796153, 8553404, 257645467, 94396519, 214506407, 59421891, 198124793, 157229030, 222163738, 301125346, 50088960, 60053229, 17475810, 220843730, 21609737, 73984730, 76050576, 333877669, 26404288, 302288396, 82572671, 91534722, 269316989, 52150534, 29560260, 267006359, 173775083, 215916577, 334574414, 197438487, 269344708, 226217455, 331734165, 4527177, 153465512, 315789837, 273366424, 155972416, 221867879, 70007704, 167187152, 287780447, 151331725, 8105016, 274989843, 95113495, 188613329, 140615771, 225197457, 216754996, 274574060, 189977325, 198322461, 206993311, 84598502, 333050022, 340870937, 108803014, 181490203, 110980062, 308866754, 20551890, 214624358, 251999383, 200468898, 24986136, 173604417, 226685461, 156938953, 133765303, 342801710, 266732680, 325470081, 319490435, 341986892, 323922071, 322312993, 326145059, 192411606, 322688405, 287797213, 337576301, 279105293, 7461062, 311892219, 188871786, 59840778, 221538607, 108886933, 1081430, 291035159, 184230055, 18943858, 279153622, 228968063, 240147161, 158490431, 323372123, 236913231, 321001553, 89256786, 281968409, 148729034, 233783407, 93762877, 88440203, 56557835, 266834205, 68299003, 80128150, 313070854, 244323736, 92509762, 252804288, 268309021, 202102209, 230167200, 290594602, 290857749, 102149309, 134854900, 276145580, 176595100, 213496980, 67657691, 95702627, 66998696, 150637090, 231198011, 209718083, 184840105, 204507423, 229293967, 306890626, 218315799, 326995641, 218690741, 32161323, 79533366, 285942709, 159846277, 254306081, 216322221, 139189083, 128285449, 162941512, 52037686, 84291759, 127741829, 108474054, 306496189, 40518813, 94256451, 158704735, 144801616, 119381492, 127146122, 188991907, 203286758, 39599405, 199159224, 337972074, 141331512, 154965018, 290490595, 180162484, 107460057, 221341376, 152177953, 198495335, 66820464, 107959954, 127924198, 220341538, 185952660, 158527728, 254859071, 34542726, 226591037, 253407037, 35661892, 38713420, 179246641, 316337362, 202350059, 283031412, 236384907, 336771354, 138690781, 337243518, 134799297, 265113663, 35647665, 77837934, 162592511, 133331867, 334582077, 177603945, 292338635, 164763879, 309492783, 73237651, 141899350, 202329986, 232262993, 148118981, 7158709, 188885721, 158207240, 246316364, 149690867, 212768972, 315946142, 199859816, 239309770, 204313553, 241253426, 305358757, 235345469, 169789259, 255160407, 165464838, 207824218, 276827014, 314574402, 265678538, 293845449, 296765368, 6988910, 296480410, 138693719, 277042632, 68995895, 70489301, 84635421, 226397855, 134808297, 72654091, 110709078, 37594348, 4772741, 150338006, 88534030, 211759221, 130706024, 146150333, 172298161, 277509741, 247795506, 50914593, 147497275, 230229458, 270819535, 149107024, 287510408, 287014142, 178970047, 143695200, 25674293, 58356171, 213709347, 252185642, 219725472, 91130312, 17950120, 174887004, 45888391, 62615520, 152871865, 43750940, 246070852, 40650752, 282966732, 38274766, 340730353, 173380019, 303433698, 129458310, 176367657, 97194158, 146149374, 279842995, 71757187, 253587706, 60657878, 173040940, 268048575, 124040399, 73732903, 127492635, 54272494, 216101315, 304344380, 132099452, 215425474, 273434591, 257258304, 2398937, 289587907, 194679745, 308163831, 171359555, 195388508, 242394048, 26825998, 118936504, 19261340, 26288538, 329792612, 282805023, 72947491, 149120214, 325652404, 170738676, 49327270, 6293671, 91370945, 293871714, 260953841, 320724504, 7469553, 180041579, 327818874, 324223131, 101069891, 315186553, 201500660, 263839429, 199699389, 315485375, 336632803, 57125658, 275502711, 59252093, 164218762, 40305822, 322469122, 251762666, 206612512, 76118356, 227096215, 228376086, 96955443, 196731690, 141852728, 167390103, 182439224, 53618185, 179355013, 267452612, 205568939, 278488672, 27069610, 306050745, 119266877, 287701079, 295879217, 226760889, 287674584, 13007866, 316749460, 43087327, 63563671, 327748105, 99247382, 292386433, 263325512, 111800709, 66471103, 230086233, 184987642, 338978604, 67590928, 289387454, 222294953, 197616453, 147831171, 255520911, 38526869, 124610020, 146671841, 339759167, 73104949, 14904348, 73175630, 313153466, 244672283, 81474849, 99692162, 221636784, 131010427, 273445983, 162930677, 190687336, 321912396, 221486971, 151208266, 171890685, 303914970, 53332758, 72664336, 293795011, 229738046, 305650017, 201461268, 89116318, 282452516, 126626212, 271135876, 25544361, 34442616, 189661310, 300515228, 75039893, 297241178, 32311890, 134865187, 87711823, 291435555, 216777091, 216777091, 239456401, 47004096, 211178831, 99975912, 153479121, 170859277, 16105165, 330113552, 122044016, 296658819, 320595612, 256027557, 54733438, 227356252, 78065108, 142795925, 269905908, 167386915, 54388241, 219227795, 39151525, 328445805, 200413143, 68877863, 74209677, 66832076, 204013316, 325974546, 27693439, 131994113, 173942792, 208586485, 234599635, 277791913, 275233847, 26591040, 280188142, 132543317, 41401927, 43488874, 96978787, 242591114, 303671571, 325500662, 263444716, 313728890, 240761476, 342620797, 36023753, 5560184, 309569203, 69840988, 304219980, 124298717, 185416091, 267781748, 37076713, 207342930, 59495824, 295273641, 98717270, 19189529, 12610522, 97844269, 239363626, 1456280, 85295117, 297397161, 92600329, 234720709, 68037156, 218176195, 50604105, 271807436, 249942057, 255824181, 179211909, 1769823, 335636187, 103960827, 112796660, 61068830, 145704697, 295286722, 341402759, 70141962, 115356513, 203605926, 225628921, 144850683, 17829823, 213708540, 268064520, 257360284, 241661803, 93755650, 182102077, 200582308, 335605693, 295274038, 132122238, 144353803, 164778367, 235003482, 11564770, 120162442, 275001701, 225729872, 160169552, 176031685, 278217090, 180048656, 88228657, 195684185, 240496729, 330322268, 194635879, 17652739, 173064547, 224809473, 205027058, 35766561, 14435735, 49484579, 18205364, 135974521, 184512334, 261458038, 254564425, 186757612, 211680436, 265331755, 113734986, 183927647, 138931320, 331273546, 260025685, 202733843, 14046473, 306799968, 254751446, 311293756, 131951878, 258572267, 82921227, 6988853, 251800508, 269043467, 217833412, 176323797, 182682808, 57603392, 292986380, 13104052, 5694806, 127395340, 257211783, 117402728, 177546571, 169372522, 210680653, 154291420, 284425965, 110907179, 151103834, 184753164, 281881503, 287329192, 175240790, 241474709, 208321446, 82183325, 167033432, 105020156, 63619796, 65456757, 254944716, 49004682, 183211085, 281992983, 316258519, 186321038, 45009804, 107506539, 322902584, 181174250, 234292491, 245344339, 51853280, 99215569, 196514517, 266849468, 227904787, 251705055, 26198852, 212792852, 150164377, 275691463, 13437358, 637558, 333194703, 338207980, 10005334, 122537003, 274958766, 88464921, 84495465, 4159244, 121490510, 319570553, 41781169, 158981477, 145980708, 171070813, 301155667, 256645257, 70235757, 23153609, 299870521, 253254497, 191678223, 13876988, 15709718, 18090931, 217576324, 142859860, 255156656, 92365318, 33172590, 326702781, 273914073, 264139899, 180396442, 283368641, 322323112, 202089022, 278763901, 88240785, 291820130, 168917122, 266687184, 216637775, 62207970, 146869220, 252419016, 114640358, 115388442, 54391738, 298093846, 317363068, 193446939, 253982471, 341656688, 299160823, 53085982, 188611559, 229844545, 108078092, 125083097, 40887843, 164002547, 8734091, 149918280, 279120574, 108758059, 150236302, 58630829, 134727145, 236119388, 155683440, 338679681, 29924124, 117790650, 220127837, 253165539, 79183005, 30160473, 4532462, 84562684, 301016877, 193011905, 103979095, 63899982, 153660569, 256908437, 218198587, 280330569, 272272571, 337773364, 212589566, 102107092, 218038802, 230059821, 305219430, 143068556, 88125611, 108484803, 319276805, 299075628, 131869972, 38563447, 331709727, 150955520, 239295922, 119560698, 162328496, 111118360, 34368383, 129713180, 316630657, 292694593, 171220657, 209244869, 27107753, 250154687, 220980715, 198741307, 148254855, 88603693, 11693759, 167225124, 183918517, 153339179, 26881077, 202053110, 310973295, 268107451, 167527652, 145074607, 149410869, 107348464, 112455699, 311353873, 15150850, 5292126, 111201687, 3537558, 188557553, 4832233, 145394135, 238513217, 126188373, 70056609, 2827852, 64899332, 292034233, 262413482, 68428294, 82634082, 148842215, 254017902, 43800275, 109987239, 183559696, 326202126, 309742736, 122739074, 251044053, 22675072, 24865253, 95882527, 140415156, 17999433, 144838771, 141776899, 206233947, 199781319, 174079972, 241640008, 285548526, 289835531, 87551697, 1285436, 6224637, 316328914, 234028458, 311895666, 145677430, 63519779, 305593590, 163010623, 329354173, 47437224, 54857449, 316535070, 169805322, 72593401, 322276311, 30792903, 308395617, 101038708, 285673499, 159357807, 33406197, 19024353, 237103733, 298635516, 300204773, 110831277, 177607284, 38001408, 24770358, 20967185, 83502795, 143608321, 51584511, 161780941, 182867228, 216656844, 218532000, 49115319, 63728679, 214035231, 46712432, 131745080, 143589927, 114968473, 288523886, 50562126, 177943486, 94119153, 325068227, 120974014, 146850779, 112758078, 72201976, 93520380, 53378738, 182530156, 64575166, 180005551, 312124347, 71392330, 134513138, 249113621, 114801513, 278236153, 68836548, 232715515, 341229925, 47000855, 177627569, 331987412, 332340029, 26858573, 207976910, 182982514, 11279424, 97049143, 68957832, 289159331, 54077017, 255105342, 189698883, 30804998, 25022076, 135748367, 284193893, 150710040, 181140640, 82166464, 162068565, 257334013, 130853256, 171770436, 220963405, 47346304, 74374506, 42291762, 58112652, 274117000, 192460099, 191996915, 148764783, 51892979, 113319582, 34395215, 264855504, 200611100, 290714541, 88959908, 66622421, 418322, 186821782, 149218346, 14728920, 90560976, 133505506, 121264776, 92865063, 121272334, 76788367, 203587333, 102067395, 282440650, 50571658, 324391062, 266839425, 50794676, 331164236, 19455971, 300341458, 93037151, 78578862, 191627469, 146053155, 164197607, 285722112, 161541117, 295600727};

static const int32_t invnthroots[1024] = {1, 0, 295600727, 182035460, 285722112, 179378970, 146053155, 151949108, 78578862, 250539426, 300341458, 324120606, 331164236, 292781901, 266839425, 19185515, 50571658, 61135927, 102067395, 139989244, 76788367, 222304243, 92865063, 222311801, 133505506, 253015601, 14728920, 194358231, 186821782, 343158255, 66622421, 254616669, 290714541, 142965477, 264855504, 309181362, 113319582, 291683598, 148764783, 151579662, 192460099, 69459577, 58112652, 301284815, 74374506, 296230273, 220963405, 171806141, 130853256, 86242564, 162068565, 261410113, 181140640, 192866537, 284193893, 207828210, 25022076, 312771579, 189698883, 88471235, 54077017, 54417246, 68957832, 246527434, 11279424, 160594063, 207976910, 316718004, 332340029, 11589165, 177627569, 296575722, 341229925, 110861062, 68836548, 65340424, 114801513, 94462956, 134513138, 272184247, 312124347, 163571026, 64575166, 161046421, 53378738, 250056197, 72201976, 230818499, 146850779, 222602563, 325068227, 249457424, 177943486, 293014451, 288523886, 228608104, 143589927, 211831497, 46712432, 129541346, 63728679, 294461258, 218532000, 126919733, 182867228, 181795636, 51584511, 199968256, 83502795, 322609392, 24770358, 305575169, 177607284, 232745300, 300204773, 44941061, 237103733, 324552224, 33406197, 184218770, 285673499, 242537869, 308395617, 312783674, 322276311, 270983176, 169805322, 27041507, 54857449, 296139353, 329354173, 180565954, 305593590, 280056798, 145677430, 31680911, 234028458, 27247663, 6224637, 342291141, 87551697, 53741046, 285548526, 101936569, 174079972, 143795258, 206233947, 201799678, 144838771, 325577144, 140415156, 247694050, 24865253, 320901505, 251044053, 220837503, 309742736, 17374451, 183559696, 233589338, 43800275, 89558675, 148842215, 260942495, 68428294, 81163095, 292034233, 278677245, 2827852, 273519968, 126188373, 105063360, 145394135, 338744344, 188557553, 340039019, 111201687, 338284451, 15150850, 32222704, 112455699, 236228113, 149410869, 198501970, 167527652, 75469126, 310973295, 141523467, 26881077, 190237398, 183918517, 176351453, 11693759, 254972884, 148254855, 144835270, 220980715, 93421890, 27107753, 134331708, 171220657, 50881984, 316630657, 213863397, 34368383, 232458217, 162328496, 224015879, 239295922, 192621057, 331709727, 305013130, 131869972, 44500949, 319276805, 235091774, 88125611, 200508021, 305219430, 113516756, 218038802, 241469485, 212589566, 5803213, 272272571, 63246008, 218198587, 86668140, 153660569, 279676595, 103979095, 150564672, 301016877, 259013893, 4532462, 313416104, 79183005, 90411038, 220127837, 225785927, 29924124, 4896896, 155683440, 107457189, 134727145, 284945748, 150236302, 234818518, 279120574, 193658297, 8734091, 179574030, 40887843, 218493480, 108078092, 113732032, 188611559, 290490595, 299160823, 1919889, 253982471, 150129638, 317363068, 45482731, 54391738, 228188135, 114640358, 91157561, 146869220, 281368607, 216637775, 76889393, 168917122, 51756447, 88240785, 64812676, 202089022, 21253465, 283368641, 163180135, 264139899, 69662504, 326702781, 310403987, 92365318, 88419921, 142859860, 126000253, 18090931, 327866859, 13876988, 151898354, 253254497, 43706056, 23153609, 273340820, 256645257, 42420910, 171070813, 197595869, 158981477, 301795408, 319570553, 222086067, 4159244, 259081112, 88464921, 68617811, 122537003, 333571243, 338207980, 10381874, 637558, 330139219, 275691463, 193412200, 212792852, 317377725, 251705055, 115671790, 266849468, 147062060, 99215569, 291723297, 245344339, 109284086, 181174250, 20673993, 107506539, 298566773, 186321038, 27318058, 281992983, 160365492, 49004682, 88631861, 65456757, 279956781, 105020156, 176543145, 82183325, 135255131, 241474709, 168335787, 287329192, 61695074, 184753164, 192472743, 110907179, 59150612, 154291420, 132895924, 169372522, 166030006, 117402728, 86364794, 127395340, 337881771, 13104052, 50590197, 57603392, 160893769, 176323797, 125743165, 269043467, 91776069, 6988853, 260655350, 258572267, 211624699, 311293756, 88825131, 306799968, 329530104, 202733843, 83550892, 331273546, 204645257, 183927647, 229841591, 265331755, 131896141, 186757612, 89012152, 261458038, 159064243, 135974521, 325371213, 49484579, 329140842, 35766561, 138549519, 224809473, 170512030, 17652739, 148940698, 330322268, 103079848, 195684185, 255347920, 180048656, 65359487, 176031685, 183407025, 225729872, 68574876, 120162442, 332011807, 235003482, 178798210, 144353803, 211454339, 295274038, 7970884, 200582308, 161474500, 93755650, 101914774, 257360284, 75512057, 213708540, 325746754, 144850683, 117947656, 203605926, 228220064, 70141962, 2173818, 295286722, 197871880, 61068830, 230779917, 103960827, 7940390, 1769823, 164364668, 255824181, 93634520, 271807436, 292972472, 218176195, 275539421, 234720709, 250976248, 297397161, 258281460, 1456280, 104212951, 97844269, 330966055, 19189529, 244859307, 295273641, 284080753, 207342930, 306499864, 267781748, 158160486, 124298717, 39356597, 69840988, 34007374, 5560184, 307552824, 342620797, 102815101, 313728890, 80131861, 325500662, 39905006, 242591114, 246597790, 43488874, 302174650, 132543317, 63388435, 26591040, 68342730, 277791913, 108976942, 208586485, 169633785, 131994113, 315883138, 325974546, 139563261, 66832076, 269366900, 68877863, 143163434, 328445805, 304425052, 219227795, 289188336, 167386915, 73670669, 142795925, 265511469, 227356252, 288843139, 256027557, 22980965, 296658819, 221532561, 330113552, 327471412, 170859277, 190097456, 99975912, 132397746, 47004096, 104120176, 216777091, 126799486, 291435555, 255864754, 134865187, 311264687, 297241178, 268536684, 300515228, 153915267, 34442616, 318032216, 271135876, 216950365, 282452516, 254460259, 201461268, 37926560, 229738046, 49781566, 72664336, 290243819, 303914970, 171685892, 151208266, 122089606, 321912396, 152889241, 162930677, 70130594, 131010427, 121939793, 99692162, 262101728, 244672283, 30423111, 73175630, 328672229, 73104949, 3817410, 146671841, 218966557, 38526869, 88055666, 147831171, 145960124, 222294953, 54189123, 67590928, 4597973, 184987642, 113490344, 66471103, 231775868, 263325512, 51190144, 99247382, 15828472, 63563671, 300489250, 316749460, 330568711, 287674584, 116815688, 295879217, 55875498, 119266877, 37525832, 27069610, 65087905, 205568939, 76123965, 179355013, 289958392, 182439224, 176186474, 141852728, 146844887, 96955443, 115200491, 227096215, 267458221, 206612512, 91813911, 322469122, 303270755, 164218762, 284324484, 275502711, 286450919, 336632803, 28091202, 199699389, 79737148, 201500660, 28390024, 101069891, 19353446, 327818874, 163534998, 7469553, 22852073, 260953841, 49704863, 91370945, 337282906, 49327270, 172837901, 325652404, 194456363, 72947491, 60771554, 329792612, 317288039, 19261340, 224640073, 26825998, 101182529, 195388508, 172217022, 308163831, 148896832, 289587907, 341177640, 257258304, 70141986, 215425474, 211477125, 304344380, 127475262, 54272494, 216083942, 73732903, 219536178, 268048575, 170535637, 60657878, 89988871, 71757187, 63733582, 146149374, 246382419, 176367657, 214118267, 303433698, 170196558, 340730353, 305301811, 282966732, 302925825, 246070852, 299825637, 152871865, 280961057, 45888391, 168689573, 17950120, 252446265, 219725472, 91390935, 213709347, 285220406, 25674293, 199881377, 178970047, 56562435, 287510408, 194469553, 270819535, 113347119, 147497275, 292661984, 247795506, 66066836, 172298161, 197426244, 130706024, 131817356, 88534030, 193238571, 4772741, 305982229, 110709078, 270922486, 134808297, 117178722, 84635421, 273087276, 68995895, 66533945, 138693719, 47096167, 6988910, 46811209, 293845449, 77898039, 314574402, 66749563, 207824218, 178111739, 255160407, 173787318, 235345469, 38217820, 241253426, 139263024, 239309770, 143716761, 315946142, 130807605, 149690867, 97260213, 158207240, 154690856, 7158709, 195457596, 232262993, 141246591, 141899350, 270338926, 309492783, 178812698, 292338635, 165972632, 334582077, 210244710, 162592511, 265738643, 35647665, 78462914, 134799297, 6333059, 138690781, 6805223, 236384907, 60545165, 202350059, 27239215, 179246641, 304863157, 35661892, 90169540, 226591037, 309033851, 254859071, 185048849, 185952660, 123235039, 127924198, 235616623, 66820464, 145081242, 152177953, 122235201, 107460057, 163414093, 290490595, 188611559, 141331512, 5604503, 199159224, 303977172, 203286758, 154584670, 127146122, 224195085, 144801616, 184871842, 94256451, 303057764, 306496189, 235102523, 127741829, 259284818, 52037686, 180635065, 128285449, 204387494, 216322221, 89270496, 159846277, 57633868, 79533366, 311415254, 218690741, 16580936, 218315799, 36685951, 229293967, 139069154, 184840105, 133858494, 231198011, 192939487, 66998696, 247873950, 67657691, 130079597, 176595100, 67430997, 134854900, 241427268, 290857749, 52981975, 230167200, 141474368, 268309021, 90772289, 92509762, 99252841, 313070854, 263448427, 68299003, 76742372, 56557835, 255136374, 93762877, 109793170, 148729034, 61608168, 89256786, 22575024, 236913231, 20204454, 158490431, 103429416, 228968063, 64422955, 18943858, 159346522, 291035159, 342495147, 108886933, 122037970, 59840778, 154704791, 311892219, 336115515, 279105293, 6000276, 287797213, 20888172, 192411606, 17431518, 322312993, 19654506, 341986892, 24086142, 325470081, 76843897, 342801710, 209811274, 156938953, 116891116, 173604417, 318590441, 200468898, 91577194, 214624358, 323024687, 308866754, 232596515, 181490203, 234773563, 340870937, 10526555, 84598502, 136583266, 198322461, 153599252, 274574060, 126821581, 225197457, 202960806, 188613329, 248463082, 274989843, 335471561, 151331725, 55796130, 167187152, 273568873, 221867879, 187604161, 273366424, 27786740, 153465512, 339049400, 331734165, 117359122, 269344708, 146138090, 334574414, 127660000, 173775083, 76570218, 29560260, 291426043, 269316989, 252041855, 82572671, 41288181, 26404288, 9698908, 76050576, 269591847, 21609737, 122732847, 17475810, 283523348, 50088960, 42451231, 222163738, 186347547, 198124793, 284154686, 214506407, 249180058, 257645467, 335023173, 187796153, 171706691, 111564452, 277213025, 29345559, 229608819, 211052972, 16772588, 128706839, 178337622, 213190377, 217130767, 89157113, 1234647, 320794310, 193654125, 163779570, 94005300, 265738498, 124705282, 76745925, 225375689, 291093134, 231029602, 101944834, 273928883, 326679389, 287853095, 2414777, 167957106, 42888295, 122004631, 183276596, 136605582, 269680462, 170878007, 312804364, 145807752, 208781989, 303710857, 58749261, 78357007, 89570472, 210958162, 177303797, 37504712, 117534032, 160751385, 321637827, 251246724, 304854525, 165581902, 247002466, 77680626, 201397764, 177454258, 178225281, 307902104, 69315579, 36641720, 326187203, 23477409, 24059662, 128775820, 271292966, 137469474, 141073720, 292763573, 147041113, 214600928, 215710474, 66294609, 15612353, 16976763, 114094683, 72761504, 175281391, 209198507, 58343183, 328223279, 331758667, 208657514, 289796272, 52529357, 341735197, 130130763, 144613436, 82887967, 233756460, 173671139, 230765503, 247977898};

static const int32_t gj[256] = {1, 88890923, 85045001, 155784997, 27338650, 160916864, 180254835, 310102440, 72009511, 161373771, 190889672, 142668801, 224600123, 311895909, 219297628, 66561745, 311987596, 37795135, 117153417, 261971914, 245468470, 236293451, 222945128, 271831692, 61248299, 318174877, 238779399, 288166042, 280356729, 129754441, 241345284, 176333293, 181397604, 250532407, 43574444, 161778452, 190888335, 173000392, 243280773, 234760582, 86639740, 343228356, 158453678, 127761294, 41878770, 271330443, 303373717, 173106749, 218508975, 25389260, 274716421, 311455255, 213328025, 261827366, 317157312, 321413963, 165631929, 342440146, 328230304, 315578669, 269400584, 109770229, 125471905, 44345055, 248379607, 259244350, 26400835, 279147437, 194920492, 153095710, 190865232, 83207052, 19489956, 98618966, 319133893, 179482504, 222690490, 89873578, 3931971, 202536057, 328223873, 25900407, 318909737, 125753054, 221821179, 40642595, 256906520, 108905618, 74461390, 85688676, 182549481, 343073646, 183403927, 80672799, 98292873, 136451813, 70537120, 27609258, 272499124, 248028252, 130572716, 238552747, 313161526, 337489508, 169856110, 285211679, 75228823, 96882831, 29081840, 205601080, 261914872, 232562070, 57037072, 199733782, 314641395, 33323143, 94520109, 71887187, 179617615, 301362331, 29962501, 82401464, 163096880, 97266640, 11251447, 93268581, 1286289, 223818340, 72210528, 17425846, 342113900, 288300085, 250658058, 289575201, 232501249, 143224961, 223094696, 126602787, 179924173, 92307187, 321139722, 75308790, 181520819, 67911669, 121350736, 273979971, 297887177, 83716059, 178531710, 12200435, 217013157, 205181400, 64080892, 312090881, 83486596, 111314120, 137132995, 62047357, 17679316, 116023627, 306685998, 161429615, 221209188, 136680774, 179466374, 157158321, 215279409, 111014539, 38859848, 116670211, 203983008, 154360503, 244747461, 199042301, 123516309, 285059136, 323926093, 185266231, 279385828, 218405756, 5696908, 52723975, 67230089, 139640576, 106504061, 295717728, 138726077, 8006261, 95624257, 316713932, 220275047, 98837722, 338357212, 299229864, 108398169, 289693716, 51660043, 89055604, 307495902, 239586427, 337001594, 89672506, 68951209, 183908852, 203832179, 214099607, 80217350, 162254241, 268352613, 291931314, 225532465, 198076789, 299603333, 25767431, 24511420, 4732033, 162193745, 40081432, 125214816, 161861063, 296748067, 28154930, 106246674, 35558534, 268350937, 79402207, 274390244, 220620337, 175710674, 36716312, 264795720, 104403594, 69631836, 315983912, 63917731, 178248179, 83293734, 214161640, 190359536, 140461439, 300690235, 296422692, 107488819, 246557879, 31717623, 152685717, 76500391, 97028904, 146879235, 285737564, 136357212, 280629049, 218173266, 229691765, 85294681, 81718591};

static const int32_t invgj[256] = {1, 261857986, 258281896, 113884812, 125403311, 62947528, 207219365, 57839013, 196697342, 246547673, 267076186, 190890860, 311858954, 97018698, 236087758, 47153885, 42886342, 203115138, 153217041, 129414937, 260282843, 165328398, 279658846, 27592665, 273944741, 239172983, 78780857, 306860265, 167865903, 122956240, 69186333, 264174370, 75225640, 308018043, 237329903, 315421647, 46828510, 181715514, 218361761, 303495145, 181382832, 338844544, 319065157, 317809146, 43973244, 145499788, 118044112, 51645263, 75223964, 181322336, 263359227, 129476970, 139744398, 159667725, 274625368, 253904071, 6574983, 103990150, 36080675, 254520973, 291916534, 53882861, 235178408, 44346713, 5219365, 244738855, 123301530, 26862645, 247952320, 335570316, 204850500, 47858849, 237072516, 203936001, 276346488, 290852602, 337879669, 125170821, 64190749, 158310346, 19650484, 58517441, 220060268, 144534276, 98829116, 189216074, 139593569, 226906366, 304716729, 232562038, 128297168, 186418256, 164110203, 206895803, 122367389, 182146962, 36890579, 227552950, 325897261, 281529220, 206443582, 232262457, 260089981, 31485696, 279495685, 138395177, 126563420, 331376142, 165044867, 259860518, 45689400, 69596606, 222225841, 275664908, 162055758, 268267787, 22436855, 251269390, 163652404, 216973790, 120481881, 200351616, 111075328, 54001376, 92918519, 55276492, 1462677, 326150731, 271366049, 119758237, 342290288, 250307996, 332325130, 246309937, 180479697, 261175113, 313614076, 42214246, 163958962, 271689390, 249056468, 310253434, 28935182, 143842795, 286539505, 111014507, 81661705, 137975497, 314494737, 246693746, 268347754, 58364898, 173720467, 6087069, 30415051, 105023830, 213003861, 95548325, 71077453, 315967319, 273039457, 207124764, 245283704, 262903778, 160172650, 502931, 161027096, 257887901, 269115187, 234670959, 86670057, 302933982, 121755398, 217823523, 24666840, 317676170, 15352704, 141040520, 339644606, 253702999, 120886087, 164094073, 24442684, 244957611, 324086621, 260369525, 152711345, 190480867, 148656085, 64429140, 317175742, 84332227, 95196970, 299231522, 218104672, 233806348, 74175993, 27997908, 15346273, 1136431, 177944648, 22162614, 26419265, 81749211, 130248552, 32121322, 68860156, 318187317, 125067602, 170469828, 40202860, 72246134, 301697807, 215815283, 185122899, 348221, 256936837, 108815995, 100295804, 170576185, 152688242, 181798125, 300002133, 93044170, 162178973, 167243284, 102231293, 213822136, 63219848, 55410535, 104797178, 25401700, 282328278, 71744885, 120631449, 107283126, 98108107, 81604663, 226423160, 305781442, 31588981, 277014832, 124278949, 31680668, 118976454, 200907776, 152686905, 182202806, 271567066, 33474137, 163321742, 182659713, 316237927, 187791580, 258531576, 254685654};

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
  int64_t x_64[PARAM_N];
  int i, index, j, m, window;
  int64_t a, sub_re, sub_img;

  for(i = 0; i < PARAM_N; ++i)
    x_64[i] = (int64_t) x[i];

  window = 1;
  m = 512;
  for(; m >= 16; m >>= 1) 
  {
    index = 0;

    for(j = 0; j < m; j += 2) 
    {
      a = gj[index];

      for(i = j; i < PARAM_N; i += (m << 1)) 
      {
        sub_re = x_64[i] - x_64[i+m];
        sub_img = x_64[i+1] - x_64[i+m+1];
        
        x_64[i] = x_64[i] + x_64[i+m];
        x_64[i+1] = x_64[i+1] + x_64[i+m+1];
        
        x_64[i+m] = reduce((int64_t)a * sub_re);
        x_64[i+m+1] = reduce((int64_t)a * sub_img);        
      }
      index += window;
    }
    window <<= 1;
  }
    
  index = 0;    
  for(j = 0; j < m; j += 2) 
  {
    a = gj[index];

    for(i = j; i < PARAM_N; i += (m << 1)) 
    {
      sub_re = x_64[i] - x_64[i+m];
      sub_img = x_64[i+1] - x_64[i+m+1];
      
      x_64[i] = barr_reduce64(x_64[i] + x_64[i+m]);
      x_64[i+1] = barr_reduce64(x_64[i+1] + x_64[i+m+1]);
      
      x_64[i+m] = reduce((int64_t)a * sub_re);
      x_64[i+m+1] = reduce((int64_t)a * sub_img);        
    }
    index += window;
  }
  window <<= 1;

  m >>= 1;
  
  index = 0;    
  for(j = 0; j < m; j += 2) 
  {
    a = gj[index];

    for(i = j; i < PARAM_N; i += (m << 1)) 
    {
      sub_re = x_64[i] - x_64[i+m];
      sub_img = x_64[i+1] - x_64[i+m+1];
      
      x_64[i] = x_64[i] + x_64[i+m];
      x_64[i+1] = x_64[i+1] + x_64[i+m+1];
      
      x_64[i+m] = reduce((int64_t)a * sub_re);
      x_64[i+m+1] = reduce((int64_t)a * sub_img);        
    }
    index += window;
  }
  window <<= 1;

  m >>= 1;
  
  index = 0;    
  for(j = 0; j < m; j += 2) 
  {
    a = gj[index];

    for(i = j; i < PARAM_N; i += (m << 1)) 
    {
      sub_re = x_64[i] - x_64[i+m];
      sub_img = x_64[i+1] - x_64[i+m+1];
      
      x_64[i] = x_64[i] + x_64[i+m];
      x_64[i+1] = x_64[i+1] + x_64[i+m+1];
      
      x_64[i+m] = reduce((int64_t)a * sub_re);
      x_64[i+m+1] = reduce((int64_t)a * sub_img);        
    }
    index += window;
  }

  for(i = 0; i < PARAM_N; ++i)
    x[i] = (int32_t) barr_reduce64(x_64[i]);

}


void idgt(poly x)
{

  int64_t x_64[PARAM_N];
  int i, index, j, m, window;
  int64_t a, mul_re, mul_img;

  for(i = 0; i < PARAM_N; ++i)
    x_64[i] = (int64_t) x[i];

  window = 256;
  for(m = 2; m <= 512; m <<= 1) 
  {
    index = 0;
    for(j = 0; j < m; j += 2) 
    {
      a = invgj[index];
      for(i = j; i < PARAM_N; i += (m << 1)) 
      {
        mul_re = reduce(x_64[i+m] * a);
        mul_img = reduce(x_64[i+m+1] * a);

        x_64[i+m] = x_64[i] - mul_re;
        x_64[i+m+1] = x_64[i+1] - mul_img;

        x_64[i] = x_64[i] + mul_re;
        x_64[i+1] = x_64[i+1] + mul_img;        
      }
      index += window;
    }
    window >>= 1;
  }  

  for(i = 0; i < PARAM_N; ++i)
    x[i] = (int32_t) barr_reduce64(x_64[i]);

}


static void poly_pointwise(poly result, const poly x, const poly y)
{ // Pointwise polynomial multiplication result = x.y

  for (int i=0; i<PARAM_N; i++)
    result[i] = reduce((int64_t)x[i]*y[i]);
}


void poly_dgt(poly x_dgt, const poly x)
{

  int i, j;

  for(i = 0, j = 0; i < PARAM_N && j < 512; i+=2, j++) {             
      x_dgt[i] = reduce(
        (int64_t)x[j] * nthroots[i] - 
        (int64_t)x[512+j] * nthroots[i+1]
      );
      
      x_dgt[i+1] = reduce(
        (int64_t)x[j] * nthroots[i+1] + 
        (int64_t)x[512+j] * nthroots[i]
      );
  } 

  dgt(x_dgt);
}


void poly_mul(poly _output, const poly _poly_a, const poly _poly_b)
{ /* It is assumed that both signals are already in the DGT domain. 
     The DGT counterpart of poly_b was computed in sign.c. */
  
  poly _mul;
  int i, j;

  /* Calculating the point-wise multiplication of input signals */
  for(i = 0, j = 0; i < PARAM_N && j < 512; i+=2, j++) {             
    _mul[i] = reduce(
      (int64_t)_poly_a[j] * _poly_b[i] -
      (int64_t)_poly_a[j+512] * _poly_b[i+1]
    );

    _mul[i+1] = reduce(
      (int64_t)_poly_a[j] * _poly_b[i+1] + 
      (int64_t)_poly_a[j+512] * _poly_b[i]
    );
  }

  /* Recovering the multiplication result in Z[x]/<x^n+1> */
  idgt(_mul);

  /* Removing the twisting factors and writing the result from the Gaussian integer to the polynomial form */
  for(i = 0, j = 0; i < PARAM_N && j < 512; i+=2, j++) {
      _output[j] = reduce(
               (int64_t)_mul[i] * invnthroots[i] -
               (int64_t)_mul[i+1] * invnthroots[i+1]);

      _output[j+512] = reduce(
               (int64_t)_mul[i] * invnthroots[i+1] + 
               (int64_t)_mul[i+1] * invnthroots[i]);
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
      result[i] = (int32_t) barr_reduce(x[i] - y[i]);
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