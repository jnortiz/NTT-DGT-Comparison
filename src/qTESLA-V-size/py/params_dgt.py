from gauss import GaussianInteger

invofkmodp = 33302449
gj_powers = [[1, 1, 1, 1, 1, 1, 1], [4406558, 19477423, 20484768, 10344459, 29493997, 27506971, 33564672], [19477423, 20484768, 10344459, 29493997, 27506971, 33564672, 1], [988369, 7037169, 26214285, 22303942, 17398315, 6057702, 33564672], [20484768, 10344459, 29493997, 27506971, 33564672, 1, 1], [30717302, 8413164, 11809804, 16978151, 4070676, 27506971, 33564672], [7037169, 26214285, 22303942, 17398315, 6057702, 33564672, 1], [30392408, 5665384, 21664476, 30033475, 16166358, 6057702, 33564672], [10344459, 29493997, 27506971, 33564672, 1, 1, 1], [14583628, 24389287, 23442917, 23220214, 29493997, 27506971, 33564672], [8413164, 11809804, 16978151, 4070676, 27506971, 33564672, 1], [7554841, 3732990, 1227325, 11260731, 17398315, 6057702, 33564672], [26214285, 22303942, 17398315, 6057702, 33564672, 1, 1], [31965169, 12975937, 7756560, 16586522, 4070676, 27506971, 33564672], [5665384, 21664476, 30033475, 16166358, 6057702, 33564672, 1], [10010313, 8491986, 25589677, 3531198, 16166358, 6057702, 33564672], [29493997, 27506971, 33564672, 1, 1, 1, 1], [29780798, 27454015, 13079905, 10344459, 29493997, 27506971, 33564672], [24389287, 23442917, 23220214, 29493997, 27506971, 33564672, 1], [255720, 8735396, 7350388, 22303942, 17398315, 6057702, 33564672], [11809804, 16978151, 4070676, 27506971, 33564672, 1, 1], [4089071, 24835380, 21754869, 16978151, 4070676, 27506971, 33564672], [3732990, 1227325, 11260731, 17398315, 6057702, 33564672, 1], [27051869, 32426145, 11900197, 30033475, 16166358, 6057702, 33564672], [22303942, 17398315, 6057702, 33564672, 1, 1, 1], [8478458, 31786565, 10121756, 23220214, 29493997, 27506971, 33564672], [12975937, 7756560, 16586522, 4070676, 27506971, 33564672, 1], [19611677, 18327945, 32337348, 11260731, 17398315, 6057702, 33564672], [21664476, 30033475, 16166358, 6057702, 33564672, 1, 1], [19452799, 9759120, 25808113, 16586522, 4070676, 27506971, 33564672], [8491986, 25589677, 3531198, 16166358, 6057702, 33564672, 1], [14033313, 2121761, 7974996, 3531198, 16166358, 6057702, 33564672], [27506971, 33564672, 1, 1, 1, 1, 1], [15781, 14087250, 20484768, 10344459, 29493997, 27506971, 33564672], [27454015, 13079905, 10344459, 29493997, 27506971, 33564672, 1], [21501702, 26527504, 26214285, 22303942, 17398315, 6057702, 33564672], [23442917, 23220214, 29493997, 27506971, 33564672, 1, 1], [8758145, 25151509, 11809804, 16978151, 4070676, 27506971, 33564672], [8735396, 7350388, 22303942, 17398315, 6057702, 33564672, 1], [21625705, 27899289, 21664476, 30033475, 16166358, 6057702, 33564672], [16978151, 4070676, 27506971, 33564672, 1, 1, 1], [20902680, 9175386, 23442917, 23220214, 29493997, 27506971, 33564672], [24835380, 21754869, 16978151, 4070676, 27506971, 33564672, 1], [19859369, 29831683, 1227325, 11260731, 17398315, 6057702, 33564672], [1227325, 11260731, 17398315, 6057702, 33564672, 1, 1], [3036860, 20588736, 7756560, 16586522, 4070676, 27506971, 33564672], [32426145, 11900197, 30033475, 16166358, 6057702, 33564672, 1], [22700705, 25072687, 25589677, 3531198, 16166358, 6057702, 33564672], [17398315, 6057702, 33564672, 1, 1, 1, 1], [3446166, 6110658, 13079905, 10344459, 29493997, 27506971, 33564672], [31786565, 10121756, 23220214, 29493997, 27506971, 33564672, 1], [1232856, 24829277, 7350388, 22303942, 17398315, 6057702, 33564672], [7756560, 16586522, 4070676, 27506971, 33564672, 1, 1], [19452428, 8729293, 21754869, 16978151, 4070676, 27506971, 33564672], [18327945, 32337348, 11260731, 17398315, 6057702, 33564672, 1], [4314075, 1138528, 11900197, 30033475, 16166358, 6057702, 33564672], [30033475, 16166358, 6057702, 33564672, 1, 1, 1], [19347624, 1778108, 10121756, 23220214, 29493997, 27506971, 33564672], [9759120, 25808113, 16586522, 4070676, 27506971, 33564672, 1], [28756497, 15236728, 32337348, 11260731, 17398315, 6057702, 33564672], [25589677, 3531198, 16166358, 6057702, 33564672, 1, 1], [30901251, 23805553, 25808113, 16586522, 4070676, 27506971, 33564672], [2121761, 7974996, 3531198, 16166358, 6057702, 33564672, 1], [21856450, 31442912, 7974996, 3531198, 16166358, 6057702, 33564672], [33564672, 1, 1, 1, 1, 1, 1], [29158115, 19477423, 20484768, 10344459, 29493997, 27506971, 33564672], [14087250, 20484768, 10344459, 29493997, 27506971, 33564672, 1], [32576304, 7037169, 26214285, 22303942, 17398315, 6057702, 33564672], [13079905, 10344459, 29493997, 27506971, 33564672, 1, 1], [2847371, 8413164, 11809804, 16978151, 4070676, 27506971, 33564672], [26527504, 26214285, 22303942, 17398315, 6057702, 33564672, 1], [3172265, 5665384, 21664476, 30033475, 16166358, 6057702, 33564672], [23220214, 29493997, 27506971, 33564672, 1, 1, 1], [18981045, 24389287, 23442917, 23220214, 29493997, 27506971, 33564672], [25151509, 11809804, 16978151, 4070676, 27506971, 33564672, 1], [26009832, 3732990, 1227325, 11260731, 17398315, 6057702, 33564672], [7350388, 22303942, 17398315, 6057702, 33564672, 1, 1], [1599504, 12975937, 7756560, 16586522, 4070676, 27506971, 33564672], [27899289, 21664476, 30033475, 16166358, 6057702, 33564672, 1], [23554360, 8491986, 25589677, 3531198, 16166358, 6057702, 33564672], [4070676, 27506971, 33564672, 1, 1, 1, 1], [3783875, 27454015, 13079905, 10344459, 29493997, 27506971, 33564672], [9175386, 23442917, 23220214, 29493997, 27506971, 33564672, 1], [33308953, 8735396, 7350388, 22303942, 17398315, 6057702, 33564672], [21754869, 16978151, 4070676, 27506971, 33564672, 1, 1], [29475602, 24835380, 21754869, 16978151, 4070676, 27506971, 33564672], [29831683, 1227325, 11260731, 17398315, 6057702, 33564672, 1], [6512804, 32426145, 11900197, 30033475, 16166358, 6057702, 33564672], [11260731, 17398315, 6057702, 33564672, 1, 1, 1], [25086215, 31786565, 10121756, 23220214, 29493997, 27506971, 33564672], [20588736, 7756560, 16586522, 4070676, 27506971, 33564672, 1], [13952996, 18327945, 32337348, 11260731, 17398315, 6057702, 33564672], [11900197, 30033475, 16166358, 6057702, 33564672, 1, 1], [14111874, 9759120, 25808113, 16586522, 4070676, 27506971, 33564672], [25072687, 25589677, 3531198, 16166358, 6057702, 33564672, 1], [19531360, 2121761, 7974996, 3531198, 16166358, 6057702, 33564672], [6057702, 33564672, 1, 1, 1, 1, 1], [33548892, 14087250, 20484768, 10344459, 29493997, 27506971, 33564672], [6110658, 13079905, 10344459, 29493997, 27506971, 33564672, 1], [12062971, 26527504, 26214285, 22303942, 17398315, 6057702, 33564672], [10121756, 23220214, 29493997, 27506971, 33564672, 1, 1], [24806528, 25151509, 11809804, 16978151, 4070676, 27506971, 33564672], [24829277, 7350388, 22303942, 17398315, 6057702, 33564672, 1], [11938968, 27899289, 21664476, 30033475, 16166358, 6057702, 33564672], [16586522, 4070676, 27506971, 33564672, 1, 1, 1], [12661993, 9175386, 23442917, 23220214, 29493997, 27506971, 33564672], [8729293, 21754869, 16978151, 4070676, 27506971, 33564672, 1], [13705304, 29831683, 1227325, 11260731, 17398315, 6057702, 33564672], [32337348, 11260731, 17398315, 6057702, 33564672, 1, 1], [30527813, 20588736, 7756560, 16586522, 4070676, 27506971, 33564672], [1138528, 11900197, 30033475, 16166358, 6057702, 33564672, 1], [10863968, 25072687, 25589677, 3531198, 16166358, 6057702, 33564672], [16166358, 6057702, 33564672, 1, 1, 1, 1], [30118507, 6110658, 13079905, 10344459, 29493997, 27506971, 33564672], [1778108, 10121756, 23220214, 29493997, 27506971, 33564672, 1], [32331817, 24829277, 7350388, 22303942, 17398315, 6057702, 33564672], [25808113, 16586522, 4070676, 27506971, 33564672, 1, 1], [14112245, 8729293, 21754869, 16978151, 4070676, 27506971, 33564672], [15236728, 32337348, 11260731, 17398315, 6057702, 33564672, 1], [29250598, 1138528, 11900197, 30033475, 16166358, 6057702, 33564672], [3531198, 16166358, 6057702, 33564672, 1, 1, 1], [14217049, 1778108, 10121756, 23220214, 29493997, 27506971, 33564672], [23805553, 25808113, 16586522, 4070676, 27506971, 33564672, 1], [4808176, 15236728, 32337348, 11260731, 17398315, 6057702, 33564672], [7974996, 3531198, 16166358, 6057702, 33564672, 1, 1], [2663422, 23805553, 25808113, 16586522, 4070676, 27506971, 33564672], [31442912, 7974996, 3531198, 16166358, 6057702, 33564672, 1], [11708223, 31442912, 7974996, 3531198, 16166358, 6057702, 33564672]]
invgj_powers = [[1, 1, 1, 1, 1, 1, 1], [33564672, 6057702, 16166358, 3531198, 7974996, 31442912, 11708223], [1, 33564672, 6057702, 16166358, 3531198, 7974996, 31442912], [33564672, 27506971, 4070676, 16586522, 25808113, 23805553, 2663422], [1, 1, 33564672, 6057702, 16166358, 3531198, 7974996], [33564672, 6057702, 17398315, 11260731, 32337348, 15236728, 4808176], [1, 33564672, 27506971, 4070676, 16586522, 25808113, 23805553], [33564672, 27506971, 29493997, 23220214, 10121756, 1778108, 14217049], [1, 1, 1, 33564672, 6057702, 16166358, 3531198], [33564672, 6057702, 16166358, 30033475, 11900197, 1138528, 29250598], [1, 33564672, 6057702, 17398315, 11260731, 32337348, 15236728], [33564672, 27506971, 4070676, 16978151, 21754869, 8729293, 14112245], [1, 1, 33564672, 27506971, 4070676, 16586522, 25808113], [33564672, 6057702, 17398315, 22303942, 7350388, 24829277, 32331817], [1, 33564672, 27506971, 29493997, 23220214, 10121756, 1778108], [33564672, 27506971, 29493997, 10344459, 13079905, 6110658, 30118507], [1, 1, 1, 1, 33564672, 6057702, 16166358], [33564672, 6057702, 16166358, 3531198, 25589677, 25072687, 10863968], [1, 33564672, 6057702, 16166358, 30033475, 11900197, 1138528], [33564672, 27506971, 4070676, 16586522, 7756560, 20588736, 30527813], [1, 1, 33564672, 6057702, 17398315, 11260731, 32337348], [33564672, 6057702, 17398315, 11260731, 1227325, 29831683, 13705304], [1, 33564672, 27506971, 4070676, 16978151, 21754869, 8729293], [33564672, 27506971, 29493997, 23220214, 23442917, 9175386, 12661993], [1, 1, 1, 33564672, 27506971, 4070676, 16586522], [33564672, 6057702, 16166358, 30033475, 21664476, 27899289, 11938968], [1, 33564672, 6057702, 17398315, 22303942, 7350388, 24829277], [33564672, 27506971, 4070676, 16978151, 11809804, 25151509, 24806528], [1, 1, 33564672, 27506971, 29493997, 23220214, 10121756], [33564672, 6057702, 17398315, 22303942, 26214285, 26527504, 12062971], [1, 33564672, 27506971, 29493997, 10344459, 13079905, 6110658], [33564672, 27506971, 29493997, 10344459, 20484768, 14087250, 33548892], [1, 1, 1, 1, 1, 33564672, 6057702], [33564672, 6057702, 16166358, 3531198, 7974996, 2121761, 19531360], [1, 33564672, 6057702, 16166358, 3531198, 25589677, 25072687], [33564672, 27506971, 4070676, 16586522, 25808113, 9759120, 14111874], [1, 1, 33564672, 6057702, 16166358, 30033475, 11900197], [33564672, 6057702, 17398315, 11260731, 32337348, 18327945, 13952996], [1, 33564672, 27506971, 4070676, 16586522, 7756560, 20588736], [33564672, 27506971, 29493997, 23220214, 10121756, 31786565, 25086215], [1, 1, 1, 33564672, 6057702, 17398315, 11260731], [33564672, 6057702, 16166358, 30033475, 11900197, 32426145, 6512804], [1, 33564672, 6057702, 17398315, 11260731, 1227325, 29831683], [33564672, 27506971, 4070676, 16978151, 21754869, 24835380, 29475602], [1, 1, 33564672, 27506971, 4070676, 16978151, 21754869], [33564672, 6057702, 17398315, 22303942, 7350388, 8735396, 33308953], [1, 33564672, 27506971, 29493997, 23220214, 23442917, 9175386], [33564672, 27506971, 29493997, 10344459, 13079905, 27454015, 3783875], [1, 1, 1, 1, 33564672, 27506971, 4070676], [33564672, 6057702, 16166358, 3531198, 25589677, 8491986, 23554360], [1, 33564672, 6057702, 16166358, 30033475, 21664476, 27899289], [33564672, 27506971, 4070676, 16586522, 7756560, 12975937, 1599504], [1, 1, 33564672, 6057702, 17398315, 22303942, 7350388], [33564672, 6057702, 17398315, 11260731, 1227325, 3732990, 26009832], [1, 33564672, 27506971, 4070676, 16978151, 11809804, 25151509], [33564672, 27506971, 29493997, 23220214, 23442917, 24389287, 18981045], [1, 1, 1, 33564672, 27506971, 29493997, 23220214], [33564672, 6057702, 16166358, 30033475, 21664476, 5665384, 3172265], [1, 33564672, 6057702, 17398315, 22303942, 26214285, 26527504], [33564672, 27506971, 4070676, 16978151, 11809804, 8413164, 2847371], [1, 1, 33564672, 27506971, 29493997, 10344459, 13079905], [33564672, 6057702, 17398315, 22303942, 26214285, 7037169, 32576304], [1, 33564672, 27506971, 29493997, 10344459, 20484768, 14087250], [33564672, 27506971, 29493997, 10344459, 20484768, 19477423, 29158115], [1, 1, 1, 1, 1, 1, 33564672], [33564672, 6057702, 16166358, 3531198, 7974996, 31442912, 21856450], [1, 33564672, 6057702, 16166358, 3531198, 7974996, 2121761], [33564672, 27506971, 4070676, 16586522, 25808113, 23805553, 30901251], [1, 1, 33564672, 6057702, 16166358, 3531198, 25589677], [33564672, 6057702, 17398315, 11260731, 32337348, 15236728, 28756497], [1, 33564672, 27506971, 4070676, 16586522, 25808113, 9759120], [33564672, 27506971, 29493997, 23220214, 10121756, 1778108, 19347624], [1, 1, 1, 33564672, 6057702, 16166358, 30033475], [33564672, 6057702, 16166358, 30033475, 11900197, 1138528, 4314075], [1, 33564672, 6057702, 17398315, 11260731, 32337348, 18327945], [33564672, 27506971, 4070676, 16978151, 21754869, 8729293, 19452428], [1, 1, 33564672, 27506971, 4070676, 16586522, 7756560], [33564672, 6057702, 17398315, 22303942, 7350388, 24829277, 1232856], [1, 33564672, 27506971, 29493997, 23220214, 10121756, 31786565], [33564672, 27506971, 29493997, 10344459, 13079905, 6110658, 3446166], [1, 1, 1, 1, 33564672, 6057702, 17398315], [33564672, 6057702, 16166358, 3531198, 25589677, 25072687, 22700705], [1, 33564672, 6057702, 16166358, 30033475, 11900197, 32426145], [33564672, 27506971, 4070676, 16586522, 7756560, 20588736, 3036860], [1, 1, 33564672, 6057702, 17398315, 11260731, 1227325], [33564672, 6057702, 17398315, 11260731, 1227325, 29831683, 19859369], [1, 33564672, 27506971, 4070676, 16978151, 21754869, 24835380], [33564672, 27506971, 29493997, 23220214, 23442917, 9175386, 20902680], [1, 1, 1, 33564672, 27506971, 4070676, 16978151], [33564672, 6057702, 16166358, 30033475, 21664476, 27899289, 21625705], [1, 33564672, 6057702, 17398315, 22303942, 7350388, 8735396], [33564672, 27506971, 4070676, 16978151, 11809804, 25151509, 8758145], [1, 1, 33564672, 27506971, 29493997, 23220214, 23442917], [33564672, 6057702, 17398315, 22303942, 26214285, 26527504, 21501702], [1, 33564672, 27506971, 29493997, 10344459, 13079905, 27454015], [33564672, 27506971, 29493997, 10344459, 20484768, 14087250, 15781], [1, 1, 1, 1, 1, 33564672, 27506971], [33564672, 6057702, 16166358, 3531198, 7974996, 2121761, 14033313], [1, 33564672, 6057702, 16166358, 3531198, 25589677, 8491986], [33564672, 27506971, 4070676, 16586522, 25808113, 9759120, 19452799], [1, 1, 33564672, 6057702, 16166358, 30033475, 21664476], [33564672, 6057702, 17398315, 11260731, 32337348, 18327945, 19611677], [1, 33564672, 27506971, 4070676, 16586522, 7756560, 12975937], [33564672, 27506971, 29493997, 23220214, 10121756, 31786565, 8478458], [1, 1, 1, 33564672, 6057702, 17398315, 22303942], [33564672, 6057702, 16166358, 30033475, 11900197, 32426145, 27051869], [1, 33564672, 6057702, 17398315, 11260731, 1227325, 3732990], [33564672, 27506971, 4070676, 16978151, 21754869, 24835380, 4089071], [1, 1, 33564672, 27506971, 4070676, 16978151, 11809804], [33564672, 6057702, 17398315, 22303942, 7350388, 8735396, 255720], [1, 33564672, 27506971, 29493997, 23220214, 23442917, 24389287], [33564672, 27506971, 29493997, 10344459, 13079905, 27454015, 29780798], [1, 1, 1, 1, 33564672, 27506971, 29493997], [33564672, 6057702, 16166358, 3531198, 25589677, 8491986, 10010313], [1, 33564672, 6057702, 16166358, 30033475, 21664476, 5665384], [33564672, 27506971, 4070676, 16586522, 7756560, 12975937, 31965169], [1, 1, 33564672, 6057702, 17398315, 22303942, 26214285], [33564672, 6057702, 17398315, 11260731, 1227325, 3732990, 7554841], [1, 33564672, 27506971, 4070676, 16978151, 11809804, 8413164], [33564672, 27506971, 29493997, 23220214, 23442917, 24389287, 14583628], [1, 1, 1, 33564672, 27506971, 29493997, 10344459], [33564672, 6057702, 16166358, 30033475, 21664476, 5665384, 30392408], [1, 33564672, 6057702, 17398315, 22303942, 26214285, 7037169], [33564672, 27506971, 4070676, 16978151, 11809804, 8413164, 30717302], [1, 1, 33564672, 27506971, 29493997, 10344459, 20484768], [33564672, 6057702, 17398315, 22303942, 26214285, 7037169, 988369], [1, 33564672, 27506971, 29493997, 10344459, 20484768, 19477423], [33564672, 27506971, 29493997, 10344459, 20484768, 19477423, 4406558]]
nthroots = [GaussianInteger(1,0), GaussianInteger(2851919,23993255), GaussianInteger(30540791,29178387), GaussianInteger(23895260,7368331), GaussianInteger(6950414,20987741), GaussianInteger(19549050,25041160), GaussianInteger(29741358,21035353), GaussianInteger(30423927,33007623), GaussianInteger(6048375,20601016), GaussianInteger(13696011,3306177), GaussianInteger(27511462,2413724), GaussianInteger(12764202,18281015), GaussianInteger(15995956,4451054), GaussianInteger(11407476,19649036), GaussianInteger(31517277,21523023), GaussianInteger(21266122,2561481), GaussianInteger(33455455,28637337), GaussianInteger(28131931,23618973), GaussianInteger(10333309,15018214), GaussianInteger(2183338,23403745), GaussianInteger(9921794,10801637), GaussianInteger(28816607,4698604), GaussianInteger(1571717,9871699), GaussianInteger(7219016,2010662), GaussianInteger(26102997,5448824), GaussianInteger(23201442,9672268), GaussianInteger(19857380,26717892), GaussianInteger(24027111,18048009), GaussianInteger(4887726,3474842), GaussianInteger(21433569,4395153), GaussianInteger(28502692,10247717), GaussianInteger(19180103,24457267), GaussianInteger(10058860,22762078), GaussianInteger(23567528,6833904), GaussianInteger(15674778,23028000), GaussianInteger(19003874,15819278), GaussianInteger(24796356,5469990), GaussianInteger(31865002,32746868), GaussianInteger(15690159,30014775), GaussianInteger(11100323,27174900), GaussianInteger(5177130,4563035), GaussianInteger(18585800,4571066), GaussianInteger(30943380,30088708), GaussianInteger(4120675,19915016), GaussianInteger(23638006,24180444), GaussianInteger(1692056,20468383), GaussianInteger(4965393,17313503), GaussianInteger(5809600,33026417), GaussianInteger(14710518,31647269), GaussianInteger(31499239,5620060), GaussianInteger(5552014,28659322), GaussianInteger(660680,8921096), GaussianInteger(3459037,26739492), GaussianInteger(28825123,31283844), GaussianInteger(14423115,8019972), GaussianInteger(11896789,24786736), GaussianInteger(3030671,30809124), GaussianInteger(27408896,8373001), GaussianInteger(18883365,24606731), GaussianInteger(25911006,5006291), GaussianInteger(12954130,502075), GaussianInteger(12139386,27697690), GaussianInteger(24569256,23524970), GaussianInteger(5286017,30596546), GaussianInteger(10118517,23446156), GaussianInteger(22633676,28959584), GaussianInteger(10745806,6636327), GaussianInteger(7807322,16350706), GaussianInteger(11097391,29692385), GaussianInteger(24839267,19819633), GaussianInteger(24539265,845508), GaussianInteger(7538218,17310535), GaussianInteger(9584998,4775846), GaussianInteger(22090141,12839464), GaussianInteger(8323977,4843085), GaussianInteger(16328900,21425253), GaussianInteger(1924691,25995311), GaussianInteger(17266850,28651849), GaussianInteger(24672981,7845055), GaussianInteger(25858820,8282102), GaussianInteger(27638383,8736072), GaussianInteger(5545714,26855830), GaussianInteger(24121222,30526487), GaussianInteger(13861993,29106990), GaussianInteger(28311258,12486311), GaussianInteger(751249,9666838), GaussianInteger(30923535,31631147), GaussianInteger(5510269,14397126), GaussianInteger(33525262,30260907), GaussianInteger(524029,19949192), GaussianInteger(2556969,19119488), GaussianInteger(13939048,6557306), GaussianInteger(10885945,11091881), GaussianInteger(23586671,31199162), GaussianInteger(2490555,1404890), GaussianInteger(15161545,8745605), GaussianInteger(6527662,33106537), GaussianInteger(31468549,17818694), GaussianInteger(2546935,18491176), GaussianInteger(7689538,6338115), GaussianInteger(15659666,32831321), GaussianInteger(12882779,7795672), GaussianInteger(3786598,32688287), GaussianInteger(31407989,33139647), GaussianInteger(8676059,30086029), GaussianInteger(6131737,2539066), GaussianInteger(16714914,25531045), GaussianInteger(22322353,3984675), GaussianInteger(20654439,4972121), GaussianInteger(456756,15144030), GaussianInteger(15152195,32145697), GaussianInteger(24977780,13378206), GaussianInteger(11107796,26997483), GaussianInteger(9234245,12206736), GaussianInteger(5380786,8880871), GaussianInteger(26458023,8242396), GaussianInteger(4367089,24507289), GaussianInteger(13338508,4703735), GaussianInteger(392788,13333718), GaussianInteger(13699097,31094914), GaussianInteger(6220127,1312369), GaussianInteger(22092205,10862237), GaussianInteger(25040741,26408790), GaussianInteger(10119159,22961710), GaussianInteger(31225162,14762539), GaussianInteger(9775324,21024207), GaussianInteger(32144103,3321960), GaussianInteger(25592996,17864258)]
invNthroots = [GaussianInteger(1,0), GaussianInteger(17864258,7971677), GaussianInteger(3321960,1420570), GaussianInteger(21024207,23789349), GaussianInteger(14762539,2339511), GaussianInteger(22961710,23445514), GaussianInteger(26408790,8523932), GaussianInteger(10862237,11472468), GaussianInteger(1312369,27344546), GaussianInteger(31094914,19865576), GaussianInteger(13333718,33171885), GaussianInteger(4703735,20226165), GaussianInteger(24507289,29197584), GaussianInteger(8242396,7106650), GaussianInteger(8880871,28183887), GaussianInteger(12206736,24330428), GaussianInteger(26997483,22456877), GaussianInteger(13378206,8586893), GaussianInteger(32145697,18412478), GaussianInteger(15144030,33107917), GaussianInteger(4972121,12910234), GaussianInteger(3984675,11242320), GaussianInteger(25531045,16849759), GaussianInteger(2539066,27432936), GaussianInteger(30086029,24888614), GaussianInteger(33139647,2156684), GaussianInteger(32688287,29778075), GaussianInteger(7795672,20681894), GaussianInteger(32831321,17905007), GaussianInteger(6338115,25875135), GaussianInteger(18491176,31017738), GaussianInteger(17818694,2096124), GaussianInteger(33106537,27037011), GaussianInteger(8745605,18403128), GaussianInteger(1404890,31074118), GaussianInteger(31199162,9978002), GaussianInteger(11091881,22678728), GaussianInteger(6557306,19625625), GaussianInteger(19119488,31007704), GaussianInteger(19949192,33040644), GaussianInteger(30260907,39411), GaussianInteger(14397126,28054404), GaussianInteger(31631147,2641138), GaussianInteger(9666838,32813424), GaussianInteger(12486311,5253415), GaussianInteger(29106990,19702680), GaussianInteger(30526487,9443451), GaussianInteger(26855830,28018959), GaussianInteger(8736072,5926290), GaussianInteger(8282102,7705853), GaussianInteger(7845055,8891692), GaussianInteger(28651849,16297823), GaussianInteger(25995311,31639982), GaussianInteger(21425253,17235773), GaussianInteger(4843085,25240696), GaussianInteger(12839464,11474532), GaussianInteger(4775846,23979675), GaussianInteger(17310535,26026455), GaussianInteger(845508,9025408), GaussianInteger(19819633,8725406), GaussianInteger(29692385,22467282), GaussianInteger(16350706,25757351), GaussianInteger(6636327,22818867), GaussianInteger(28959584,10930997), GaussianInteger(23446156,23446156), GaussianInteger(30596546,28278656), GaussianInteger(23524970,8995417), GaussianInteger(27697690,21425287), GaussianInteger(502075,20610543), GaussianInteger(5006291,7653667), GaussianInteger(24606731,14681308), GaussianInteger(8373001,6155777), GaussianInteger(30809124,30534002), GaussianInteger(24786736,21667884), GaussianInteger(8019972,19141558), GaussianInteger(31283844,4739550), GaussianInteger(26739492,30105636), GaussianInteger(8921096,32903993), GaussianInteger(28659322,28012659), GaussianInteger(5620060,2065434), GaussianInteger(31647269,18854155), GaussianInteger(33026417,27755073), GaussianInteger(17313503,28599280), GaussianInteger(20468383,31872617), GaussianInteger(24180444,9926667), GaussianInteger(19915016,29443998), GaussianInteger(30088708,2621293), GaussianInteger(4571066,14978873), GaussianInteger(4563035,28387543), GaussianInteger(27174900,22464350), GaussianInteger(30014775,17874514), GaussianInteger(32746868,1699671), GaussianInteger(5469990,8768317), GaussianInteger(15819278,14560799), GaussianInteger(23028000,17889895), GaussianInteger(6833904,9997145), GaussianInteger(22762078,23505813), GaussianInteger(24457267,14384570), GaussianInteger(10247717,5061981), GaussianInteger(4395153,12131104), GaussianInteger(3474842,28676947), GaussianInteger(18048009,9537562), GaussianInteger(26717892,13707293), GaussianInteger(9672268,10363231), GaussianInteger(5448824,7461676), GaussianInteger(2010662,26345657), GaussianInteger(9871699,31992956), GaussianInteger(4698604,4748066), GaussianInteger(10801637,23642879), GaussianInteger(23403745,31381335), GaussianInteger(15018214,23231364), GaussianInteger(23618973,5432742), GaussianInteger(28637337,109218), GaussianInteger(2561481,12298551), GaussianInteger(21523023,2047396), GaussianInteger(19649036,22157197), GaussianInteger(4451054,17568717), GaussianInteger(18281015,20800471), GaussianInteger(2413724,6053211), GaussianInteger(3306177,19868662), GaussianInteger(20601016,27516298), GaussianInteger(33007623,3140746), GaussianInteger(21035353,3823315), GaussianInteger(25041160,14015623), GaussianInteger(20987741,26614259), GaussianInteger(7368331,9669413), GaussianInteger(29178387,3023882), GaussianInteger(23993255,30712754)]
