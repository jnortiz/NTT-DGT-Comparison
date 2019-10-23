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

# Computing the powers of the k-th root of i mod p
kth_root_of_i = [ pow(kthroots[N//2], i) % p for i in range(N//2) ]
inv_kth_root_of_i = [ pow(invkthroots[N//2], i) % p for i in range(N//2) ]

# Scalar k^-1 mod p used in to scale the output signal in backward DGT
invkmodp = modinv(k, p)

# for i in range(k//2):
# 	print((pow(g,i,p)) % p, end=", ")

# print("\n-----------------------------------------\n")

# for i in range(k//2):
# 	print((pow(g_inv,i,p)) % p, end=", ")
