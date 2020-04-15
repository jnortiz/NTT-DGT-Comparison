from gaussian import GaussianInteger

modinv = lambda y,p:pow(y,p-2,p)

# qTESLA-p-I
p = 343576577
N = 1024
k = N//2

# k-th root of i mod p
if p == 12289: # NewHope n = 1024 and p = 12,289
	kthroots = {
    512:GaussianInteger(8458, 7974)
	}
	invkthroots = {
    512: GaussianInteger(1330, 3875)
	}

elif p == 8380417: # Dilithium n = 256 and p = 8,380,417
	kthroots = {
		128:GaussianInteger(2089168, 5159577)
	}

	invkthroots = {
		128:GaussianInteger(4799789, 4692735)
	}

elif p == 343576577: # qTESLA-p-I n = 1024 and 343,576,577
	kthroots = {
		512:GaussianInteger(214546212, 116281851)
	}

	invkthroots = {
		512:GaussianInteger(237202985, 175356116)
	}
elif p == 856145921: # qTESLA-p-III n = 2048 and 853,145,921
	kthroots = {
		1024:GaussianInteger(22345081, 298508328)
	}

	invkthroots = {
		1024:GaussianInteger(673284886, 470366539)
	}	
else:	
	raise Exception()

# Primitive root of p
PROOTS = {
    12289:11,
	8380417:10,
	343576577:3,
	856145921:3
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