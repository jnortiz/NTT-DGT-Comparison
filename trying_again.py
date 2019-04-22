import sys

N = 512
q = 12889
psi = 10302

invMod = lambda y,q:pow(y,q-2,q)

def integer_reversal(num):
	x = num
	y = 0
	while x != 0:
		y <<= 1
		y |= x & 1
		x >>= 1
	return y