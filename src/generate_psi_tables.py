import sys
import itertools, numbers
import random

N = 2048
q = 1129725953

# psi is a n-th primitive root of unity such that psi^{n} = -1 mod q and psi^{2n} = 1 mod q
def find_root():
	root = 0
	for j in range(q):
		if(pow(j,N,q) == (q-1)):
			root = j
			break
	return root

# Returns the multiplicative inverse of n modulo mod. The inverse x has the property that
# 0 <= x < mod and (x * n) % mod = 1. The inverse exists if and only if gcd(n, mod) = 1.
def reciprocal(n, mod):
	check_int(n)
	check_int(mod)
	if not (0 <= n < mod):
		raise ValueError()
	x, y = mod, n
	a, b = 0, 1
	while y != 0:
		a, b = b, a - x // y * b
		x, y = y, x % y
	if x == 1:
		return a % mod
	else:
		raise ValueError("Reciprocal does not exist")

# Returns floor(sqrt(n)) for the given integer n >= 0.
def sqrt(n):
	check_int(n)
	if n < 0:
		raise ValueError()
	i = 1
	while i * i <= n:
		i *= 2
	result = 0
	while i > 0:
		if (result + i)**2 <= n:
			result += i
		i //= 2
	return result

# Returns a list of unique prime factors of the given integer in
# ascending order. For example, unique_prime_factors(60) = [2, 3, 5].
def unique_prime_factors(n):
	check_int(n)
	if n < 1:
		raise ValueError()
	result = []
	i = 2
	end = sqrt(n)
	while i <= end:
		if n % i == 0:
			n //= i
			result.append(i)
			while n % i == 0:
				n //= i
			end = sqrt(n)
		i += 1
	if n > 1:
		result.append(n)
	return result

# Returns silently if the given value is an integer, otherwise raises a TypeError.
def check_int(n):
	if not isinstance(n, numbers.Integral):
		raise TypeError()

# Tests whether val is a primitive degree-th root of unity modulo mod.
# In other words, val^degree % mod = 1, and for each 1 <= k < degree, val^k % mod != 1.
def is_primitive_root(val, degree, mod):
	check_int(val)
	check_int(degree)
	check_int(mod)
	if not (0 <= val < mod):
		raise ValueError()
	if not (1 <= degree < mod):
		raise ValueError()
	pf = unique_prime_factors(degree)
	return pow(val, degree, mod) == 1 and all((pow(val, degree // p, mod) != 1) for p in pf)

def bitrev_shuffle(x):
    N = len(x)
    j = 0
    for i in range(1, N):
        b = N >> 1
        while j >= b:
            j -= b
            b >>= 1
        j += b
        if j > i:
            x[i], x[j] = x[j], x[i]
    return x

psi = find_root()
print(is_primitive_root(psi,N*2,q))
print(psi)

_str = []
for j in range(N):
	_str.append(pow(psi, j, q))

psi_inv = reciprocal(psi,q)

_str_inv = [0]*N
for j in range(N):
	_str_inv[j] = pow(psi_inv, j, q)

for i in range(N):
	assert((_str[i]*_str_inv[i])%q == 1)

print(_str)
print("\n")
print(_str_inv)

bitrev_shuffle(_str)
bitrev_shuffle(_str_inv)

print(_str)
print("\n")
print(_str_inv)