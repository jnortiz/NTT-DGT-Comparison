import sys
import unittest
import time
from math import log,floor
import generate_prime as Prime

invMod = lambda y,q:pow(y,q-2,q)
N = 4

#################################################################
# Apply schoolbook polynomial multiplication and reduces by Z_q // <x^N + 1>
def mul(a, b, q):
    assert len(a) == len(b)
    N = len(a)
    c = [0]*N

    # Mul and reduce
    for i in range(N):
        for j in range(N):
            v = a[i]*b[j]*(-1)**(int((i+j)//float(N)))

            c[(i+j) % N] = (c[(i+j) % N] + v) % q

    return c

#################################################################
def generate_k(n,not_approved=[]):
    k = 1
    while (not Prime.is_prime(k*n+1)) or (k in not_approved):
        k+=1
    return k

#################################################################
def generate_primitive_root(n,q):
    not_approved = []
    r = None

    while r is None:
        k = generate_k(n,not_approved)

        for i in range(q):
            if Prime.is_prime(i):
                s = set()
                for j in range(q):
                    s.add(i**j % q)
                if s == set(range(q))-{0}:
                    r = i
        if r is None:
            not_approved.append(k)

    return r,k

#################################################################
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

#################################################################
def ntt_in_place(x,wN,q):
    N = len(x)
    bitrev_shuffle(x)
    trans_size = 2
    for trans_size in [2**i for i in range(1,int(log(N,2))+1)]:
        wb = 1
        wb_step = wN**(N//trans_size)

        for t in range(trans_size >> 1):
            for trans in range(N // trans_size):
                i = trans * trans_size + t
                j = i + (trans_size >> 1)
                a = x[i] % q
                b = x[j] * wb
                x[i] = a + b % q
                x[j] = a - b % q
            wb *= wb_step

def ntt(x,q):
    x = list(x)
    r,k = generate_primitive_root(len(x),q)  
    wN = r**k 
    ntt_in_place(x,wN,q)
    return  [y % q for y in x]

def intt_in_place(x,wN,q):
    N = len(x)
    bitrev_shuffle(x)
    trans_size = 2
    for trans_size in [2**i for i in range(1,int(log(N,2))+1)]:
        wb = 1
        # wb_step = exp(2j * pi / trans_size)
        wb_step = invMod(wN**(N//trans_size) % q,q)

        for t in range(trans_size >> 1):
            for trans in range(N // trans_size):
                i = trans * trans_size + t
                j = i + (trans_size >> 1)
                a = x[i]
                b = x[j] * wb
                x[i] = a + b % q
                x[j] = a - b % q
            wb *= wb_step

def intt(x,q):
    N = len(x)
    x = list(x)
    r,k = generate_primitive_root(len(x),q)
    wN = r**k
    intt_in_place(x,wN,q)
    return  [y*invMod(len(x),q)%q for y in x]

def ntt_mul(a, b, q):
	
    a_ntt = ntt(a,q)
    assert(intt(a_ntt, q) == a)
    b_ntt = ntt(b,q)
    assert(intt(b_ntt, q) == b)

    mult = [(x*y)%q for x, y in zip(a_ntt, b_ntt)]
    
    c = intt(mult,q)

    return c

class TestNTT(unittest.TestCase):

    def test_transformation(self):
        print("Testing iterative NTT")
        x = [x for x in range (N)]
        q = Prime.generate_large_prime(log(N,2)+1)

        start_time = time.time()
        intt(ntt(x,q),q)
        end_time = time.time()
        print("----------", end_time - start_time, "s. ----------")

        self.assertEqual(x, intt(ntt(x,q),q))

    def test_mul(self):
        print("Testing polynomial multiplication using NTT")
        
        x = [x for x in range (N)]
        y = [y for y in range (N)]

        q = Prime.generate_large_prime(log(N,2)+1)

        start_time = time.time()
        ntt_mul(x, y, q)
        end_time = time.time()
        print("----------", end_time - start_time, "s. ----------")

        print(mul(x,y,q))
        print(ntt_mul(x,y,q))

        self.assertEqual(mul(x,y,q), ntt_mul(x,y,q))

if __name__ == '__main__':
    unittest.main()