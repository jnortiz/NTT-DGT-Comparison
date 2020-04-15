from gaussian import GaussianInteger
from params import p, g, g_inv

bit_reverse = lambda n, width: int('{:0{width}b}'.format(n, width=width)[::-1], 2)
log_2 = lambda n: n.bit_length()

def dgt_gentleman_sande(A):
    n = len(A)

    m = n//2
    while m >= 1:
        for j in range(m):
            a = pow(g, (j * n)//(2 * m), p)
            i = j
            while i < n:
                xi = A[i]
                xim = A[i + m]
                
                A[i] = (xi + xim) % p
                A[i + m] = (a * (xi - xim)) % p

                i = i + 2*m
        m = m//2

    return A

def idgt_cooley_tukey(A):

    n = len(A)

    m = 1
    while m < n:
        for j in range(m):
            a = pow(g_inv, int((j * n)/(2 * m)), p)

            i = j
            while i < n:
                xi = A[i]
                xim = A[i + m]

                A[i] = (xi + a * xim) % p
                A[i + m] = (xi - a * xim) % p

                i = i + 2*m
        m = 2*m 

    return A

def dgt_cooley_tukey(A):

    n = len(A)

    m = 1
    Distance = n//2    
    while(m < n):        
        for k in range(m):            
            JFirst = 2 * k * Distance
            JLast = JFirst + Distance - 1
            a = pow(g, bit_reverse(k, log_2(n//2)-1))
            for j in range(JFirst, JLast+1, 1):
                Temp = a * A[j + Distance]
                A[j + Distance] = (A[j] - Temp) % p
                A[j] = (A[j] + Temp) % p        
        m = m * 2
        Distance = Distance//2

    return A

def idgt_gentleman_sande(A):

    n = len(A)

    ProblemSize = n
    Distance = 1
    while(ProblemSize > 1):        
        for JFirst in range(Distance):
            j = JFirst
            r = 0
            while(j < n):
                a = pow(g_inv, bit_reverse(r, log_2(n//2)-1))
                #print(r, end=", ")
                Temp = A[j]
                A[j] = (Temp + A[j + Distance]) % p
                A[j + Distance] = ((Temp - A[j + Distance]) * a) % p
                r += 1
                j = j + 2 * Distance

        ProblemSize = ProblemSize//2
        Distance = Distance * 2

    return A