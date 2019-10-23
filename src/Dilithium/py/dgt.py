from gaussian import GaussianInteger
from params import N, p, kth_root_of_i, inv_kth_root_of_i, g, g_inv, invkmodp

fold   = lambda a: [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
unfold = lambda a: [a[j].re for j in range(N//2)] + [a[j].imag for j in range(N//2)]

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

def dgt(x):

    # Folding the input signal into N//2 Gaussian integers
    x_folded = fold(x)
    # print(x_folded)
    # Applying the right-angle convolution
    x_folded = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]
    
    # print("x_folded: ", end="")
    # print(x_folded)
    # print("\n\n")

    k = len(x_folded)

    m = k//2
    window = 1
    while m >= 1:
        index = 0
        for j in range(m):
            a = pow(g, index, p)
            i = j
            while i < k:
                xi = x_folded[i]
                xim = x_folded[i + m]
                
                x_folded[i] = (xi + xim) % p
                x_folded[i + m] = (a * (xi - xim)) % p

                i = i + (m << 1)
            index += window
        m >>= 1
        window <<= 1

    return x_folded

def idgt(x_folded):

    k = len(x_folded)
    x = list(x_folded)

    m = 1
    while m < k:
        for j in range(m):
            a = pow(g_inv, int((j * k)/(2 * m)), p)

            i = j
            while i < k:
                xi = x[i]
                xim = x[i + m]

                x[i] = (xi + a * xim) % p
                x[i + m] = (xi - a * xim) % p

                i = i + 2*m
        m = 2*m 

    # Removing the twisting factors and scaling by k^-1
    x = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(x)]

    # Unfolding the output signal
    return unfold(x)
