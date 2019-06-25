from dgt import *
from ntt import *

N_RUN = 10**3

def compare_ntt_dgt():
    print("Comparing DGT and NTT")
    
    status = True
    for _ in range(N_RUN):    
        a = gen_polynomial_modp(N)
        b = gen_polynomial_modp(N)
        
        if(ntt_mul(a, b) != dgt_gentlemansande_mul(a, b))
            status = False
            
    return status

print(compare_ntt_dgt())